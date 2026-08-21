library(tidyverse)
library(viridis)
library(zoo)

if (!dir.exists("OUTPUT4")) dir.create("OUTPUT4")

cat("[checkpoint 1] wd =", getwd(), "\n")
cat("[checkpoint 1] input files exist:",
    file.exists("all_reads.tsv"), file.exists("cell_barcode_sample_cluster.csv"),
    file.exists("EF523390.1.masked.gb"), "\n")

# 1. Load Data
all_reads <- read_tsv("all_reads.tsv", show_col_types = FALSE) %>%
  mutate(sample_num = as.numeric(stringr::str_extract(Sample, "\\d+")))
cat("[checkpoint 2] all_reads rows:", nrow(all_reads), "\n")

clusters <- read_csv("cell_barcode_sample_cluster.csv", show_col_types = FALSE) %>%
  mutate(
    cluster = sprintf("%02d", as.numeric(cluster)),
    sample_num = as.numeric(stringr::str_extract(Sample, "\\d+"))
  )
cat("[checkpoint 3] clusters rows:", nrow(clusters), "unique clusters:",
    paste(sort(unique(clusters$cluster)), collapse=","), "\n")

# 1b. Cluster label mapping (cluster ID -> descriptive cell-type name)
cluster_labels <- c(
  "7" = "Cycling CD4-associated T cells I",
  "9" = "Conventional CD4 T cells",
  "20" = "Cycling CD4-associated T cells II",
  "18" = "TCF7+ early activated T cells",
  "19" = "CD44high activated/memory-like T cells",
  "21" = "Candidate MDV-transformed T cells",
  "12" = "Cycling T cells",
  "13" = "Activated cytotoxic T cells",
  "1" = "Thrombocytes I",
  "16" = "Thrombocytes II",
  "2" = "Plasma cells",
  "5" = "Antigen-presenting cells",
  "8" = "Conventional dendritic cells",
  "10" = "Metabolically active MHC-II+ immune cells",
  "14" = "Testicular stromal/supporting cells",
  "22" = "Activated mature B cells",
  "3" = "Innate-like T cells",
  "4" = "Antigen-presenting B cells",
  "15" = "TCF7+ IL7R+ resting/early T cells",
  "17" = "TCF7+ TOX+ transitional T cells",
  "6" = "GNLY+ XCL1+ cytotoxic lymphocytes",
  "11" = "Mature B cells"
)

# The `clusters` table zero-pads cluster IDs (e.g. "07", "09", "20"), but the
# mapping above uses unpadded keys, so build a padded version that will
# actually match the `cluster` column used for plotting.
cluster_labels_padded <- setNames(cluster_labels, sprintf("%02d", as.numeric(names(cluster_labels))))

# 2. Parse GenBank
parse_genbank_cds <- function(gb_file) {
  lines <- readLines(gb_file, warn = FALSE)
  feat_start <- grep("^FEATURES", lines)
  lines <- lines[feat_start[1]:length(lines)]
  feat_idx <- grep("^ {5}\\w+", lines)

  df_list <- lapply(seq_along(feat_idx), function(i) {
    start_idx <- feat_idx[i]
    end_idx <- if (i < length(feat_idx)) feat_idx[i+1] - 1 else length(lines)
    chunk <- lines[start_idx:end_idx]
    chunk_str <- paste(chunk, collapse = " ")
    if (!grepl("^ {5}CDS", chunk[1])) return(NULL)
    nums <- as.numeric(unlist(stringr::str_extract_all(chunk[1], "\\d+")))
    locus_tag <- stringr::str_match(chunk_str, '/locus_tag="([^"]+)"')[, 2]
    note <- stringr::str_match(chunk_str, '/note="([^"]+)"')[, 2]
    product <- stringr::str_match(chunk_str, '/product="([^"]+)"')[, 2]

    best_name <- coalesce(note, product, locus_tag, "CDS")
    raw_label <- paste0(best_name, " [", coalesce(locus_tag, "NA"), "]")
    raw_label <- str_squish(gsub("\t", " ", raw_label))

    tibble(start = min(nums), end = max(nums), raw_name = raw_label)
  })

  bind_rows(df_list) %>%
    arrange(start) %>%
    group_by(raw_name) %>%
    mutate(name = if(n() > 1) paste0(raw_name, "_", row_number()) else raw_name) %>%
    ungroup()
}

genes_gb <- parse_genbank_cds("EF523390.1.masked.gb")
genomic_order <- genes_gb$name
cat("[checkpoint 4] genes_gb rows:", nrow(genes_gb), "\n")
cat("[checkpoint 4] genomic_order:", paste(genomic_order, collapse=" | "), "\n")

# 3. Sliding Window Calculations
step_size <- 200
k_window <- 5
genomic_max <- max(genes_gb$end)

get_windowed_data <- function(df, group_col) {
  df %>%
    mutate(bin = floor(pos / step_size)) %>%
    group_by(!!sym(group_col), bin) %>%
    summarise(raw_umis = n_distinct(UB), .groups = "drop") %>%
    group_by(!!sym(group_col)) %>%
    complete(bin = 0:floor(genomic_max/step_size), fill = list(raw_umis = 0)) %>%
    mutate(
      smooth_depth = zoo::rollsum(raw_umis, k = k_window, fill = 0, align = "center"),
      pos = bin * step_size
    )
}

# Mapping Table
s_map <- tibble(sample_num = c(3, 5, 2, 12, 4, 6),
                label = c("Kidney-1", "Kidney-2", "Spleen-1", "Spleen-2", "Testes-1", "Testes-2"))

# Process
cluster_win <- get_windowed_data(all_reads %>% inner_join(clusters, by=c("Barcode","sample_num")), "cluster")
sample_win  <- get_windowed_data(all_reads %>% filter(sample_num %in% s_map$sample_num), "sample_num")
cat("[checkpoint 5] cluster_win rows:", nrow(cluster_win), " sample_win rows:", nrow(sample_win), "\n")

# Map function
map_to_genes <- function(win_df, group_col, is_sample = FALSE) {
  res <- map_df(genomic_order, function(g_name) {
    g_row <- genes_gb %>% filter(name == g_name)
    win_df %>%
      filter(pos >= g_row$start, pos <= g_row$end) %>%
      group_by(!!sym(group_col)) %>%
      summarise(peak_depth = max(smooth_depth, na.rm = TRUE), .groups = "drop") %>%
      mutate(gene_name = factor(g_name, levels = genomic_order),
             log2_val = ifelse(peak_depth > 0, log2(peak_depth), 0))
  })

  if (is_sample) {
    res <- res %>% rename(sample_num = group_col) %>%
           inner_join(s_map, by = "sample_num") %>%
           mutate(Sample = factor(label, levels = s_map$label)) %>%
           select(-sample_num, -label)
  }
  res
}

cluster_map <- map_to_genes(cluster_win, "cluster")
sample_map  <- map_to_genes(sample_win, "sample_num", is_sample = TRUE)
cat("[checkpoint 6] cluster_map rows:", nrow(cluster_map), " sample_map rows:", nrow(sample_map), "\n")

# Export Tables
write_tsv(cluster_map, "OUTPUT4/table_cluster_heatmap.tsv")
write_tsv(sample_map, "OUTPUT4/table_sample_heatmap.tsv")
cat("[checkpoint 7] tables written:",
    file.exists("OUTPUT4/table_cluster_heatmap.tsv"),
    file.exists("OUTPUT4/table_sample_heatmap.tsv"), "\n")
cluster_order_padded <- sprintf("%02d", as.numeric(names(cluster_labels)))

make_heatmap <- function(df, x_col, title, x_labels = NULL, x_order = NULL) {
  p <- ggplot(df, aes(x = !!sym(x_col), y = gene_name, fill = log2_val)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_viridis_c(name = "Log2(UMI)", option = "magma") +
    scale_y_discrete(limits = genomic_order) +
    theme_minimal() +
    theme( axis.text.x = element_text(size = 14, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 10, face = "bold"),
      panel.grid = element_blank(),
      legend.key.height = unit(10, "line"),
      legend.key.width = unit(3, "line"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18, face = "bold")
    ) + labs(x = "", y = "")

  if (!is.null(x_labels)) {
    p <- p + scale_x_discrete(limits = x_order, labels = function(x) paste0(x, " - ", x_labels[x]))
  } else if (!is.null(x_order)) {
    p <- p + scale_x_discrete(limits = x_order)
  }

  p
}
make_heatmap <- function(df, x_col, title, x_labels = NULL, x_order = NULL) {
  p <- ggplot(df, aes(x = !!sym(x_col), y = gene_name, fill = log2_val)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_viridis_c(name = "Log2(UMI)", option = "magma") +
    scale_y_discrete(limits = genomic_order) +
    theme_minimal() +
    theme( axis.text.x = element_text(size = 14, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 10, face = "bold"),
      panel.grid = element_blank(),
      legend.key.height = unit(10, "line"),
      legend.key.width = unit(3, "line"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18, face = "bold")
    ) + labs(x = "", y = "")

  if (!is.null(x_labels)) {
    p <- p + scale_x_discrete(limits = x_order, labels = function(x) paste0(x, " - ", x_labels[x]))
  } else if (!is.null(x_order)) {
    p <- p + scale_x_discrete(limits = x_order)  }
  p }
p_cluster <- make_heatmap(cluster_map, "cluster", "By Cluster",
                           x_labels = cluster_labels_padded,
                           x_order  = cluster_order_padded)

cat("[checkpoint 8] building cluster heatmap plot...\n")
cat("[checkpoint 9] saving cluster heatmap...\n")
ggsave("OUTPUT4/heatmap_by_cluster.png", p_cluster, width = 12, height = 20)
cat("[checkpoint 10] cluster PNG exists:", file.exists("OUTPUT4/heatmap_by_cluster.png"), "\n")
cat("[checkpoint 11] building sample heatmap plot...\n")
p_sample <- make_heatmap(sample_map, "Sample", "By Sample")
cat("[checkpoint 12] saving sample heatmap...\n")
ggsave("OUTPUT4/heatmap_by_sample.png", p_sample, width = 9, height = 18)
cat("[checkpoint 13] sample PNG exists:", file.exists("OUTPUT4/heatmap_by_sample.png"), "\n")
cat("[checkpoint 14] DONE - final OUTPUT4 listing:\n")
print(list.files("OUTPUT4", full.names = TRUE))