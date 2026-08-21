library(tidyverse)
library(cowplot)
library(scales)
library(zoo)

# ------------------------------------------------------------
# Infection threshold (single source of truth)
# ------------------------------------------------------------
UMI_THRESHOLD <- 10

clusters <- read_csv(  "cell_barcode_sample_cluster.csv",
  show_col_types = F)

clusters <- clusters %>%  mutate(
    sample_id = stringr::str_extract(Sample, "\\d+")  ) %>%
  rename(    cell_barcode = Barcode  ) %>%
  mutate(    cluster = gsub("Cluster ", "", cluster),
   cluster = as.character(cluster)  )

# Fixed biological order
cluster_order <- c(
  "18", "20", "12", "7", "21", "5", "22", "2", "9", "13",
  "19", "17",  "10",  "14","16",
  "15","11","8", "4", "6", "1","3")

cluster_cols <- c(  "7"  = "#67000D",
  "9"  = "#67000D",
  "20" = "#67000D",   "18" = "#67000D",
  "19" = "#67000D",  "21" = "#67000D",
  "12" = "#67000D",
  "13" = "#67000D",

  "1"  = "black",
  "2"  = "black",
  "5"  = "black",
  "8"  = "black",
  "10" = "black",
  "14" = "black",
  "16" = "black",
  "22" = "black",

  "3"  = "#08519C",
  "4"  = "#08519C",
  "15" = "#08519C",
  "17" = "#08519C",  "6"  = "#08519C",
  "11" = "#08519C")

cluster_labels <- c(
  "7"  = "Cycling CD4-associated T cells I",
  "9"  = "Conventional CD4 T cells",
  "20" = "Cycling CD4-associated T cells II",
  "18" = "TCF7+ early activated T cells",
  "19" = "CD44high activated/memory-like T cells",
  "21" = "Candidate MDV-transformed T cells",
  "12" = "Cycling T cells",
  "13" = "Activated cytotoxic T cells",

  "1"  = "Thrombocytes I",
  "2"  = "Plasma cells",
  "5"  = "Antigen-presenting cells",
  "8"  = "Conventional dendritic cells",
  "10" = "Metabolically active MHC-II+ immune cells",
  "14" = "Testicular stromal/supporting cells",
  "16"  = "Thrombocytes II",
  "22" = "Activated mature B cells",

  "3"  = "Innate-like T cells",
  "4"  = "Antigen-presenting B cells",
  "15" = "TCF7+ IL7R+ resting/early T cells",
  "17" = "TCF7+ TOX+ transitional T cells",
  "6"  = "GNLY+ XCL1+ cytotoxic lymphocytes",
  "11" = "Mature B cells"
)

# ------------------------------------------------------------
# Load reads
# ------------------------------------------------------------

samples <- c(  "02","03","04","05",
  "06","07","11","12")

all_reads <- map_df(samples, function(s){

  file <- paste0(
    "mdv_filtered/MDV_",    s,
    ".reads.tsv"  )

  read_tsv(
    file,
    col_names = c(
      "CB","UB","chrom",
      "start","end","NH"    ),
    show_col_types = FALSE  ) %>%
    mutate(sample_id = s) %>%
    filter(
      CB != "NA",
      UB != "NA"    )})

all_reads <- all_reads %>%
  rename(cell_barcode = CB) %>%
  inner_join(    clusters,
    by = c("cell_barcode","sample_id") )

# ------------------------------------------------------------
# Cell-level UMI counts
# ------------------------------------------------------------

cell_stats <- all_reads %>%
  group_by(
    sample_id,
    cell_barcode,
    cluster  ) %>%  summarise(
    n_umis = n_distinct(UB),
    .groups = "drop"   )

# ------------------------------------------------------------
# Cluster summary
# A cell is only "infected" if it has >= UMI_THRESHOLD UMIs
# ------------------------------------------------------------

cluster_summary <- cell_stats %>%
  group_by(cluster) %>%
  summarise(
    total_cells = n(),
    infected_cells = sum(n_umis >= UMI_THRESHOLD),
    fraction_infected = infected_cells / total_cells,
    mean_umis = mean(
      n_umis[n_umis >= UMI_THRESHOLD],      na.rm = TRUE ),
    .groups = "drop"  ) %>%
  replace_na(    list(mean_umis = 0)  )

infected_cells <- cell_stats %>% filter(n_umis >= UMI_THRESHOLD)

all_reads <- all_reads %>%
  inner_join(
    infected_cells %>%
      select(
        sample_id,
        cell_barcode      ),
    by = c(      "sample_id",
      "cell_barcode" )  )

# ------------------------------------------------------------
# Genomic positions
# ------------------------------------------------------------

all_reads <- all_reads %>%
  mutate(
    pos = round((start + end) / 2)
  ) %>%
  filter(
    !is.na(pos)
  ) %>%
  distinct(
    cluster,
    pos,
    cell_barcode,
    UB
  )

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

window_size <- 500
step_size <- 100
k <- window_size / step_size

# ------------------------------------------------------------
# Coverage calculation
# Positions with depth below UMI_THRESHOLD are set to NA so they
# are not drawn - the coverage track only shows signal that meets
# the same infection criterion used for the cell-level metrics.
# ------------------------------------------------------------

depth <- map_df(cluster_order, function(cl){

  df <- all_reads %>%
    filter(as.character(cluster) == cl)

  if(nrow(df) == 0){    return(tibble())  }

  df %>%    mutate(
      bin = floor(pos / step_size)
    ) %>%    count(
      bin,      name = "depth"
    ) %>%    arrange(bin) %>%
    mutate(      depth = zoo::rollsum(
        depth,        k = k,        fill = NA,
        align = "center"      ),
      depth = replace_na(depth, 0),
      pos = bin * step_size,
      # apply the 10-UMI infection threshold to the coverage track itself
      depth_thresholded = ifelse(depth >= UMI_THRESHOLD, depth, NA_real_),
      log_depth = log10(depth_thresholded),
      cluster = cl    )

})

# ------------------------------------------------------------
# Annotation file
# ------------------------------------------------------------

annot_file <- if (
  file.exists("EF523390.1.masked.gtf")
) {
  "EF523390.1.masked.gtf"
} else {
  "EF523390.1.masked.gff"
}

annot <- read_tsv(
  annot_file,
  comment = "#",
  col_names = c(
    "chrom","source","feature",
    "start","end","score",
    "strand","frame","attribute"
  ),
  show_col_types = FALSE
)

genes <- annot %>%
  filter(feature == "CDS") %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  ) %>%
  mutate(
    ymin = if_else(strand == "+", 0.52, 0.18),
    ymax = if_else(strand == "+", 0.82, 0.48)
  )
str(genes) # ------------------------------------------------------

highlight_regions <- tibble(
  xmin = c(134367,143796),
  xmax = c(136386,150769),
  colour = c("red","darkgreen")
)

genome_end <- max(
  c(depth$pos, genes$end),
  na.rm = TRUE
)

xlims <- c(0, genome_end)

xbreaks <- seq(
  0,
  ceiling(genome_end / 10000) * 10000,
  by = 10000
)

# y-axis floor is now the infection threshold (log10 scale), not zero
y_floor <- log10(1) # UMI_THRESHOLD)

ymax <- max(
  depth$log_depth,
  na.rm = TRUE
)

make_plot <- function(df, cl){
  stats <- cluster_summary %>%
    filter(cluster == cl)
  label_text <- paste0(
    "Cluster ", cl,
    " | ",
    unname(cluster_labels[cl]),
    "\n",
    round(stats$fraction_infected * 100, 1),
    "% infected cells | ",
    round(stats$mean_umis, 2),    " UMIs/cell"  )

  ggplot(df, aes(pos, log_depth)) +
    geom_rect(
      data = highlight_regions,
      inherit.aes = FALSE,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = colour
      ),
      alpha = 0.20
    ) +
    geom_point(
      colour = unname(cluster_cols[cl]),
      size = 1,
      alpha = 0.6,
      na.rm = TRUE
    ) +
    annotate(
      "text",
      x = xlims[1] + 0.02 * diff(xlims),
      y = ymax * 1.05,
      label = label_text,
      colour = "black",
      hjust = 0,
      vjust = 1,
      size = 4.5,
      fontface = "bold"
    ) +
    scale_fill_identity() +
    scale_x_continuous(
      limits = xlims,
      breaks = xbreaks,
      labels = function(x) x/1000,
      expand = c(0,0)    ) +
    coord_cartesian(
      xlim = xlims,
      ylim = c(y_floor, ymax * 1.05),
      expand = F    ) +
    scale_y_continuous(
      limits = c(y_floor, ymax * 1.05),
      expand = c(0,0),
      labels = function(x) round(10^x)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )
}
# ------------------------------------------------------------

plot_list <- lapply(  cluster_order,
  function(cl){
    make_plot(
      depth %>% filter(cluster == cl), cl ) })

p_annot <- ggplot(genes) +
  geom_rect(
    aes( xmin = start,      xmax = end,
      ymin = ymin,      ymax = ymax    ),
    fill = "grey60",    colour = "black"  ) +
  scale_x_continuous(
    limits = xlims,
    breaks = xbreaks,
    labels = function(x) x/1000,
    expand = c(0,0)  ) +
  coord_cartesian(    xlim = xlims,    expand = F  ) +
  theme_minimal() +
  labs(    x = "Genomic position (Kb)",
    y = "Genes"   )

combined <- plot_grid(  plotlist = c(
    plot_list,   list(p_annot) ,  list(p_annot)  ),  ncol = 2,
  align = "hv",     # align both horizontally and vertically
  axis  = "tblr"    # match top/bottom/left/right margins across all panels
)
ggsave(  "plot_cluster_infection_sorted.png",  combined,
  width = 14,  height = 14,  dpi = 300 )

ggsave(  "plot_cluster_infection_sorted.pdf",
  combined,  width = 14,  height = 14 )