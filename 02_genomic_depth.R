library(tidyverse)
library(cowplot)
library(scales)
library(zoo)

# ------------------------------------------------------------
# 0. Sample mapping
# ------------------------------------------------------------
sample_map <- tibble(
  sample = c("07","11","03","05","02","12","04","06"),
  label = c("Control-1","Control-2",
            "Kidney-1","Kidney-2",
            "Spleen-1","Spleen-2",
            "Testes-1","Testes-2") )

plot_order <- sample_map$sample

# ------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------
window_size <- 1000
step_size   <- 200
k <- window_size / step_size

# ------------------------------------------------------------
# 1. UMI SLIDING WINDOW (FROM reads.tsv, WITH ≥2 UMI FILTER)
# ------------------------------------------------------------
depth <- map_df(plot_order, function(s) {

  file <- paste0("mdv_filtered/MDV_", s, ".reads.tsv")

  cat("Reading:", file, "\n")

  if (!file.exists(file)) {
    stop(paste("Missing file:", file))
  }

  df <- read_tsv(
    file,
    col_names = c("CB","UB","chrom","start","end","NH"),
    col_types = cols(
      CB = col_character(),
      UB = col_character(),
      chrom = col_character(),
      start = col_double(),
      end = col_double(),
      NH = col_double()
    )
  )

  if (nrow(df) == 0) {
    warning(paste("EMPTY FILE:", file))
    return(tibble())
  }

  # ------------------------------------------------------------
  # ✅ STEP 1: Clean + compute UMIs per cell
  # ------------------------------------------------------------
  df <- df %>%
    filter(CB != "NA", UB != "NA")

  cell_stats <- df %>%
    group_by(CB) %>%
    summarise(n_umis = n_distinct(UB), .groups = "drop")

  # ------------------------------------------------------------
  # ✅ STEP 2: KEEP only cells with ≥10 UMIs
  # ------------------------------------------------------------
  valid_cells <- cell_stats %>%
    filter(n_umis >= 10)
  df <- df %>%    inner_join(valid_cells, by = "CB")
  # safety check
  if (nrow(df) == 0) {
    warning(paste("All cells filtered out in:", file))
    return(tibble())
  }

  # ------------------------------------------------------------
  # ✅ STEP 3: midpoint + collapse
  # ------------------------------------------------------------
  df <- df %>%
    mutate(pos = round((start + end) / 2)) %>%
    filter(!is.na(pos)) %>%    distinct(pos, CB, UB)

  df <- df %>%    mutate(bin = floor(pos / step_size))
  counts <- df %>%     count(bin, name = "depth")

  # ------------------------------------------------------------
  # ✅ STEP 5: smoothing
  # ------------------------------------------------------------
  counts <- counts %>%
    arrange(bin) %>%
    mutate(
      depth = zoo::rollsum(depth, k = k, fill = NA, align = "center"),
      depth = replace_na(depth, 0),
      pos   = bin * step_size,
      log_depth = log10(depth + 1),
      sample = s,
      label = sample_map$label[match(s, sample_map$sample)]    )
  counts })

# ------------------------------------------------------------
# 2. Annotation
# ------------------------------------------------------------
annot_file <- if (file.exists("MDV.gtf")) "MDV.gtf" else "MDV.gff"

annot <- read_tsv(
  annot_file,
  comment = "#",
  col_names = c("chrom","source","feature","start","end",
                "score","strand","frame","attribute"),
  col_types = cols(.default = col_character())
)

genes <- annot %>%
  filter(feature == "gene") %>%
  mutate(start = as.numeric(start),
         end   = as.numeric(end)) %>%
  mutate(
    ymin = if_else(strand == "+", 0.52, 0.18),
    ymax = if_else(strand == "+", 0.82, 0.48)  )

# ------------------------------------------------------------
# 3. Highlight regions
# ------------------------------------------------------------
highlight_regions <- tibble(
  xmin = c(134367,137565,137804,138038,143796),
  xmax = c(136386,138428,138939,138242,150769),
  label = c("meq",
            "CxC chemokine","CxC chemokine","CxC chemokine",
            "ICP4"),
  colour = c("red","blue",
             "blue","blue","darkgreen")
)

label_regions <- highlight_regions %>%
  filter(label %in% c("meq","ICP4"))

# ------------------------------------------------------------
# 4. Axes
# ------------------------------------------------------------
genome_end <- max(c(depth$pos, genes$end), na.rm = TRUE)

xlims <- c(0, genome_end)
xbreaks <- seq(0, ceiling(genome_end / 10000) * 10000, by = 10000)

ymax <- max(depth$log_depth, na.rm = TRUE)

# ------------------------------------------------------------
# 5. Plot function
# ------------------------------------------------------------
make_plot <- function(df, label_text) {

  ggplot(df, aes(pos, log_depth)) +

    geom_rect(
      data = highlight_regions,
      inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax,
          ymin = -Inf, ymax = Inf,
          fill = colour),
      alpha = 0.2 ) +

    geom_point(size =1, alpha = 0.6) +

    geom_text(
      data = label_regions,
      inherit.aes = FALSE,
      aes(x = (xmin + xmax)/2,
          y = ymax - 0.5,
          label = label),
      size = 4,
      fontface = "bold"
    ) +

    annotate("text",
      x = xlims[1] + 0.45 * xlims[2],
      y = ymax - 0.5,
      label = label_text,
      size = 5,
      fontface = "bold"
    ) +

    scale_fill_identity() +

    scale_x_continuous(
      limits = xlims,
      breaks = xbreaks,
      labels = function(x) x / 1000,
      expand = c(0,0)
    ) +

    coord_cartesian(xlim = xlims, expand = FALSE) +

    scale_y_continuous(
      limits = c(1, ymax),
      expand = c(0,0),
      breaks = pretty_breaks(),
      labels = function(x) round(10^x)
    ) +

    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

# ------------------------------------------------------------
# 6. Build plots
# ------------------------------------------------------------
plot_list <- list()

for (i in seq_along(plot_order)) {

  s <- plot_order[i]
  lab <- sample_map$label[sample_map$sample == s]

  df <- depth %>% filter(sample == s)

  plot_list[[i]] <- make_plot(df, lab)
}

# ------------------------------------------------------------
# 7. Annotation panel
# ------------------------------------------------------------
p_annot <- ggplot(genes) +
  geom_rect(
    aes(xmin = start, xmax = end,
        ymin = ymin, ymax = ymax),
    fill = "grey60", colour = "black"
  ) +
  scale_x_continuous(
    limits = xlims,
    breaks = xbreaks,
    labels = function(x) x / 1000,
    expand = c(0,0)
  ) +
  coord_cartesian(xlim = xlims, expand = FALSE) +
  labs(x = "Genomic position (Kb)", y = "Genes") +
  theme_minimal()

# ------------------------------------------------------------
# 8. Combine
# ------------------------------------------------------------
combined <- plot_grid(
  plotlist = c(plot_list, list(p_annot)),
  ncol = 1,
  align = "v",
  axis = "lr" )
ggsave("plot_grouped_UMI_sliding_filtered.png",
       combined, width = 9, height = 11)
ggsave("plot_grouped_UMI_sliding_filtered.pdf",
       combined, width = 9, height = 11)