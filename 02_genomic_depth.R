library(tidyverse)
library(cowplot)
library(scales)

samples <- c("01","02","03","04","05","06","07","08")

depth_list <- lapply(samples, function(s) {
  file <- paste0("MDV_", s, ".depth")
  
  read_tsv(file,
    col_names = c("chrom","pos","depth"),
    col_types = cols(.default = col_character())
  ) %>%
    mutate(
      pos = as.numeric(pos),
      depth = as.numeric(depth),
      sample = s
    ) %>%
    filter(!is.na(pos), !is.na(depth), depth > 0) %>%
    mutate(log_depth = log10(depth))
})

depth <- bind_rows(depth_list)

# ------------------------------------------------------------
# 2. Annotation (UNCHANGED)
# ------------------------------------------------------------
annot_file <- if (file.exists("MDV.gtf")) {
  "MDV.gtf"
} else if (file.exists("MDV.gff")) {
  "MDV.gff"
}

annot <- read_tsv(
  annot_file,
  comment = "#",
  col_names = c("chrom","source","feature","start","end",
                "score","strand","frame","attribute"),
  col_types = cols(.default = col_character())
)

genes <- annot %>%
  filter(feature == "gene") %>%
  mutate(
    start = as.numeric(start),
    end   = as.numeric(end)
  ) %>%
  filter(!is.na(start), !is.na(end))

if (nrow(genes) == 0) {
  genes <- annot %>%
    filter(feature %in% c("CDS","exon")) %>%
    mutate(
      start = as.numeric(start),
      end   = as.numeric(end)
    ) %>%
    filter(!is.na(start), !is.na(end))
}

genes <- genes %>%
  mutate(
    ymin = if_else(strand == "+", 0.52, 0.18),
    ymax = if_else(strand == "+", 0.82, 0.48)
  )

# ------------------------------------------------------------
# 3. Highlight regions (UNCHANGED)
# ------------------------------------------------------------
highlight_regions <- tibble(
  xmin = c(134367, 137565, 143796),
  xmax = c(136386, 139376, 150769),
  label = c("meq", "CxC chemokine", "ICP4"),
  colour = c("red", "blue", "darkgreen")
)

# axes
genome_end <- max(c(depth$pos, genes$end, 178000), na.rm = TRUE)
xlims <- c(0, genome_end)
xbreaks <- seq(0, ceiling(genome_end/10000)*10000, by = 10000)
ymax <- max(depth$log_depth, na.rm = TRUE)

# ------------------------------------------------------------
# 4. MULTI-SAMPLE DEPTH PLOT
# ------------------------------------------------------------
p_depth <- ggplot(depth, aes(x = pos, y = log_depth, colour = sample)) +

  geom_rect(
    data = highlight_regions,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax,
        ymin = -Inf, ymax = Inf,
        fill = colour),
    alpha = 0.2
  ) +

  geom_line(linewidth = 0.35) +

  geom_text(
    data = highlight_regions,
    inherit.aes = FALSE,
    aes(x = (xmin + xmax)/2,
        y = c(ymax - 1, ymax - 0.5, ymax - 1),
        label = label),
    size = 4,
    fontface = "bold",
    colour = "black"
  ) +

  scale_fill_identity() +

  scale_colour_manual(values = scales::hue_pal()(length(samples))) +

  scale_x_continuous(
    limits = xlims,
    breaks = xbreaks,
    labels = function(x) x / 1000,
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    name = "Read depth (log10)",
    expand = expansion(mult = c(0, 0.1))
  ) +

  labs(
    x = NULL,
    colour = "Sample"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 0, 5.5)
  )

# ------------------------------------------------------------
# 5. Annotation plot (UNCHANGED)
# ------------------------------------------------------------
p_annot <- ggplot(genes) +

  geom_rect(
    data = highlight_regions,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax,
        ymin = -Inf, ymax = Inf,
        fill = colour),
    alpha = 0.2
  ) +

  geom_rect(
    aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
    fill = "grey60",
    colour = "black",
    linewidth = 0.5
  ) +

  scale_fill_identity() +

  scale_x_continuous(
    limits = xlims,
    breaks = xbreaks,
    labels = function(x) x / 1000,
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    limits = c(0,1),
    breaks = c(0.25,0.75),
    labels = c("",""),
    expand = c(0,0)
  ) +

  labs(
    x = "Genomic position (Kb)",
    y = "Genes"
  ) +

  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(0, 5.5, 5.5, 5.5)
  )

# ------------------------------------------------------------
# 6. Combine
# ------------------------------------------------------------
combined <- plot_grid(
  p_depth,
  p_annot,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(2,1)
)

# ------------------------------------------------------------
# 7. Save
# ------------------------------------------------------------
ggsave("plot_combined_8samples.pdf", combined, width = 12, height = 5)
ggsave("plot_combined_8samples.png", combined, width = 12, height = 5, dpi = 300)
