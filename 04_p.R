library(tidyverse)

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------
df <- read_tsv("MDV_cluster_sample_full_summary.tsv",
               show_col_types = FALSE)

# make cluster a factor (ensures consistent colours)
df <- df %>%
  mutate(
    cluster = factor(cluster),
    sample_id = factor(sample_id)
  )

# ------------------------------------------------------------
# 2. Plot: fraction infected vs UMIs per cell
# ------------------------------------------------------------
p <- ggplot(df, aes(
  x = fraction_infected,
  y = umis_per_infected_cell,
  colour = cluster, size=total_cells,
)) +

  geom_point(size = 3, alpha = 0.9) +

  # optional: label extreme clusters
  geom_text(
    data = df %>% filter(umis_per_infected_cell > 30 | fraction_infected > 0.9),
    aes(label = cluster),
    size = 3,
    vjust = -0.6,
    show.legend = FALSE
  ) +

  facet_wrap(~ sample_id, ncol = 4) +

  scale_colour_manual(
    values = scales::hue_pal()(length(unique(df$cluster)))
  ) +

  labs(
    x = "Fraction of infected cells",
    y = "UMIs per infected cell",
    colour = "Cluster",
    title = "MDV infection: prevalence vs transcriptional intensity"
  ) +

  theme_minimal(base_size = 14) +

  theme(
    panel.grid.minor = element_blank(),

    # BIG LEGEND (as requested)
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.8, "cm"),

    # nicer layout
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)  )

ggsave("fraction_vs_UMIs_by_sample.pdf", p, width = 12, height = 8)
ggsave("fraction_vs_UMIs_by_sample.png", p, width = 12, height = 8, dpi = 300)
