library(tidyverse)
library(viridis)
library(ggrepel)

# ------------------------------------------------------------
# 1. Load cluster mapping
# ------------------------------------------------------------

clusters <- read_csv(
  "cell_barcode_sample_cluster.csv",
  show_col_types = FALSE
)

clusters <- clusters %>%
  mutate(sample_id = stringr::str_extract(Sample, "\\d+")) %>%
  rename(
    cell_barcode = Barcode,
    cluster = cluster
  )

# ------------------------------------------------------------
# 2. Viral reads
# ------------------------------------------------------------

samples <- c("02","03","04","05","06","07","11","12")

viral <- map_df(samples, function(s) {

  file <- paste0("mdv_filtered/MDV_", s, ".reads.tsv")

  if (!file.exists(file))
    stop(paste("Missing file:", file))

  read_tsv(
    file,
    col_names = c(
      "CB","UB","chrom",
      "start","end","NH"
    ),
    show_col_types = FALSE
  ) %>%
    mutate(sample_id = s) %>%
    filter(CB != "NA", UB != "NA") %>%
    group_by(
      sample_id,
      cell_barcode = CB
    ) %>%
    summarise(
      n_unique_umis = n_distinct(UB),
      n_reads = n(),
      .groups = "drop"
    )
})

# ------------------------------------------------------------
# 3. Total cells per cluster
# ------------------------------------------------------------

total_cells <- clusters %>%
  filter(sample_id %in% samples) %>%
  group_by(sample_id, cluster) %>%
  summarise(
    total_cells = n(),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 4. Infected cells
# ------------------------------------------------------------

infected_cells <- viral %>%
  inner_join(
    clusters,
    by = c("cell_barcode", "sample_id") ) %>%
  filter(n_unique_umis > 10) %>%
  group_by(sample_id, cluster) %>%
  summarise(
    infected_cells = n(),
    total_umis = sum(n_unique_umis),
    total_reads = sum(n_reads),
    .groups = "drop"  )

# ------------------------------------------------------------
# 5. Final summary table
# ------------------------------------------------------------

final <- total_cells %>%
  left_join(
    infected_cells,
    by = c("sample_id", "cluster")  ) %>%
  mutate(
    infected_cells = replace_na(infected_cells, 0),
    total_umis = replace_na(total_umis, 0),
    total_reads = replace_na(total_reads, 0),

    fraction_infected =
      infected_cells / total_cells,

    umis_per_infected_cell =
      if_else(
        infected_cells > 0,
        total_umis / infected_cells,
        NA_real_
      )
  )

write_tsv(
  final,
  "MDV_cluster_sample_full_summary.tsv"
)

# ------------------------------------------------------------
# Sample labels
# ------------------------------------------------------------

sample_map <- tibble(
  sample_id = c(
    "02","03","04","05",
    "06","12"#07","11","12"
  ),
  label = c(
    "Spleen-1",
    "Kidney-1",
    "Testes-1",
    "Kidney-2",
    "Testes-2",
 #   "Control-1",
 #  "Control-2",
    "Spleen-2"
  )
)

label_order <- c(
#  "Control-1",
 # "Control-2",
  "Kidney-1",
  "Kidney-2",
  "Spleen-1",
  "Spleen-2",
  "Testes-1",
  "Testes-2"
)

# ------------------------------------------------------------
# Cluster order
# ------------------------------------------------------------

cluster_order <- c(
  "7","9","20","18","19","21","12","13",
  "1","2","5","8","10","14","16", "22",
  "3","4","15","17","6","11"
)

# ------------------------------------------------------------
# Cluster colours
# ------------------------------------------------------------

cluster_cols <- c(
  "7"  = "#67000D",
  "9"  = "#8B0000",
  "20" = "#A50F15",
  "18" = "#CB181D",
  "19" = "#EF3B2C",
  "21" = "#FB6A4A",
  "12" = "#FC9272",
  "13" = "#FCAE91",

  "1"  = "#525252",
  "2"  = "#525252",
  "5"  = "#737373",
  "8"  = "#737373",
  "10" = "#737373",
  "14" = "#737373",
  "16" = "#737373",
  "22" = "#737373",

  "3"  = "#6BAED6",
  "4"  = "#6BAED6",
  "15" = "#6BAED6",
  "17" = "#6BAED6",
  "6"  = "#6BAED6",
  "11" = "#6BAED6"
)

# ------------------------------------------------------------
# Cluster legend labels
# ------------------------------------------------------------

cluster_labels <- c(
  "7"  = "7  Cycling CD4-associated T cells I",
  "9"  = "9  Conventional CD4 T cells",
  "20" = "20 Cycling CD4-associated T cells II",
  "18" = "18 TCF7+ early activated T cells",
  "19" = "19 CD44high activated/memory-like T cells",
  "21" = "21 Candidate MDV-transformed T cells",
  "12" = "12 Cycling T cells",
  "13" = "13 Activated cytotoxic T cells",

  "1"  = "1  Thrombocytes I",
  "2"  = "2  Plasma cells",
  "5"  = "5  Antigen-presenting cells",
  "8"  = "8  Conventional dendritic cells",
  "10" = "10 Metabolically active MHC-II+ immune cells",
  "14" = "14 Testicular stromal/supporting cells",
    "16" = "Thrombocytes II",
  "22" = "22 Activated mature B cells",

  "3"  = "3  Innate-like T cells",
  "4"  = "4  Antigen-presenting B cells",
  "15" = "15 TCF7+ IL7R+ resting/early T cells",
  "17" = "17 TCF7+ TOX+ transitional T cells",
  "6"  = "6  GNLY+ XCL1+ cytotoxic lymphocytes",
  "11" = "11 Mature B cells"
)

# ------------------------------------------------------------
# Plotting dataframe
# ------------------------------------------------------------

df <- final %>%
  mutate(sample_id = as.character(sample_id)) %>%
  left_join(sample_map, by = "sample_id") %>%
  mutate(
    sample = factor(
      label,
      levels = label_order
    ),
    cluster = factor(
      cluster,
      levels = cluster_order
    )
  )

df_labels <- df %>%
  filter(
    !(sample %in%
        c("Control-1", "Control-2"))
  )

df2 <- subset(df, sample!="Control-1" & sample!="Control-2")
str(df2)
df2$sample <- factor(df2$sample, levels = label_order)
# replace NA in umis_per_infected_cell with 0 for plotting purposes
df2$umis_per_infected_cell[is.na(df2$umis_per_infected_cell)] <- 0

p <- ggplot(
  df2,
  aes(
    x = fraction_infected,
    y = umis_per_infected_cell,
    colour = cluster,
    size = total_cells
  )
) +
  geom_point(
    data = df2 %>%
      filter(
        !is.na(umis_per_infected_cell),
        fraction_infected >= -0.02,
        fraction_infected <= 1.1,
        umis_per_infected_cell >= 4,
        umis_per_infected_cell <= 65
      ),
    alpha = 0.7
  ) +
  geom_text_repel(
    data = df2 %>%
      filter(
        !is.na(umis_per_infected_cell),
        umis_per_infected_cell >=1    ),
    aes(label = cluster),
    size = 7,
    fontface = "bold",
    segment.color = "grey40",
    segment.size = 0.5,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~sample,
    ncol = 3,
    drop = FALSE
  ) +
  scale_colour_manual(
    values = cluster_cols,
    breaks = cluster_order,
    labels = cluster_labels,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(-0.02, 1.1),
    expand = c(0, 0),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
  ) +
  scale_y_log10(
    limits = c(9.9, 51),
    breaks = c( 10, 20, 50),
    labels = c( "10", "20", "50")
  ) +
  labs(
    x = "Fraction of infected cells",
    y = "UMIs per infected cell",
    colour = "Cluster annotation",
    size = "Total cells"
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 7,
        alpha = 1
      ),
      ncol = 1
    )
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 21, face = "bold"),
    legend.text = element_text(size = 19),
    legend.key.height = unit(0.6, "cm"),
    strip.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )

ggsave(  "fraction_vs_UMIs_by_sample_LABELLED.pdf",
  p,  width = 16,  height = 8 )
ggsave(  "fraction_vs_UMIs_by_sample_LABELLED.png",
  p,  width = 18,  height = 8,  dpi = 300 )
