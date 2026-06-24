library(tidyverse)

# ------------------------------------------------------------
# 1. Load cluster mapping
# ------------------------------------------------------------
clusters <- read_csv("cell_barcode_sample_cluster.csv",
                     show_col_types = FALSE)

clusters <- clusters %>%
  mutate(
    sample_id = stringr::str_extract(Sample, "\\d+")
  ) %>%
  rename(
    cell_barcode = Barcode,
    cluster = cluster
  )

# ------------------------------------------------------------
# 2. Load viral per-cell data (all samples)
# ------------------------------------------------------------
samples <- c("02","03","04","05","06","07","11","12")
samples <- c("02")

viral <- map_df(samples, function(s) {

  file <- paste0("MDV_", s, "_cell_barcode_viral_reads.tsv")

  read_tsv(file, show_col_types = FALSE) %>%
    mutate(sample_id = s)

})

# ------------------------------------------------------------
# 3. TOTAL CELLS PER CLUSTER (all cells)
# ------------------------------------------------------------
total_cells <- clusters %>%
  filter(sample_id %in% samples) %>%
  group_by(sample_id, cluster) %>%
  summarise(
    total_cells = n(),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 4. INFECTED CELLS (cells with viral reads)
# ------------------------------------------------------------
infected_cells <- viral %>%
  inner_join(clusters, by = c("cell_barcode", "sample_id")) %>%
  group_by(sample_id, cluster) %>%
  summarise(
    infected_cells = n(),
    total_umis = sum(n_unique_umis, na.rm = TRUE),
    total_reads = sum(n_reads, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 5. COMBINE FULL TABLE
# ------------------------------------------------------------
final <- total_cells %>%
  left_join(infected_cells,
            by = c("sample_id", "cluster")) %>%
  mutate(
    infected_cells = replace_na(infected_cells, 0),
    total_umis = replace_na(total_umis, 0),
    total_reads = replace_na(total_reads, 0),
    fraction_infected = infected_cells / total_cells,
    umis_per_infected_cell = if_else(infected_cells > 0,
                                     total_umis / infected_cells,
                                     0)
  ) %>%
  arrange(sample_id, cluster)

# ------------------------------------------------------------
# 6. SAVE OUTPUT
# ------------------------------------------------------------
write_tsv(final, "MDV_cluster_sample_full_summary.tsv")
