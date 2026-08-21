library(tidyverse)
library(cowplot)
library(scales)
library(zoo)
library(ggplot2)
library(viridis)
library(ggrepel)

dir.create("OUTPUT2", showWarnings = FALSE)

# ------------------------------------------------------------
# Load clusters
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
  ) %>%
  mutate(cluster = gsub("Cluster ", "", cluster))

# ------------------------------------------------------------
# Load reads (EXCLUDE controls 07, 11)
# ------------------------------------------------------------
samples <- c("02","03","04","05","06","12")

all_reads <- map_df(samples, function(s) {

  file <- paste0("mdv_filtered/MDV_", s, ".reads.tsv")

  read_tsv(
    file,
    col_names = c("CB","UB","chrom","start","end","NH"),
    show_col_types = FALSE
  ) %>%
    mutate(sample_id = s) %>%
    filter(CB != "NA", UB != "NA")
})

# ------------------------------------------------------------
# Attach clusters
# ------------------------------------------------------------
all_reads <- all_reads %>%
  rename(cell_barcode = CB) %>%
  inner_join(clusters, by = c("cell_barcode","sample_id"))

# ------------------------------------------------------------
# Compute UMIs per cell
# ------------------------------------------------------------
cell_stats <- all_reads %>%
  group_by(sample_id, cell_barcode, cluster) %>%
  summarise(n_umis = n_distinct(UB), .groups = "drop")

# ------------------------------------------------------------
# KEEP infected cells (>=2 UMIs)
# ------------------------------------------------------------
infected_cells <- cell_stats %>%
  filter(n_umis >= 2)

all_reads <- all_reads %>%
  inner_join(
    infected_cells %>% select(sample_id, cell_barcode),
    by = c("sample_id","cell_barcode")
  )

# ------------------------------------------------------------
# Compute positions
# ------------------------------------------------------------
all_reads <- all_reads %>%
  mutate(pos = round((start + end) / 2)) %>%
  filter(!is.na(pos)) %>%
  distinct(sample_id, cluster, pos, cell_barcode, UB)

# ------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------
window_size <- 500
step_size   <- 100
k <- window_size / step_size

# ------------------------------------------------------------
# Depth (cluster level)
# ------------------------------------------------------------
depth <- all_reads %>%
  mutate(bin = floor(pos / step_size)) %>%
  count(cluster, bin, name = "depth") %>%
  arrange(bin) %>%
  group_by(cluster) %>%
  mutate(
    depth = zoo::rollsum(depth, k = k, fill = NA, align = "center"),
    depth = replace_na(depth, 0),
    pos   = bin * step_size,
    log_depth = log10(depth + 1)
  ) %>%
  ungroup()

# ============================================================
# ✅ PCA 1 — SAMPLES
# ============================================================
depth_sample <- all_reads %>%
  mutate(bin = floor(pos / step_size)) %>%
  count(sample_id, bin, name = "depth") %>%
  arrange(bin) %>%
  group_by(sample_id) %>%
  mutate(
    depth = zoo::rollsum(depth, k = k, fill = NA, align = "center"),
    depth = replace_na(depth, 0),
    pos = bin * step_size,
    log_depth = log10(depth + 1)
  ) %>%
  ungroup()

mat_sample <- depth_sample %>%
  select(sample_id, pos, log_depth) %>%
  pivot_wider(names_from = pos, values_from = log_depth, values_fill = 0) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

# remove zero-variance columns
mat_sample <- mat_sample[, apply(mat_sample, 2, var) > 0]

pca_sample <- prcomp(mat_sample, scale. = TRUE)

var_exp_sample <- (pca_sample$sdev^2) / sum(pca_sample$sdev^2)

pca_sample_df <- as.data.frame(pca_sample$x) %>%
  rownames_to_column("sample")

p1 <- ggplot(pca_sample_df,
             aes(PC1, PC2, colour = sample)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = sample)) +
  scale_colour_viridis_d() +
  theme_minimal() +
  labs(
    title = "PCA: Samples",
    x = paste0("PC1 (", round(100 * var_exp_sample[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * var_exp_sample[2], 1), "%)")
  ) +
  theme(legend.position = "none")

ggsave("OUTPUT2/PCA_samples.png", p1, width = 6, height = 5)
ggsave("OUTPUT2/PCA_samples.pdf", p1, width = 6, height = 5)

# ============================================================
# ✅ PCA 2 — CLUSTERS
# ============================================================
mat_cluster <- depth %>%
  select(cluster, pos, log_depth) %>%
  pivot_wider(names_from = pos, values_from = log_depth, values_fill = 0) %>%
  column_to_rownames("cluster") %>%
  as.matrix()

mat_cluster <- mat_cluster[, apply(mat_cluster, 2, var) > 0]

pca_cluster <- prcomp(mat_cluster, scale. = TRUE)

var_exp_cluster <- (pca_cluster$sdev^2) / sum(pca_cluster$sdev^2)

pca_cluster_df <- as.data.frame(pca_cluster$x) %>%
  rownames_to_column("cluster")

p2 <- ggplot(pca_cluster_df,
             aes(PC1, PC2, colour = cluster)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = cluster)) +
  scale_colour_viridis_d() +
  theme_minimal() +
  labs(
    title = "PCA: Clusters",
    x = paste0("PC1 (", round(100 * var_exp_cluster[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * var_exp_cluster[2], 1), "%)")
  ) +
  theme(legend.position = "none")

ggsave("OUTPUT2/PCA_clusters.png", p2, width = 6, height = 5)
ggsave("OUTPUT2/PCA_clusters.pdf", p2, width = 6, height = 5)

# ============================================================
# ✅ PCA 3 — SAMPLE × CLUSTER
# ============================================================
depth_sc <- all_reads %>%
  mutate(bin = floor(pos / step_size)) %>%
  count(sample_id, cluster, bin, name = "depth") %>%
  arrange(bin) %>%
  group_by(sample_id, cluster) %>%
  mutate(
    depth = zoo::rollsum(depth, k = k, fill = NA, align = "center"),
    depth = replace_na(depth, 0),
    pos = bin * step_size,
    log_depth = log10(depth + 1)
  ) %>%
  ungroup()

mat_sc <- depth_sc %>%
  select(sample_id, cluster, pos, log_depth) %>%
  unite(id, sample_id, cluster, sep = "_") %>%
  pivot_wider(names_from = pos, values_from = log_depth, values_fill = 0) %>%
  column_to_rownames("id") %>%
  as.matrix()

mat_sc <- mat_sc[, apply(mat_sc, 2, var) > 0]

pca_sc <- prcomp(mat_sc, scale. = TRUE)

var_exp_sc <- (pca_sc$sdev^2) / sum(pca_sc$sdev^2)

pca_sc_df <- as.data.frame(pca_sc$x) %>%
  rownames_to_column("id") %>%
  separate(id, into = c("sample","cluster"), sep = "_")

p3 <- ggplot(pca_sc_df,
             aes(PC1, PC2, colour = cluster)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_text_repel(aes(label = paste0(sample, "_", cluster)),
                  size = 2.5, max.overlaps = 50) +
  scale_colour_viridis_d() +
  theme_minimal() +
  labs(
    title = "PCA: Sample × Cluster",
    x = paste0("PC1 (", round(100 * var_exp_sc[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * var_exp_sc[2], 1), "%)")
  )

ggsave("OUTPUT2/PCA_sample_cluster.png", p3, width = 7, height = 6)
ggsave("OUTPUT2/PCA_sample_cluster.pdf", p3, width = 7, height = 6)