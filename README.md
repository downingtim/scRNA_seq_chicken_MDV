# MDV Single-Cell Viral Genome Analysis Pipeline

Pipeline for analysing Marek's Disease Virus (MDV) reads from single-cell sequencing datasets.

The workflow:

1. Extracts MDV reads from aligned BAM files
2. Generates cell barcode and UMI level viral read tables
3. Realigns viral reads against an MDV reference genome
4. Produces whole-genome viral coverage plots
5. Quantifies infection across cell clusters
6. Generates cluster-level viral genome profiles
7. Performs PCA on viral genomic coverage patterns
8. Creates gene-level heatmaps across samples and clusters

## Repository Contents

- `01_processing.sh`
- `02_genomic_depth.R`
- `03_cluster.R`
- `04_cluster_samples.R`
- `05_PCA.R`
- `06_clusters.R`

## Requirements

### Command-line software

- bash
- samtools
- bwa
- awk

### R packages

- tidyverse
- cowplot
- scales
- zoo
- viridis
- ggrepel
- ggplot2

Install required packages:

```r
install.packages(c(
  "tidyverse",
  "cowplot",
  "scales",
  "zoo",
  "viridis",
  "ggrepel",
  "ggplot2"
))
```

## Required Input Files

### BAM files

The pipeline expects position-sorted BAM files:

- MDV_02.bam
- MDV_03.bam
- MDV_04.bam
- MDV_05.bam
- MDV_06.bam
- MDV_07.bam
- MDV_11.bam
- MDV_12.bam

Index files will be generated automatically if missing.

### Viral reference

- EF523390.1.masked.fa

BWA index files must already exist.

### Annotation files

Depending on the analysis step:

- MDV.gtf
- MDV.gff
- EF523390.1.masked.gtf
- EF523390.1.masked.gff
- EF523390.1.masked.gb

### Cell cluster annotations

- cell_barcode_sample_cluster.csv

Expected columns:

- Barcode
- Sample
- cluster

## Step 1: Read Processing

Script:

```bash
bash 01_processing.sh
```

This script:

- Checks BAM integrity
- Creates BAM indexes if required
- Extracts MDV-aligned reads
- Extracts CB, UB and NH tags
- Generates per-read tables
- Produces depth files
- Converts BAM to FASTQ
- Realigns reads to the MDV reference
- Produces sorted and indexed viral BAM files

Output directories:

- mdv_filtered/
- viral_realign/

Key outputs:

- mdv_filtered/*.MDV.bam
- mdv_filtered/*.reads.tsv
- mdv_filtered/*.depth
- mdv_filtered/*.fastq
- viral_realign/*.sorted.bam

## Step 2: Genome Coverage by Sample

Script:

```r
source("02_genomic_depth.R")
```

Generates sliding-window viral genome coverage plots for each sample.

Features:

- UMI deduplication
- Cell filtering
- Minimum infection threshold of 10 UMIs per cell
- Sliding window smoothing
- Viral gene annotation track
- Highlighted meq, CXC chemokine and ICP4 regions

Outputs:

- plot_grouped_UMI_sliding_filtered.png
- plot_grouped_UMI_sliding_filtered.pdf

## Step 3: Cluster Infection Analysis

Script:

```r
source("03_cluster.R")
```

Calculates:

- Total cells per cluster
- Infected cells per cluster
- Fraction infected
- Total viral UMIs
- Total viral reads
- UMIs per infected cell

Infection threshold:

```text
> 10 unique viral UMIs
```

Outputs:

- MDV_cluster_sample_full_summary.tsv
- fraction_vs_UMIs_by_sample_LABELLED.png
- fraction_vs_UMIs_by_sample_LABELLED.pdf

## Step 4: Cluster-Level Genome Profiles

Script:

```r
source("04_cluster_samples.R")
```

Builds whole-genome viral coverage tracks for each cell cluster.

Key parameters:

```text
UMI threshold = 10
Window size   = 500 bp
Step size     = 100 bp
```

Outputs:

- plot_cluster_infection_sorted.png
- plot_cluster_infection_sorted.pdf

Each panel reports:

- Cluster identity
- Fraction of infected cells
- Mean viral UMIs per infected cell

## Step 5: Principal Component Analysis

Script:

```r
source("05_PCA.R")
```

Performs PCA on viral genomic coverage profiles.

Control samples 07 and 11 are excluded.

Generated analyses:

### PCA by sample

Outputs:

- OUTPUT2/PCA_samples.png
- OUTPUT2/PCA_samples.pdf

### PCA by cluster

Outputs:

- OUTPUT2/PCA_clusters.png
- OUTPUT2/PCA_clusters.pdf

### PCA by sample-cluster combination

Outputs:

- OUTPUT2/PCA_sample_cluster.png
- OUTPUT2/PCA_sample_cluster.pdf

## Step 6: Gene-Level Heatmaps

Script:

```r
source("06_clusters.R")
```

Parses CDS features from:

```text
EF523390.1.masked.gb
```

For each viral gene the pipeline:

- Calculates rolling-window UMI coverage
- Determines peak signal
- Computes log2-transformed abundance
- Generates heatmaps across samples and clusters

Outputs:

- OUTPUT4/table_cluster_heatmap.tsv
- OUTPUT4/table_sample_heatmap.tsv
- OUTPUT4/heatmap_by_cluster.png
- OUTPUT4/heatmap_by_sample.png

## Sample Metadata

- 02 = Spleen-1
- 12 = Spleen-2
- 03 = Kidney-1
- 05 = Kidney-2
- 04 = Testes-1
- 06 = Testes-2
- 07 = Control-1
- 11 = Control-2

## Typical Workflow

Run processing:

```bash
bash 01_processing.sh
```

Generate sample genome plots:

```r
source("02_genomic_depth.R")
```

Calculate infection metrics:

```r
source("03_cluster.R")
```

Generate cluster genome profiles:

```r
source("04_cluster_samples.R")
```

Run PCA analyses:

```r
source("05_PCA.R")
```

Generate heatmaps:

```r
source("06_clusters.R")
```

## Notes

- Viral reads are assigned using MDV reference contigs matching `EF*`.
- Cell barcodes are extracted from the `CB` BAM tag.
- UMIs are extracted from the `UB` BAM tag.
- Most analyses define an infected cell as having at least 10 unique viral UMIs.
- Genomic coverage is based on midpoint positions of deduplicated viral UMIs.
- Viral coordinates are based on MDV reference accession `EF523390.1`.