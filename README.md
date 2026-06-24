This pipeline extracts Marek’s disease virus (MDV) reads from single-cell RNA-seq (Cell Ranger) BAM files and quantifies viral transcription at the level of cells and genes.

# [1] The 01_processing.sh ##

Inputs: MDV_{ID}.bam (Cell Ranger BAM files), MDV.fasta (viral reference genome), MDV.gtf or MDV.gff (annotation file)

# The 01_processing.sh script does the following:
- extraction of MDV reads from BAM files
- extraction of cell barcodes (CB) and UMIs (UB)
- assignment of reads to genes using midpoint coordinates
- removal of PCR duplicates using barcode and UMI
- calculation of per-cell viral abundance
- calculation of global viral gene expression
- Each read is assigned a single genomic position: midpoint = floor((start + end) / 2). This prevents one read from being counted multiple times across overlapping features.
-  PCR duplicate removal: Reads are collapsed using: cell barcode + UMI - This ensures that each original RNA molecule is counted once.

Key outputs per sample
- MDV__reads_with_CB_UB.bed   # Contains each read with genomic position, barcode, and UMI.
- MDV__reads_midpoints.bed    # Contains one position per read (midpoint).
- MDV__midpoints_to_genes.tsv # Contains mapping of each read (via midpoint) to a gene.
- MDV__cell_barcode_viral_reads.tsv

Per-cell viral counts:
- cell_barcode
- n_reads (total viral reads)
- n_unique_umis (PCR-collapsed molecules)

MDV_{ID}_gene_expression.tsv: 
- gene_id
- total_UMIs (PCR-corrected)
- fraction (relative expression)

# [2] Rscript 02_genomic_depth.R

# requires MDV_1.depth from 01_processing.sh, and GFF and GTF for MDV
# plots depth per sample

# [3] Rscript 03_cluster.R

# links viral counts to cell types, samples, UMIs, barcode
# gives you the numbers of infected cells per sample per cell type
# makes MDV_cluster_sample_full_summary.tsv

# [4] Rscript 04_p.R

# reads in MDV_cluster_sample_full_summary.tsv
# viusalises fraction_cells vs UMI rate