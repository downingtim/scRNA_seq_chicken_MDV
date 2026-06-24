#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# MDV SC-RNA PIPELINE (FINAL STABLE VERSION)
# ============================================================

SAMPLES="01 02 03 04 05 06 07 08"
SAMPLES="02"
ANNOT="MDV.gtf"

echo "======================================"
echo "Running MDV pipeline"
echo "======================================"

for ID in $SAMPLES; do

    echo ""
    echo "=============================="
    echo "Processing sample $ID"
    echo "=============================="

    BAM_IN="MDV_${ID}.bam"
    BAM_OUT="MDV.${ID}.bam"
    PREFIX="MDV_${ID}"

    # --------------------------------------------------------
    # 1. Extract viral reads
    # --------------------------------------------------------
    echo "[1/6] Extracting viral reads"

    samtools index "$BAM_IN"

    samtools idxstats "$BAM_IN" | cut -f1 | grep "alli" | grep -v '^\*$' > tmp.refs.txt
    refs=$(cat tmp.refs.txt)

    samtools view -b "$BAM_IN" $refs > "$BAM_OUT"
    samtools index "$BAM_OUT"

    # --------------------------------------------------------
    # 2. Correct depth (ONLY MDV contig)
    # --------------------------------------------------------
    echo "[2/6] Generating depth file"

    MDV_CHROM=$(head -n1 tmp.refs.txt)

    echo "Using contig: $MDV_CHROM"

    samtools depth -aa -r "$MDV_CHROM" "$BAM_OUT" > "${PREFIX}.depth"

    echo "Depth lines:"
    wc -l "${PREFIX}.depth"

    # --------------------------------------------------------
    # 3. Build BED with CB and UB
    # --------------------------------------------------------
    echo "[3/6] Extracting CB/UB"

    samtools view -F 2308 "$BAM_OUT" | awk 'BEGIN{OFS="\t"}
    {
        chrom=$3
        start=$4-1
        end=start+length($10)

        cb="NA"; ub="NA"

        for(i=12;i<=NF;i++){
            if($i ~ /^CB:Z:/) cb=substr($i,6)
            if($i ~ /^UB:Z:/) ub=substr($i,6)
        }

        if(cb=="NA" || ub=="NA") next

        print chrom, start, end, $1, $5, ".", cb, ub
    }' | sort -k1,1 -k2,2n > "${PREFIX}_reads_with_CB_UB.bed"

    echo "Reads with CB+UB:"
    wc -l "${PREFIX}_reads_with_CB_UB.bed"

    # --------------------------------------------------------
    # 4. Midpoints
    # --------------------------------------------------------
    echo "[4/6] Creating midpoints"

    awk 'BEGIN{OFS="\t"}
    {
        mid=int(($2+$3)/2)
        print $1, mid, mid+1, $4, $5, $6, $7, $8
    }' "${PREFIX}_reads_with_CB_UB.bed" \
    > "${PREFIX}_reads_midpoints.bed"

    BAM_CHROM=$(head -n1 "${PREFIX}_reads_midpoints.bed" | cut -f1)

    echo "Using chromosome for genes: $BAM_CHROM"

    # --------------------------------------------------------
    # 5. Gene assignment (FIXED)
    # --------------------------------------------------------
    echo "[5/6] Assigning reads to genes"

    awk -v chrom="$BAM_CHROM" 'BEGIN{FS=OFS="\t"}
    $3=="gene"{
        start=$4-1
        gene="NA"

        if(match($9,/gene_id[ ="]+([^";]+)/,a)) gene=a[1]
        else if(match($9,/Name=([^;]+)/,a)) gene=a[1]
        else if(match($9,/ID=([^;]+)/,a)) gene=a[1]

        print chrom, start, $5, gene, ".", $7
    }' "$ANNOT" > "${PREFIX}_genes.bed"

    bedtools intersect -wa -wb \
        -a "${PREFIX}_reads_midpoints.bed" \
        -b "${PREFIX}_genes.bed" \
    > "${PREFIX}_midpoints_to_genes.tsv"

    echo "Gene assignments:"
    wc -l "${PREFIX}_midpoints_to_genes.tsv"

    # --------------------------------------------------------
    # 6a. Per-cell viral reads
    # --------------------------------------------------------
    echo "[6/6] Calculating summaries"

    awk 'BEGIN{OFS="\t"}
    {
        cb=$7; ub=$8

        read[cb]++

        key=cb OFS ub
        if(!(key in seen)){
            seen[key]=1
            umi[cb]++
        }
    }
    END{
        for(cb in read){
            print cb, read[cb], umi[cb]
        }
    }' "${PREFIX}_reads_with_CB_UB.bed" \
    | sort -k2,2nr > "${PREFIX}_cell_barcode_viral_reads.tsv"

    sed -i '1i cell_barcode\tn_reads\tn_unique_umis' \
    "${PREFIX}_cell_barcode_viral_reads.tsv"

    # --------------------------------------------------------
    # 6b. Global gene expression (FINAL FIX)
    # --------------------------------------------------------
    awk 'BEGIN{OFS="\t"}
    {
        gene=$12
        cb=$7
        ub=$8

        if(gene=="" || gene==".") next

        key=gene OFS cb OFS ub

        if(!(key in seen)){
            seen[key]=1
            gene_count[gene]++
            total++
        }
    }
    END{
        if(total==0){
            print "WARNING: no gene counts" > "/dev/stderr"
        }

        for(g in gene_count){
            frac=gene_count[g]/total
            print g, gene_count[g], frac
        }
    }' "${PREFIX}_midpoints_to_genes.tsv" \
    | sort -k2,2nr > "${PREFIX}_gene_expression.tmp"

    echo -e "gene_id\ttotal_UMIs\tfraction" \
    > "${PREFIX}_gene_expression.tsv"

    cat "${PREFIX}_gene_expression.tmp" \
    >> "${PREFIX}_gene_expression.tsv"

    rm -f "${PREFIX}_gene_expression.tmp"

    echo "✅ Sample $ID complete"

done

rm -f tmp.refs.txt

echo ""
echo "======================================"
echo "✅ ALL SAMPLES COMPLETE"
echo "======================================"

echo "Outputs per sample:"
echo "  *_cell_barcode_viral_reads.tsv"
echo "  *_gene_expression.tsv"
echo "  *.depth (correct ~178k rows)"
