#!/usr/bin/env bash
set -euo pipefail

############################################
# CONFIG
############################################
REF="EF523390.1.masked.fa"
THREADS=4

mkdir -p mdv_filtered
mkdir -p viral_realign

# ✅ Input BAMs (pos-sorted)
BAMS=(
MDV_02.bam
MDV_03.bam
MDV_04.bam
MDV_05.bam
MDV_06.bam
MDV_07.bam
MDV_11.bam
MDV_12.bam
)

############################################
# STEP 0 — CHECK + INDEX BAMS
############################################
echo "Checking and indexing BAMs..."

for BAM in "${BAMS[@]}"
do
    echo "Checking $BAM"

    if [ ! -f "$BAM" ]; then
        echo "ERROR: $BAM not found"
        exit 1
    fi

    # Check integrity
    samtools quickcheck "$BAM"

    # ✅ Ensure index exists
    if [ ! -f "${BAM}.bai" ]; then
        echo "Index missing → indexing $BAM"
        samtools index "$BAM"
    fi
done

echo "All BAMs OK and indexed"
echo "----------------------------------"

############################################
# MAIN LOOP
############################################
for BAM in "${BAMS[@]}"
do
    SAMPLE=$(basename "$BAM" .bam)

    echo ""
    echo "=================================="
    echo "Processing $SAMPLE"
    echo "=================================="

    ############################################
    # STEP 1 — Extract MDV reads
    ############################################
    echo "Extracting MDV references..."

    refs=$(samtools idxstats "$BAM" | cut -f1 | grep -E 'EF' | grep -v '^\*$' || true)

    if [ -z "$refs" ]; then
        echo "WARNING: No MDV refs found for $SAMPLE → skipping"
        continue
    fi

    echo "Found refs:"
    echo "$refs"

    samtools view -b "$BAM" $refs > mdv_filtered/${SAMPLE}.MDV.bam
    samtools index mdv_filtered/${SAMPLE}.MDV.bam

    ############################################
    # STEP 2 — Extract CB / UB / NH tags
    ############################################
    echo "Extracting CB/UB/NH..."

    samtools view mdv_filtered/${SAMPLE}.MDV.bam | \
    awk '{
        cb="NA"; ub="NA"; nh="NA";
        for(i=12;i<=NF;i++){
            if($i ~ /^CB:Z:/){cb=substr($i,6)}
            if($i ~ /^UB:Z:/){ub=substr($i,6)}
            if($i ~ /^NH:i:/){nh=substr($i,6)}
        }
        end=$4+length($10)-1;
        print cb"\t"ub"\t"$3"\t"$4"\t"end"\t"nh
    }' > mdv_filtered/${SAMPLE}.reads.tsv

    ############################################
    # STEP 3 — Summarise CB counts
    ############################################
    echo "Summarising CB usage..."

    cut -f1,2 mdv_filtered/${SAMPLE}.reads.tsv | \
        sort -u | \
        cut -f1 | \
        sort | \
        uniq -c > mdv_filtered/${SAMPLE}_cb_counts.txt

    ############################################
    # STEP 4 — Depth calculation
    ############################################
    echo "Computing depth..."

    samtools depth -aa mdv_filtered/${SAMPLE}.MDV.bam | \
        grep -E 'EF' > mdv_filtered/${SAMPLE}.depth

    ############################################
    # STEP 5 — BAM → FASTQ
    ############################################
    echo "Generating FASTQ..."

    samtools fastq mdv_filtered/${SAMPLE}.MDV.bam > mdv_filtered/${SAMPLE}.fastq
    ls -lt  mdv_filtered/${SAMPLE}.fastq

    ############################################
    # STEP 6 — Realign to viral reference
    ############################################
    echo "Realigning..."

    bwa mem -t ${THREADS} ${REF} \
        mdv_filtered/${SAMPLE}.fastq \
        > viral_realign/${SAMPLE}.sam

    ############################################
    # STEP 7 — Convert + sort + index
    ############################################
    samtools view -b viral_realign/${SAMPLE}.sam \
        > viral_realign/${SAMPLE}.bam

    samtools sort -o viral_realign/${SAMPLE}.sorted.bam \
        viral_realign/${SAMPLE}.bam

    samtools index viral_realign/${SAMPLE}.sorted.bam

    ############################################
    # CLEANUP (optional but recommended)
    ############################################
    rm viral_realign/${SAMPLE}.sam
    rm viral_realign/${SAMPLE}.bam

    echo "DONE $SAMPLE"
done

echo ""
echo "✅ ALL SAMPLES COMPLETE"
