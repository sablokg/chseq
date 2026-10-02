#!/usr/bin/env bash
#===============================================================================
# ChIP-seq Analysis Pipeline
# Steps: QC -> Trimming -> Alignment -> Filtering/Dedup -> Peak Calling ->
#        Signal Tracks -> QC Metrics -> (optional) Motif & Annotation
#
# Requirements (install via conda/mamba recommended):
#   fastqc, multiqc, trim_galore (or fastp), bowtie2 or bwa, samtools,
#   picard, macs2, deeptools, bedtools, meme (optional), R + ChIPseeker (optional)
#
#   conda create -n chipseq -c bioconda -c conda-forge \
#       fastqc multiqc trim-galore bowtie2 bwa samtools picard macs2 \
#       deeptools bedtools meme
#
# Usage:
#   ./chipseq_pipeline.sh
#   (edit the CONFIG section below, or export the variables before running)
#===============================================================================

set -euo pipefail
set -x  # comment out to reduce verbosity

#-------------------------------------------------------------------------------
# 0. CONFIG — edit these for your experiment
#-------------------------------------------------------------------------------
THREADS=8

# Sample layout: one ChIP sample + one Input/IgG control (paired-end example)
SAMPLE="chip_sample"
CONTROL="input_control"

# Raw FASTQ files (paired-end; set R2 variables to "" for single-end)
CHIP_R1="raw/${SAMPLE}_R1.fastq.gz"
CHIP_R2="raw/${SAMPLE}_R2.fastq.gz"
CTRL_R1="raw/${CONTROL}_R1.fastq.gz"
CTRL_R2="raw/${CONTROL}_R2.fastq.gz"

# Genome / index
GENOME_FASTA="ref/genome.fa"
BOWTIE2_INDEX="ref/genome"        # prefix for bowtie2-build
EFFECTIVE_GENOME_SIZE="2.7e9"     # e.g. human hs=2.7e9, mouse mm=1.87e9

# Blacklist regions (ENCODE blacklist bed, recommended)
BLACKLIST="ref/blacklist.bed"

# MACS2 settings
MACS2_GENOME="hs"                 # hs, mm, ce, dm, or effective size number
IS_NARROW=true                    # true=narrowPeak (TF), false=broadPeak (histone marks)

# Output dirs
OUTDIR="results"
mkdir -p "$OUTDIR"/{qc/raw,qc/trimmed,trimmed,bam,peaks,bigwig,logs,multiqc}

#-------------------------------------------------------------------------------
# 1. RAW READ QC
#-------------------------------------------------------------------------------
fastqc -t "$THREADS" -o "$OUTDIR/qc/raw" \
    "$CHIP_R1" "$CHIP_R2" "$CTRL_R1" "$CTRL_R2"

#-------------------------------------------------------------------------------
# 2. ADAPTER / QUALITY TRIMMING (Trim Galore, paired-end)
#-------------------------------------------------------------------------------
trim_galore --paired --cores "$THREADS" --fastqc \
    -o "$OUTDIR/trimmed" "$CHIP_R1" "$CHIP_R2" \
    2> "$OUTDIR/logs/trim_${SAMPLE}.log"

trim_galore --paired --cores "$THREADS" --fastqc \
    -o "$OUTDIR/trimmed" "$CTRL_R1" "$CTRL_R2" \
    2> "$OUTDIR/logs/trim_${CONTROL}.log"

CHIP_TR1="$OUTDIR/trimmed/${SAMPLE}_R1_val_1.fq.gz"
CHIP_TR2="$OUTDIR/trimmed/${SAMPLE}_R2_val_2.fq.gz"
CTRL_TR1="$OUTDIR/trimmed/${CONTROL}_R1_val_1.fq.gz"
CTRL_TR2="$OUTDIR/trimmed/${CONTROL}_R2_val_2.fq.gz"

#-------------------------------------------------------------------------------
# 3. ALIGNMENT (Bowtie2 — standard for ChIP-seq)
#-------------------------------------------------------------------------------
align_sample () {
    local NAME=$1 R1=$2 R2=$3
    bowtie2 -p "$THREADS" -x "$BOWTIE2_INDEX" \
        -1 "$R1" -2 "$R2" \
        --very-sensitive -X 2000 \
        2> "$OUTDIR/logs/bowtie2_${NAME}.log" \
        | samtools sort -@ "$THREADS" -o "$OUTDIR/bam/${NAME}.sorted.bam" -
    samtools index "$OUTDIR/bam/${NAME}.sorted.bam"
}

align_sample "$SAMPLE" "$CHIP_TR1" "$CHIP_TR2"
align_sample "$CONTROL" "$CTRL_TR1" "$CTRL_TR2"

#-------------------------------------------------------------------------------
# 4. FILTER (mapping quality, unmapped, mitochondrial) + MARK/REMOVE DUPLICATES
#-------------------------------------------------------------------------------
filter_dedup () {
    local NAME=$1
    local IN="$OUTDIR/bam/${NAME}.sorted.bam"
    local FILT="$OUTDIR/bam/${NAME}.filt.bam"
    local DEDUP="$OUTDIR/bam/${NAME}.dedup.bam"

    # Remove unmapped, secondary, supplementary, low MAPQ, and chrM reads
    samtools view -@ "$THREADS" -b -q 30 -F 0x904 "$IN" \
        $(samtools view -H "$IN" | grep -oP '(?<=SN:)[^\s]+' | grep -v -i 'chrM\|MT') \
        > "$FILT"
    samtools index "$FILT"

    # Mark duplicates and remove them
    picard MarkDuplicates \
        I="$FILT" O="$DEDUP" \
        M="$OUTDIR/logs/${NAME}_dup_metrics.txt" \
        REMOVE_DUPLICATES=true

    samtools index "$DEDUP"

    # Remove ENCODE blacklist regions if provided
    if [[ -f "$BLACKLIST" ]]; then
        bedtools intersect -v -abam "$DEDUP" -b "$BLACKLIST" \
            > "$OUTDIR/bam/${NAME}.final.bam"
    else
        cp "$DEDUP" "$OUTDIR/bam/${NAME}.final.bam"
    fi
    samtools index "$OUTDIR/bam/${NAME}.final.bam"
    samtools flagstat "$OUTDIR/bam/${NAME}.final.bam" > "$OUTDIR/logs/${NAME}_flagstat.txt"
}

filter_dedup "$SAMPLE"
filter_dedup "$CONTROL"

CHIP_BAM="$OUTDIR/bam/${SAMPLE}.final.bam"
CTRL_BAM="$OUTDIR/bam/${CONTROL}.final.bam"

#-------------------------------------------------------------------------------
# 5. FRAGMENT / LIBRARY QC (cross-correlation, complexity — optional but recommended)
#-------------------------------------------------------------------------------
# NSC/RSC via phantompeakqualtools (install separately if desired):
# run_spp.R -c=$CHIP_BAM -savp -out="$OUTDIR/logs/${SAMPLE}_spp.out"

#-------------------------------------------------------------------------------
# 6. PEAK CALLING (MACS2)
#-------------------------------------------------------------------------------
if [[ "$IS_NARROW" == true ]]; then
    macs2 callpeak \
        -t "$CHIP_BAM" -c "$CTRL_BAM" \
        -f BAMPE -g "$MACS2_GENOME" \
        -n "$SAMPLE" --outdir "$OUTDIR/peaks" \
        -q 0.05 --keep-dup all \
        2> "$OUTDIR/logs/macs2_${SAMPLE}.log"
else
    macs2 callpeak \
        -t "$CHIP_BAM" -c "$CTRL_BAM" \
        -f BAMPE -g "$MACS2_GENOME" \
        -n "$SAMPLE" --outdir "$OUTDIR/peaks" \
        --broad --broad-cutoff 0.1 --keep-dup all \
        2> "$OUTDIR/logs/macs2_${SAMPLE}.log"
fi

#-------------------------------------------------------------------------------
# 7. SIGNAL TRACKS (bigWig, CPM-normalized) for genome browser visualization
#-------------------------------------------------------------------------------
bamCoverage -b "$CHIP_BAM" -o "$OUTDIR/bigwig/${SAMPLE}.bw" \
    --normalizeUsing CPM --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" \
    -p "$THREADS" --extendReads

bamCoverage -b "$CTRL_BAM" -o "$OUTDIR/bigwig/${CONTROL}.bw" \
    --normalizeUsing CPM --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" \
    -p "$THREADS" --extendReads

# Log2 ratio ChIP vs Input track (nice for visualization)
bigwigCompare -b1 "$OUTDIR/bigwig/${SAMPLE}.bw" -b2 "$OUTDIR/bigwig/${CONTROL}.bw" \
    --operation log2 -o "$OUTDIR/bigwig/${SAMPLE}_vs_${CONTROL}_log2ratio.bw" \
    -p "$THREADS"

#-------------------------------------------------------------------------------
# 8. PEAK QC: FRiP score (Fraction of Reads in Peaks)
#-------------------------------------------------------------------------------
PEAKFILE="$OUTDIR/peaks/${SAMPLE}_peaks.$( [[ "$IS_NARROW" == true ]] && echo narrowPeak || echo broadPeak )"
TOTAL_READS=$(samtools view -c "$CHIP_BAM")
READS_IN_PEAKS=$(bedtools intersect -a "$CHIP_BAM" -b "$PEAKFILE" -u -bed | wc -l)
FRIP=$(echo "scale=4; $READS_IN_PEAKS / $TOTAL_READS" | bc)
echo -e "Sample\tTotalReads\tReadsInPeaks\tFRiP\n${SAMPLE}\t${TOTAL_READS}\t${READS_IN_PEAKS}\t${FRIP}" \
    > "$OUTDIR/logs/${SAMPLE}_FRiP.txt"

#-------------------------------------------------------------------------------
# 9. (OPTIONAL) MOTIF ENRICHMENT
#-------------------------------------------------------------------------------
# bedtools getfasta -fi "$GENOME_FASTA" -bed "$PEAKFILE" -fo "$OUTDIR/peaks/${SAMPLE}_peaks.fa"
# meme-chip -oc "$OUTDIR/peaks/${SAMPLE}_memechip" -db motif_databases/JASPAR/JASPAR2022_CORE.meme \
#     "$OUTDIR/peaks/${SAMPLE}_peaks.fa"

#-------------------------------------------------------------------------------
# 10. (OPTIONAL) PEAK ANNOTATION — run separately in R
#-------------------------------------------------------------------------------
# library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# peaks <- readPeakFile("results/peaks/chip_sample_peaks.narrowPeak")
# peakAnno <- annotatePeak(peaks, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, tssRegion=c(-3000,3000))
# plotAnnoPie(peakAnno)

#-------------------------------------------------------------------------------
# 11. AGGREGATE QC REPORT
#-------------------------------------------------------------------------------
multiqc "$OUTDIR" -o "$OUTDIR/multiqc"

echo "Pipeline complete. Peaks: $PEAKFILE"
echo "Signal tracks: $OUTDIR/bigwig/"
echo "QC report: $OUTDIR/multiqc/multiqc_report.html"
