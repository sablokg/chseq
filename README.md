# chseq

- my entire chipseq annotation and peak calling pipeline. 
- coded earlier while working at Poland and releasing now.

```
QC – FastQC on raw reads
Trimming – Trim Galore (adapters + quality)
Alignment – Bowtie2 (paired-end, sorted BAM)
Filtering/dedup – MAPQ≥30, remove chrM/unmapped/secondary reads, Picard MarkDuplicates, ENCODE blacklist removal
Peak calling – MACS2 (narrow for TFs, broad for histone marks — toggle with IS_NARROW)
Signal tracks – CPM-normalized bigWig + log2 ChIP/Input ratio track (deepTools)
QC metrics – flagstat, duplicate rate, FRiP score
MultiQC aggregate report at the end

```

Gaurav Sablok \
gsablok@proton.me
