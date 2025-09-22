#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48G
#SBATCH --array=1-55
#SBATCH --output=logs/log%J.out
#SBATCH --error=logs/log%J.err

#Script for detecting STMN2 exon2a coverage in RNA-Seq data of XDP postmortem brains

config=$1
file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

filename="${file##*/}"
name="${filename%%.Aligned.sortedByCoord.out.bam}"

module load samtools
module load bedtools2

WorkDir=/data/talkowski/Samples/XDP/Brain_RNASeq/v2/alignments
QCDir=/data/talkowski/Samples/XDP/Brain_RNASeq/v2/STMN2

samtools view -b $file 8:79611117-79666158 > ${name}_STMN2.bam
samtools index ${name}_STMN2.bam
bedtools coverage -F 0.5 -a exon2a.bed -b ${name}_STMN2.bam > ${name}_exon2a.tab
