#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-54
#SBATCH --output=logs/log%J.out 
#SBATCH --error=logs/log%J.err 

#Script for RNASeqQC of alignments

config=$1
file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

filename="${file##*/}"
name="${filename%%.Aligned.sortedByCoord.out.bam}"

/data/talkowski/tools/bin/rnaseqc.v2.4.2.linux -s ${name} /data/talkowski/tools/ref/RNA-Seq/human/ref_components/GRCh38.112_svacorrected/Homo_sapiens.GRCh38.112.sva_corrected.genes.gtf ${file} .
