#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48G
#SBATCH --array=1-47
#SBATCH --output=logs/log%J.out
#SBATCH --error=logs/log%J.err

#Script for QC of alignments

config=$1
file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

filename="${file##*/}"
name="${filename%%.Aligned.sortedByCoord.out.bam}"

module load Java/1.8.0_191

WorkDir=/data/talkowski/Samples/XDP/Brain_RNASeq/v2/alignments
QCDir=/data/talkowski/Samples/XDP/Brain_RNASeq/v2/QC/CollectRNASeqMetrics

/data/talkowski/tools/bin/gatk-4.1.8.1/gatk --java-options "-Xmx48G" CollectRnaSeqMetrics \
	-I ${file} \
	-O ${QCDir}/${SN}.RNASeqMetrics.txt \
	--REF_FLAT /data/talkowski/tools/ref/RNA-Seq/human/ref_components/GRCh38.112_svacorrected/Homo_sapiens.GRCh38.112.sva_corrected.refFlat \
	-STRAND SECOND_READ_TRANSCRIPTION_STRAND
