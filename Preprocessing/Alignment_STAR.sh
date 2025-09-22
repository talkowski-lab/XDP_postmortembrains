#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --output=logs/log%J.out
#SBATCH --error=logs/log%J.err

#STAR alignment script used for RNA-Seq in XDP postmortem brains

name=$1
ALIGNDIR=/data/talkowski/Samples/XDP/Brain_RNASeq/v2/alignments
DATADIR=/data/talkowski/Samples/XDP/Brain_RNASeq/v2/data_trimmed
WORKDIR=$ALIGNDIR/${name}
mkdir -p $WORKDIR

module load STAR/2.7.10a
module load samtools

echo $name

STAR --runThreadN 8 \
     --genomeDir /data/talkowski/tools/ref/RNA-Seq/human/ref_components/GRCh38.112_svacorrected/STAR \
     --twopassMode Basic \
     --outFilterMultimapNmax 1 \
     --outFilterMismatchNoverLmax 0.05 \
     --outSAMtype BAM Unsorted \
     --outReadsUnmapped Fastx \
     --readFilesCommand zcat \
     --quantMode GeneCounts \
     --alignEndsType EndToEnd \
     --outFileNamePrefix $WORKDIR/${name}. \
     --readFilesIn $DATADIR/${name}/${name}.R1.fq.gz $DATADIR/${name}/${name}.R2.fq.gz; \
samtools sort -o $WORKDIR/${name}.Aligned.sortedByCoord.out.bam $WORKDIR/${name}*.Aligned.out.bam; \
samtools index $WORKDIR/${name}.Aligned.sortedByCoord.out.bam
