#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/log%J.out
#SBATCH --error=logs/log%J.err

#Script to check coverage of TAF1 exons

module load Java/1.8.0_191

ID=$1

REF=/data/talkowski/tools/ref/RNA-Seq/human/ref_components/GRCh38.112_svacorrected/Homo_sapiens.GRCh38.dna.primary_assembly.ercc.sva_corrected.fa
Target=/data/talkowski/Samples/XDP/Brain_RNASeq/v2/scripts/TAF1_exons_GRCh38.intervals

java -jar /apps/lab/miket/picard/2.18.11/picard.jar CollectHsMetrics BAIT_INTERVALS=$Target TARGET_INTERVALS=$Target INPUT=/data/talkowski/Samples/XDP/Brain_RNASeq/v2/alignments/${ID}/${ID}.Aligned.sortedByCoord.out.bam OUTPUT=${ID}_HsMetrics.tab PER_TARGET_COVERAGE=${ID}_perexoncov.tab  REFERENCE_SEQUENCE=${REF} VALIDATION_STRINGENCY=LENIENT
