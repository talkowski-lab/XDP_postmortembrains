#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type BEGIN
#SBATCH --mail-user "adomingo1@mgh.harvard.edu"
#SBATCH --output=logs/log%J.out 
#SBATCH --error=logs/log%J.err 

#Script for RMATS pre-processing of RNA-Seq data of XDP postmortem brains

bfile=$1
module load rmats

gtf="/data/talkowski/tools/ref/RNA-Seq/human/ref_components/GRCh38.112_svacorrected/Homo_sapiens.GRCh38.112.sva_corrected.gtf"
WORKDIR="/data/talkowski/Samples/XDP/Brain_RNASeq/v2/RMATS"

mkdir $WORKDIR/tmp/$bfile

rmats --b1 $bfile --gtf $gtf -t paired --readLength 150 --libType fr-firststrand --nthread 8 --od $WORKDIR/rmats_pre --tmp $WORKDIR/tmp/$bfile --task prep
