#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type BEGIN
#SBATCH --mail-user "adomingo1@mgh.harvard.edu"
#SBATCH --output=logs/log%J.out 
#SBATCH --error=logs/log%J.err 

#Script for RMATS post-processing

bfile1=$1
bfile2=$2
module load rmats

gtf="/data/talkowski/tools/ref/RNA-Seq/human/ref_components/GRCh38.112_svacorrected/Homo_sapiens.GRCh38.112.sva_corrected.gtf"
WORKDIR="/data/talkowski/Samples/XDP/Brain_RNASeq/v2/RMATS"

rmats --b1 $bfile1 --b2 $bfile2 --gtf $gtf -t paired --readLength 150 --nthread 8 --od $WORKDIR/rmats_post --tmp $WORKDIR/rmats_pre --task post --novelSS
