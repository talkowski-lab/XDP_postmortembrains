#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=64G
#SBATCH --output=logs/log%J.out
#SBATCH --error=logs/log%J.err

#Script for running IRFinder in RNA-Seq data from XDP postmortem brains

module load SAMTools/1.10 
module load STAR/2.7.10a
module load bedtools/2.20.1

file=$1
name=$2

REF=/data/talkowski/tools/ref/RNA-Seq/human/ref_components/GRCh38.112_svacorrected/IRFinder

mkdir ${name}
echo ${name}

/data/talkowski/dg520/app/IRFinder/bin/IRFinder -m BAM -r ${REF} -d ${name} \
  ${file}
