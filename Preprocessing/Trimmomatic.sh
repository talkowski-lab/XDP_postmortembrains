#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=32G
#SBATCH --output=logs/log%J.out 
#SBATCH --error=logs/log%J.err 

#Script for trimming fastqs used for XDP postmortem brains RNA-Seq data

R1_file=$1
R2_file=$2
name=$3

module load java
mkdir -p ${name}

echo $name
java -jar /apps/lab/miket/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
    -phred33 \
    ${R1_file} ${R2_file} ${name}/${name}.R1.fq.gz ${name}/${name}.R1.unpaired.fq.gz ${name}/${name}.R2.fq.gz ${name}/${name}.R2.unpaired.fq.gz \
    ILLUMINACLIP:/apps/lab/miket/Trimmomatic-0.36/adapters/TruSeq3-PE-miket.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
