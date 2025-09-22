# XDP postmortembrains
Repository for preprocessing and analysis scripts for RNA-Seq data for the paper, "X-linked dystonia-parkinsonism is a genetic four-repeat tauopathy"

## RNA extraction, sequencing, and preprocessing						
RNA was extracted from brain pieces following bead-based homogenization using Trizol. Following quantification and quality control for RNA integrity using a Bioanalyzer, polyA-enriched mRNA-Seq libraries were constructed using the Illumina TruSeq Stranded mRNA kit. Libraries were sequenced using an Illumina NovaSeqX, with each library receiving a minimum of 25 million and a maximum of 107 million paired-end 150-bp reads. Raw fastq sequences were obtained and qualities were assessed using FastQC and RNA-SeqQC, adapters were trimmed using trimmomatic v0.36, then alignment to the human genome reference (GRCh38) was performed using STAR v2.7.10. Gene counts were obtained directly from STAR. Intron retention in TAF1 intron 32 was measured using IRFinder. [Link to scripts used for preprocessing](/Preprocessing).

## RNA-Seq differential expression analysis, enrichment analysis, and coexpression analysis
Differential expression analysis was performed in R using DESeq2 v1.46.0 for each brain region. Prior to obtaining DEGs, features were filtered for low expression (only genes with >1 cpm were included). Normalization was performed using DESeq median-of-ratios method, then principal component analyses were used to determine the suitability of models that incorporated surrogate variables (obtained via the SVAseq v3.54.0 package) to regress out unwanted technical variables from the differential expression analyses. DEGs were obtained at a predetermined false-discovery rate (FDR)-adjusted p-value cutoff of <0.01 and log2(fold-change) cutoff >0.58 (for upregulated genes) or log2(fold-change) <-0.58 (for downregulated genes). AnnotationDB v1.6.8 was used to convert features to gene symbols.

To determine direction-dependent overlaps between DEGs across brain regions, we used genes with adj-p<0.01 (regardless of fold-change) and tested for overlap via hypergeometric tests. For the analysis of genes that are dysregulated according to (CCCTCT)n repeat length, we first subsetted to genes that were already determined to be DEGs in the striatum. A linear model was then constructed using (CCCTCT)n repeat length (among XDP samples) and normalized expression as variables to determine the DEGs that show greater upregulation/downregulation with increasing (CCCTCT)n repeat length. Cutoff p-value was set to p<0.05. All functional enrichment analyses of gene lists (DEGs, overlapping DEGs, DEGs dysregulated according to repeat-number, module genes) were performed using the R package for gprofiler2 v0.2.3, using the query sets "GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC", "WP". Background was set to genes analyzed per brain region and p-values were FDR-adjusted. 				

Coexpression analyses were performed using the R package WGCNA v1.3 (Weighted Gene Network Coexpression Analysis). Normalized count matrices were log-transformed and then used to create signed networks from adjacency and topological overlap matrices, at soft power threshold of 20. Minimum gene module size was set to 50 genes. These modules were then refined to include only genes with module membership p<0.05 in their original assigned modules; otherwise, genes were reassigned to the unclassified module (M0). Module eigengenes were then correlated with traits (ME-trait correlations) and enrichment of DEGs within modules was tested using hypergeometric testing. [Link to analysis scripts](/Analysis)

## Bootstrapping
To avoid spurious results from having unbalanced XDP and control samples for differential expression analysis, results were overlapped with DEG discovery from ~50,000 bootstraps. [Link to bootstrapping analysis](/Bootstrapping)


## Data availability
Transcriptomic data are deposited in dbGap Project: phs001525.v2.p1. 
<img width="468" height="642" alt="image" src="https://github.com/user-attachments/assets/a2954bdc-30b1-489b-9df7-192a84c9a8fc" />

