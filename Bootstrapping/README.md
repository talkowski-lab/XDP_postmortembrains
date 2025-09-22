# Bootstrapping Results

The main commands to do this analysis: 

```bash
# Generate bootstrap sample sets
$ Rscript ./scipts/pick_bootstrap_samples.R 
    --tissues CAU 
    --chunk_size 200 
    --n_bootstraps 50000
    --output_dir ./results/bootstraps/CAU/SVA.only/ 
# Calculate bootstrap DESeq2 results
```

### Bootstrapping dependencies

Install R dependencies for running `DESeq2`
```
# R dependencies
install.packages(
    c(
        'glue',
        'optparse',
        'future',
        'furrr',
        'tictoc',
        'magrittr',
        'BiocManager',
        'tidyverse'
    )
)
BiocManager::install(
    c(
        'sva',
        'DESeq2'
    )
)
```

### Model Forumulations

All the relevant model biulding parameters are in the file .`/models/SVA.only.model`
```
model.name="SVA.only"
effect.variable="Genotype"      # Main effect variable in the model (XDP vs CON)
include.SVs=TRUE                # Estimate Surrogate Variable via SVA + include in model
covariate.variables=""          # No covariates for the model (use SVs)
filter.variable="Genotype"      # Remove genes with low expression in either Genotype
CPM.cutoff=0.1                  
min.frac.samples.expressed=0.5
intercept=1                     # Include intercept term in DESeq2 models
n.samples.per.group=6           # Each bootstrap has 6 samples per group (Genotype)
seed=9                          # Random seed for picking 
```

## Output Format

For each model we pick a number of and total bootstraps and split them into chunks i.e. 100 bootstraps -> 5 chunks of 20 bootstraps with all bootstraps in a single file e.g. `chunk1,tsv`. 
The file tree looks like this 

```{bash}
$ tree results/bootstraps/standard/
results/bootstraps/standard/
├── bootstrap.parameters.txt    # pacakge session info + params 
├── this.model                  # model parameters file
├── bootstrap.samples.lists     # list of all samples in each bootstrap
│   ├── chunk1.tsv
│   ├── chunk2.tsv
│   ├── chunk3.tsv
│   ├── chunk4.tsv
│   └── chunk5.tsv
└── DEG.results                 # all DESeq2 results for all bootstraps, grouped by chunk
    ├── Genotype-XDP-Control    # specific contrast results
    │   ├── chunk1.tsv          # DESeq2 results for all bootstraps in chunk1
    │   ├── chunk2.tsv
    │   ├── chunk3.tsv
    │   ├── chunk4.tsv
    │   └── chunk5.tsv
    ├── SV1-NA-NA               # DESeq2 results for the Surrogate Variable contrasts (unused)
    │   ├── chunk1.tsv
    │   └── ...
    ├── ...
    └── ...
        
```

For each contrast the results are also organized per chunk, combinng bootstraps results into a single table with 1 column per bootstrap i.e. the file `results/bootstraps/standard/DEG.results/Genotype-XDP-Control/chunk1.tsv` looks like this:

| GeneName | statistic | bs1  | bs2  | bs3  | ... | bsN  |
| -------- | --------- | ---- | ---- | ---- | --- | ---- |
| gene_A   | LFC       | 1.1  | 1.2  | 0.7  | ... | 1.1  |
| gene_A   | pvalue    | 0.1  | 0.7  | 0.01 | ... | 0.04 |
| gene_B   | LFC       | 2.1  | 2.2  | 1.7  | ... | 1.1  |
| gene_B   | pvalue    | 0.05 | 0.01 | 0.01 | ... | 0.4  |
| ...      | ...       | ...  | ...  | ...  | ... | ...  |
| gene_Z   | pvalue    | 0.15 | 0.41 | 0.09 | ... | 0.01 |
| gene_Z   | LFC       | 1.5  | 1.4  | 1.7  | ... | 1.3  |
| gene_B   | pvalue    | 0.05 | 0.01 | 0.01 | ... | 0.4  |

These files can easily be aggregated by column-joining all the chunks together

Given all the individual bootstraps we can calculate the frequency of expression patterns across all bootstraps e.g. how frequently is a gene nominally upregulated i.e. log2FoldChange > 0 && pvalue < 0.05.

