# Generate a bunch of tsv files, each line in each tsv corresponds to a pair
# of control vs case samples
################
# Parse Arguments
library(magrittr)
library(optparse)
parsed_args <- 
    OptionParser() %>%
    add_option(
        c('-t', '--tissues'),
        help='Tissue samples to use',
        action='store',
        type='character',
        dest='TISSUES',
        default=NA
    ) %>% 
    add_option(
        c('-s', '--sample_metadata'),
        help='sample metadata file name',
        action='store',
        type='character',
        dest='SAMPLE_METADATA_FILE',
        default=NA
    ) %>% 
    add_option(
        c('-c', '--counts'),
        help='counts matrix file',
        action='store',
        type='character',
        dest='COUNTS_FILE',
        default=NA
    ) %>% 
    add_option(
        c('-o', '--output_dir'),
        help='Where to save files of DEG results',
        action='store',
        type='character',
        dest='OUTPUT_DIR',
        default=NA
    ) %>% 
    add_option(
        c('-f', '--force_redo'),
        help='Where to save files of DEG results',
        action='store_true',
        type='logical',
        dest='FORCE_REDO',
        default=FALSE
    ) %>% 
    add_option(
        c('-g', '--gene_id_column'),
        help='',
        action='store',
        type='character',
        dest='GENE_ID_COLUMN',
        default='EnsemblID'
    ) %>% 
    add_option(
        c('-i', '--sample_id_column'),
        help='',
        action='store',
        type='character',
        dest='SAMPLE_ID_VAR',
        default='SampleID'
    ) %>% 
    parse_args(positional_arguments=TRUE)
options <- parsed_args$options
model.files <- parsed_args$args
if (length(model.files) == 0) {
    model.files <- NULL
}
#################
#################
# Dependencies
suppressPackageStartupMessages({
    library(here)
    here::i_am('scripts/pick_bootstrap_samples.R')
    BASE_DIR <- here()
    SCRIPT_DIR <- here('scripts')
    source(file.path(SCRIPT_DIR, 'locations.R'))
    source(file.path(SCRIPT_DIR, 'bootstrap_utils.R'))
    library(future)
    library(furrr)
    library(tidyverse)
})
################
# Load input Data
# Load sample metadata
sample.metadata <- 
    if (is.na(options$SAMPLE_METADATA_FILE)) {
        load_sample_metadata()  # in bootstrap.utils.R
    } else {
        options$SAMPLE_METADATA_FILE %>% 
        read_tsv(show_col_types=FALSE)
    }
# Filter tissues
if (is.na(options$TISSUES)) {
    tissues <- unique(sample.metadata$Region)
    options$OUTPUT_DIR <- file.path(options$OUTPUT_DIR, 'All')
} else {
    tissues <- str_split(options$TISSUES, ',') %>% unlist()
    options$OUTPUT_DIR <- 
        file.path(
            options$OUTPUT_DIR, 
            str_replace_all(options$TISSUES, ',', '+')
        )
}
sample.metadata <- sample.metadata %>% filter(Region %in% tissues)
sample.metadata %>% dplyr::count(Region)
message('loaded sample metadata')
print(glue('Building model with {nrow(sample.metadata)} total samples from {options$TISSUES}'))
# Load counts
counts.matrix <- 
    {
        if (is.na(options$COUNTS_FILE)) {
            load_counts()  # in bootstrap.utils.R
        } else {
            options$COUNTS_FILE %>% 
            read_tsv(show_col_types=FALSE)
        }
    } %>%
    dplyr::select(
        c(
            options$SAMPLE_ID_VARIABLE,
            sample.metadata$SampleID
        )
    )
print(dim(sample.metadata))
print(dim(counts.matrix))
# Load model params
# Should be 1 row per model, each column specifying a parameter/threshold
model.parameters <- 
    load_model_specifications(
        filepaths=model.files,
        output_dir=options$OUTPUT_DIR,
        sample.metadata=sample.metadata
    )
################
# Compute DEG results using all samples for each model
# Build DESeq2 model + get DEG results for each model using all samples at once
model.parameters %>%
    # filter(model.name %in% c('blank', 'standard', 'SVA.only')) %>% 
    mutate(results_file=file.path(model.dir, 'full.model.results.tsv')) %>% 
    rowwise() %>% 
    pmap(
        .l=.,
        .f=check_cached_results,
        results_fnc=
            function(sample_metadata, counts_matrix, ...) {
                deseq_pipeline(
                    sample_metadata=sample_metadata,
                    counts_matrix=counts_matrix,
                    ...
                ) %>% 
                mutate(contrast=glue('{variable}-{numerator}-{denominator}')) %>% 
                select(-c(name)) %>% 
                # rename('contrast'=name) %>% 
                unnest(deg.results)
            },
        sample_metadata=sample.metadata,
        counts_matrix=counts.matrix,
        force_redo=options$FORCE_REDO,
        return_data=FALSE,
        sampleID.variable=options$SAMPLE_ID_VAR,
        gene_id_column=options$GENE_ID_COLUMN,
        .progress=TRUE
    )
