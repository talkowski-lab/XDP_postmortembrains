# Generate a bunch of tsv files, each line in each tsv corresponds to a pair
# of control vs case samples
################
# Parse Arguments
library(optparse)
library(stringr)
library(magrittr)
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
        c('-o', '--output_dir'),
        help='Where to save files of samples lists',
        action='store',
        type='character',
        dest='OUTPUT_DIR',
        default=NA
    ) %>% 
    add_option(
        c('-b', '--n_bootstraps'),
        help='Total number of bootstraps to generate',
        action='store',
        type='integer',
        default=50000,
        dest='N_BOOTSTRAPS'
    ) %>% 
    add_option(
        c('-c', '--chunk_size'),
        help='How many bootstraps to write per output file',
        action='store',
        type='integer',
        default=200,
        dest='CHUNK_SIZE'
    ) %>% 
    add_option(
        c('-l', '--sample_id_column'),
        help='',
        action='store',
        type='character',
        dest='SAMPLE_ID_VAR',
        default='SampleID'
    ) %>% 
    parse_args(positional_arguments=TRUE)
options <- parsed_args$options
if (is.na(options$OUTPUT_DIR)) {
    stop('Must specify --output_dir')
}
model.files <- parsed_args$args
if (length(model.files) == 0) {
    stop('Must specify .model files as positional arguments')
}
# options=list(SAMPLE_METADATA_FILE=NA, OUTPUT_DIR='/home/sidreed/TalkowskiLab/Projects/XDP.Brains/remote/results/bootstraps', CHUNK_SIZE=15, N_BOOTSTRAPS=100, SAMPLE_ID_VAR='SampleID'); model.files=c("/home/sidreed/TalkowskiLab/Projects/XDP.Brains/remote/models/blank.model","/home/sidreed/TalkowskiLab/Projects/XDP.Brains/remote/models/RS.log.model","/home/sidreed/TalkowskiLab/Projects/XDP.Brains/remote/models/RS.raw.model","/home/sidreed/TalkowskiLab/Projects/XDP.Brains/remote/models/standard.model","/home/sidreed/TalkowskiLab/Projects/XDP.Brains/remote/models/SVA.only.model")
message('finished parsing arguments')
#################
# Dependencies
suppressPackageStartupMessages({
    library(here)
    here::i_am('scripts/pick_bootstrap_samples.R')
    BASE_DIR <- here()
    SCRIPT_DIR <- here('scripts')
    source(file.path(SCRIPT_DIR, 'locations.R'))
    source(file.path(SCRIPT_DIR, 'bootstrap_utils.R'))
    library(tidyverse)
})
################
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
print(glue('Picking bootstraps from {nrow(sample.metadata)} total samples from {options$TISSUES}'))
print(glue('Saving bootstraps files to {options$OUTPUT_DIR}'))
################
# Pick samples for each bootstrap + save lists to files
model.files %>%
    load_model_specifications(
        output_dir=options$OUTPUT_DIR,
        sample.metadata=sample.metadata
    ) %>% 
    add_column(
        chunk_size=options$CHUNK_SIZE,
        n_bootstraps=options$N_BOOTSTRAPS,
        sampleID.variable=options$SAMPLE_ID_VAR
    ) %>% 
    # {.} -> model.specifications
    generate_bootstraps(sample.metadata=sample.metadata)
message('generated lists of sample bootstraps for all models')
