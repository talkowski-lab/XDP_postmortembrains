# Generate a bunch of tsv files, each line in each tsv corresponds to a pair
# of control vs case samples
################
# Parse Arguments
library(magrittr)
library(optparse)
parsed_args <- 
    OptionParser() %>%
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
        c('-m', '--model_file'),
        help='.model file specifying model parameters',
        action='store',
        type='character',
        dest='MODEL_FILE',
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
        c('-t', '--threads'),
        help='How many threads to use',
        action='store',
        type='integer',
        dest='THREADS',
        default=NA
    ) %>% 
    # add_option(
    #     c('-l', '--log_dir'),
    #     help='dir to save log files to',
    #     action='store',
    #     type='character',
    #     dest='LOG_DIR',
    #     default=NA
    # ) %>% 
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
chunk.files <- parsed_args$args
# options=list(FORCE_REDO=TRUE, SAMPLE_METADATA_FILE=NA, COUNTS_FILE=NA, MODEL_FILE='/home/sidreed/TalkowskiLab/Projects/XDP.Brains/remote/results/bootstraps/standard/this.model', GENE_ID_COLUMN='EnsemblID', OUTPUT_DIR='/home/sidreed/TalkowskiLab/Projects/XDP.Brains/remote/results/bootstraps/standard/DEG.results', SAMPLE_ID_VAR='SampleID'); chunk.files=c("/home/sidreed/TalkowskiLab/Projects/XDP.Brains/remote/results/bootstraps/standard/bootstrap.samples.lists/chunk1.tsv")
################
# Validate stuff
if (is.na(options$MODEL_FILE)) {
    stop("Must specify --model_file path")
}
if (is.na(options$OUTPUT_DIR)) {
    stop("Must specify --output_dir path")
}
# if (is.na(options$LOG_DIR)) {
#     stop("Must specify --log_dir path")
# }
if (length(chunk.files) == 0){
    stop('Must specify chunk of bootstraps to run')
}
# options=list(); chunk.files=list()
#################
# Dependencies
suppressPackageStartupMessages({
    library(here)
    here::i_am('scripts/run_bootstraps.R')
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
# Load counts
counts.matrix <- 
    if (is.na(options$COUNTS_FILE)) {
        load_counts()  # in bootstrap.utils.R
    } else {
        options$COUNTS_FILE %>% 
        read_tsv(show_col_types=FALSE)
    }
# Load model params
# Should be 1 row tibble (1 model), each column is a model parameter
model.parameters <- 
    c(options$MODEL_FILE) %>% 
    load_model_specifications(
        output_dir="",
        sample.metadata=sample.metadata
    ) %>%
    mutate(results.dir=options$OUTPUT_DIR)
# List all chunks of bootstraps to be processed (usually just one)
chunks.list <- 
    chunk.files %>%
    load_bootstrap_chunks()
################
# Compute DEG results in parallel on all bootstraps for each chunk (in serial)
library(furrr)
workers <- 
    if (is.na(options$THREADS)) {
        length(availableWorkers())
    } else {
        options$THREADS
    }
message(glue('Using {workers} threads'))
plan(multicore, workers=workers)
chunks.list %>%
    pmap(
        run_pipeline_on_chunk,
        model.parameters=model.parameters,
        sample.metadata=sample.metadata,
        counts.matrix=counts.matrix,
        sampleID.variable=options$SAMPLE_ID_VAR,
        gene_id_column=options$GENE_ID_COLUMN,
        force_redo=options$FORCE_REDO,
        .progress=TRUE
    )
