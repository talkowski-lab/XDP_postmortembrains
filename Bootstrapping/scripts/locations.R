# Locations
# BASE_DIR <- "/data/talkowski/Samples/XDP/Brain_RNASeq/v2/bootstrapDE"
SCRIPT_DIR                <- file.path(BASE_DIR, "scripts")
BOOTSTRAP_DIR             <- file.path(BASE_DIR, "results", "bootstraps")
MODEL_SPECIFICATIONS_DIR  <-  file.path(BASE_DIR, "models")
# InputData
SAMPLE_METADATA_FILE      <- file.path(BASE_DIR, "All.sample_metadata.tsv")
RAW_COUNTS_FILE           <- file.path(BASE_DIR, "All.counts.tsv")
# SAMPLE_METADATA_FILE      <- file.path(BASE_DIR, "CAU.sample_metadata.tsv")
# RAW_COUNTS_FILE           <- file.path(BASE_DIR, "CAU.counts.tsv")
MODEL_SPECIFICATIONS_FILE <- file.path(BASE_DIR, "models.tsv")
BOOTSTRAP_MODEL_SPECIFICATIONS <- file.path(BOOTSTRAP_DIR, 'bootstrap.model.specifications.tsv')
