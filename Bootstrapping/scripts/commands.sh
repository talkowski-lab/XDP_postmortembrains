# Fixed Locations
ROOT_DIR="/data/talkowski/Samples/XDP/Brain_RNASeq/v2/bootstrapDE"
SCRIPT_DIR="${ROOT_DIR}/scripts"
MODELS_DIR="${ROOT_DIR}/models"
RESULTS_DIR="${ROOT_DIR}/results/bootstraps"
MODEL_DIR="${RESULTS_DIR}/${MODEL_NAME}"
BOOTSTRAP_LISTS_DIR="${MODEL_DIR}/bootstrap.sample.lists"
DEG_RESULTS_DIR="${MODELS_DIR}/DEG.results"
# Input Data
SAMPLE_METADATA_FILE="${ROOT_DIR}/CAU.sample_metadata.tsv"
COUNTS_FILE="${ROOT_DIR}/CAU.counts.tsv"
# Model Specific Output Directories
MODEL_LOG_DIR="${OUTPUT_DIR}/logs"
# Pick samples for each bootstrap
pick_bootstrap_sample_sets() {
    # Rscript scripts/pick_bootstrap_samples.R -o ./results/bootstraps/ -b 50000 -c 200 ./models/*.model
    Rscript ${SCRIPT_DIR}/pick_bootstrap_samples.R \
        --output_dir      ${RESULTS_DIR}           \
        --n_bootstraps    50000                    \
        --chunk_size      200                      \
        ${MODELS_DIR}/*.model
}
# Launch 1 job locally for specific bootstrap file
launch_job_local() {
    bootstraps_chunk_file="$1"
    ${SCRIPT_DIR}/run_bootstraps.sh           \
        --results_dir       ${RESULTS_DIR}    \
        --log_dir           ${MODEL_LOG_DIR}  \
        ${bootstraps_chunk_file}
        
}
# Launch 1 job per set of bootstraps 
launch_all_jobs_slurm() {
    # Launch 1 job per bootstrap file in bootstraps_dir, each job computes in parallel 
    ${SCRIPT_DIR}/run_bootstraps.sh          \
        --results_dir       ${RESULTS_DIR}   \
        --log_dir           ${MODEL_LOG_DIR} \
        --r_libs            "${R_LIBS}"      \
        # --ntasks-per-node 4 \
        --cpus-per-task 4
        --mem 30000                          \
        --partition 'normal'                 \
        ${RESULTS_DIR}/results/bootstraps/${MODELS_NAME}/bootstrap.sample.lists/
}

