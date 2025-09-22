#!/bin/bash
# Author: Siddharth Reed, Talkowski Lab
# Default CLI Options
# Locations
BASE_DIR="$(pwd)"
SCRIPT_DIR="${BASE_DIR}/scripts"
# model params
R_LIBS_DIR="${R_LIBS}"
FORCE_REDO=0
# Job submission params
PARALELL=0
THREADS=8
MEMORY='30G'
QUEUE='normal'
# Args
help() {
    echo "Compute DEG results with DESeq2 for all the bootstrapped sets of samples specified in \$bootstraps_dir
Usage: $0 [OPTIONS]
        -p|--parallel
            use lsf or dont to parallelize jobs
        -q|--queue
            lsf queue to submit jobs to
        -m|--memory
            memory for each lsf job
        -n|--workers
            workers per lsf job
        -h     
            Print this message
"
}

run_chunk() {
    chunk_file="${1}"
    chunk_id="$(basename ${chunk_file%%.tsv})"
    echo ${chunk_id}
    Rscript ${SCRIPT_DIR}/run_bootstraps.R   \
        --threads ${THREADS}                 \
        --model_file ${MODEL_DIR}/this.model \
        --output_dir ${DEG_RESULTS_DIR}      \
        "${chunk_file}" |& tee ${LOG_DIR}/${chunk_id}.log
}

make_cmd() {
    chunk_file="${1}"
    chunk_id="$(basename ${chunk_file%%.tsv})"
    echo ${chunk_id}
    echo "module load R/4.4.0
export R_LIBS=${R_LIBS_DIR}
Rscript ${SCRIPT_DIR}/run_bootstraps.R
    --threads ${THREADS}
    --log_file   ${LOG_DIR}/${chunk_id}.log
    --model_file ${MODEL_DIR}/this.model
    --output_dir ${DEG_RESULTS_DIR}
    ${chunk_file}"
}

launch_job() {
    chunk_file="${1}"
    chunk_id="$(basename ${chunk_file%%.tsv})"
    # echo "sbatch -o ${LOG_DIR}/${chunk_id}.out -e ${LOG_DIR}/${chunk_id}.err -J Bootstrapping_${chunk_id} --partition ${QUEUE} --cpus-per-task ${WORKERS} --mem ${MEMORY} --wrap=$(make_cmd ${chunk_file})"
    sbatch \
        -o "${LOG_DIR}/${chunk_id}.out" \
        -e "${LOG_DIR}/${chunk_id}.err" \
        -J "Bootstrapping_${chunk_id}" \
        --partition "${QUEUE}" \
        --cpus-per-task "${THREADS}" \
        --mem "${MEMORY}" \
        --wrap="$(make_cmd ${chunk_file})"
}
# Handle CLI arg
[[ $? -ne 0 ]] && echo "No Args" && exit 1
VALID_ARGS=$(getopt -o r:pq:m:t:h --long rlibs:,parallel,partition:,memory:,threads:help -- "$@")
eval set -- "$VALID_ARGS"
while [ : ]; do
    case "$1" in
        -r|--rlibs)
            R_LIBS_DIR="${2}"
            shift 2
        ;;
        -p|--parallel)
            PARALELL=1
            shift 1
        ;;
        -q|--partition)
            QUEUE="${2}"
            shift 2
        ;;
        -m|--memory)
            MEMORY="${2}"
            shift 2
        ;;
        -t|--threads)
            THREADS="${2}"
            shift 2
        ;;
        -f|--force_redo)
            FORCE_REDO=1
            shift 1
        ;;
        -h|--help) 
            help 
            exit 0 
            ;;
        --)
            shift 
            break
            ;;
    esac
done
# Set up paths
MODEL_DIR="${1}"
BOOTSTRAP_LIST_DIR="${MODEL_DIR}/bootstrap.samples.lists"
MODEL_FILE="${MODEL_DIR}/this.model"
LOG_DIR="${MODEL_DIR}/logs"
mkdir -p "${LOG_DIR}"
DEG_RESULTS_DIR="${MODEL_DIR}/DEG.results"
mkdir -p "${DEG_RESULTS_DIR}"
# 
for chunk_file in ${BOOTSTRAP_LIST_DIR}/chunk*.tsv; do
    chunk_file="$(readlink -e "${chunk_file}")"
    chunk_id="$(basename ${chunk_file%%.tsv})"
    chunk_result_file="${DEG_RESULTS_DIR}/Genotype-XDP-Control/${chunk_id}.tsv"
    if [ -s  ${chunk_result_file} ] && [ ${FORCE_REDO} == 0 ]; then
        echo ${chunk_id} result file detected, skipping
    else
        [[ ${PARALELL} == 1 ]] && launch_job "${chunk_file}" || run_chunk "${chunk_file}"
    fi
    # if [[ PARALELL == 1 ]]; then
    #     launch_job "${chunk_file}"
    # else
    #     run_chunk "${chunk_file}"
    # fi
done
