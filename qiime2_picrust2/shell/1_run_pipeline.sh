#!/usr/bin/env bash

export LINE="=========================================="
export _LOG_COUNTER=1

function log {
    printf "\n${LINE}\n\n[$(date '+%d-%m-%Y %H:%M:%S.%N')][Pipeline][OP#$(printf "%02d" ${_LOG_COUNTER})] ${*}\n\n${LINE}\n\n"
    _LOG_COUNTER=$((_LOG_COUNTER + 1))
}


force_curl () {
    while true
    do
      if curl -fsSL "${1}" -o "${2}"
      then
        return
      fi
    done
}


force_docker_pull () {
    while true
    do
      if docker pull "${1}"
      then
        return
      fi
    done
}


get_latest_quay_tag() {
    echo "$(
        curl -fsSL "https://quay.io/api/v1/repository/${1}/${2}" \
        | grep \
            --only-matching \
            --perl-regexp \
            '(?<=\"tags\": {\")[^\"]*(?=\":)'
    )"
}


compose_quay_img() {
    echo "quay.io/${1}/${2}:$(get_latest_quay_tag ${1} ${2})"
}


# Required variables begin
export ROOT_DIR="$(realpath "${ROOT_DIR}")/"
export SAMPLEDATA_DIR="$(realpath "${SAMPLEDATA_DIR}")/"
export SCRIPT_DIR="$(realpath "${SCRIPT_DIR}")/"
# Required variables end

log "Create environment in ${ROOT_DIR}"

export LOG_DIR="${ROOT_DIR}logs/"
export QZV_DIR="${ROOT_DIR}qzv/"
export SAMPLEDATA_CSV="${SAMPLEDATA_DIR}qiime2_sample_data.csv"
export METADATA_TSV="${SAMPLEDATA_DIR}qiime2_meta_data.tsv"

export QIIME2_DIR="${ROOT_DIR}qiime2/"
export QIIME2_FEATURES_BIOM="${QIIME2_DIR}feature-table.biom"
export QIIME2_FEATURES_FASTA="${QIIME2_DIR}dna-sequences.fasta"
export QIIME2_SCRIPT_1="${SCRIPT_DIR}qiime2_1.sh"
export QIIME2_SCRIPT_2="${SCRIPT_DIR}qiime2_2.sh"

export REFERENCE_NAME="SILVA"
export REFERENCE_VERSION="138.1"

export REFERENCE_DIR="/data/reference/${REFERENCE_NAME}/${REFERENCE_NAME}_v${REFERENCE_VERSION}/"
export TAXA_REFERENCE_FEATURES="${REFERENCE_DIR}${REFERENCE_NAME}-${REFERENCE_VERSION}-full-length-seq-taxonomy.qza"
export TAXA_REFERENCE_CLASSIFIER="${REFERENCE_DIR}${REFERENCE_NAME}-${REFERENCE_VERSION}-SSURef-full-length-classifier.qza"
export TAXA_REFERENCE_SEQUENCES="${REFERENCE_DIR}${REFERENCE_NAME}-${REFERENCE_VERSION}-SSURef-Full-Seqs.qza"
export TAXA_REFERENCE_HEADER="${REFERENCE_DIR}${REFERENCE_NAME}_${REFERENCE_VERSION}_taxonomy_headed.tsv"

export PICRUST2_DIR="${ROOT_DIR}picrust2/"
export PICRUST2_SCRIPT="${PICRUST2_DIR}picrust2.sh"

cd "${ROOT_DIR}" || exit 1

mkdir \
    --mode 0777 \
    --parents \
    --verbose \
    "${LOG_DIR}" \
    "${PICRUST2_DIR}" \
    "${QIIME2_DIR}" \
    "${QZV_DIR}" \
    "${SCRIPT_DIR}"



log "Check QIIME2 sampledata"

if [ ! -s "${SAMPLEDATA_CSV}" ] && [ ! -s "${METADATA_TSV}" ]
    then
    log "Create sampledata in ${SAMPLEDATA_DIR}"
    export IMG="ivasilyev/curated_projects:latest"
    force_docker_pull "${IMG}"
    docker run \
        --env RAW_DIR="${RAW_DIR}" \
        --env SAMPLEDATA_DIR="${SAMPLEDATA_DIR}" \
        --net host \
        --rm \
        --volume /data:/data \
        --volume /data1:/data1 \
        --volume /data2:/data2 \
        --volume /data03:/data03 \
        --volume /data04:/data04 \
        "${IMG}" \
            bash -c '
                git pull --quiet;
                python3 ./meta/scripts/qiime2_sample_data.py \
                    --extension ".fastq.gz" \
                    --input "${RAW_DIR}" \
                    --output "${SAMPLEDATA_DIR}"
            ' \
    |& tee "${LOG_DIR}qiime2_sample_data.log"
else
    log "QIIME2 sampledata does exist"
fi



log "Deploy QIIME2 script"

force_curl \
    "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/qiime2_picrust2/shell/2a_run_qiime2_dada2.sh" \
    "${QIIME2_SCRIPT_1}"

force_curl \
    "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/qiime2_picrust2/shell/2b_run_qiime2_vsearch.sh" \
    "${QIIME2_SCRIPT_2}"

cd "${QIIME2_DIR}" || exit 1

log "Run QIIME2"

export IMG="$(compose_quay_img qiime2 core)"

force_docker_pull "${IMG}"

docker run \
    --env QIIME2_DIR="${QIIME2_DIR}" \
    --env QIIME2_SCRIPT_1="${QIIME2_SCRIPT_1}" \
    --env QIIME2_SCRIPT_2="${QIIME2_SCRIPT_2}" \
    --env QIIME2_FEATURES_BIOM="${QIIME2_FEATURES_BIOM}" \
    --env QIIME2_FEATURES_FASTA="${QIIME2_FEATURES_FASTA}" \
    --env QIIME2_OTU_ASV_MAPPER="${QIIME2_OTU_ASV_MAPPER}" \
    --env QIIME2_OTU_TABLE="${QIIME2_OTU_TABLE}" \
    --env SAMPLEDATA_CSV="${SAMPLEDATA_CSV}" \
    --env METADATA_TSV="${METADATA_TSV}" \
    --env TAXA_REFERENCE_FEATURES="${TAXA_REFERENCE_FEATURES}" \
    --env TAXA_REFERENCE_CLASSIFIER="${TAXA_REFERENCE_CLASSIFIER}" \
    --env TAXA_REFERENCE_SEQUENCES="${TAXA_REFERENCE_SEQUENCES}" \
    --env TAXA_REFERENCE_HEADER="${TAXA_REFERENCE_HEADER}" \
    --net host \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    --workdir="${QIIME2_DIR}" \
    "${IMG}" \
    bash -c '
        bash "${QIIME2_SCRIPT_1}" \
        && bash "${QIIME2_SCRIPT_2}"
    ' \
|& tee "${LOG_DIR}qiime2.log"

rm \
    --force \
    --verbose \
    "${QIIME2_SCRIPT_1}" \
    "${QIIME2_SCRIPT_2}"

cd "${ROOT_DIR}" || exit 1



log "Deploy PICRUSt2 script"

force_curl \
    "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/qiime2_picrust2/shell/3_run_picrust2.sh" \
    "${PICRUST2_SCRIPT}"

cd "${PICRUST2_DIR}" || exit 1

log "Run PICRUSt2"

export IMG="$(compose_quay_img biocontainers picrust2)"

export TOTAL_RAM="$(
    free --giga \
    | awk '{ print $2 }' \
    | tail +2 \
    | awk '{ sum += $1 } END { print sum }'
)"

force_docker_pull "${IMG}"

docker run \
    --cpus "$(nproc)" \
    --env QIIME2_FEATURES_BIOM="${QIIME2_FEATURES_BIOM}" \
    --env QIIME2_FEATURES_FASTA="${QIIME2_FEATURES_FASTA}" \
    --env PICRUST2_DIR="${PICRUST2_DIR}" \
    --env PICRUST2_RESULTS_DIR="${PICRUST2_RESULTS_DIR}" \
    --memory "${TOTAL_RAM}g" \
    --memory-swappiness 100 \
    --net host \
    --oom-kill-disable \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    --workdir="${PICRUST2_DIR}" \
    "${IMG}" \
    bash "${PICRUST2_SCRIPT}" \
|& tee "${LOG_DIR}$(basename "${PICRUST2_SCRIPT}").log"

rm -f "${PICRUST2_SCRIPT}"

cd "${ROOT_DIR}" || exit 1

rm -f "$(realpath "${0}")"
