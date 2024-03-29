#!/usr/bin/env bash

export LINE="=========================================="
export _LOG_COUNTER=1

function log {
    printf "\n${LINE}\n\n[$(date '+%d-%m-%Y %H:%M:%S.%N')][Pipeline][OP#${_LOG_COUNTER}] $@\n\n${LINE}\n\n"
    _LOG_COUNTER=$((_LOG_COUNTER + 1))
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
export SAMPLEDATA_CSV="${SAMPLEDATA_DIR}qiime2_sample_data.csv"
export METADATA_TSV="${SAMPLEDATA_DIR}qiime2_meta_data.tsv"

export QIIME2_DIR="${ROOT_DIR}qiime2/"
export QIME2_FEATURES_BIOM="${QIIME2_DIR}bioms/feature-table.biom"
export QIME2_FEATURES_FASTA="${QIIME2_DIR}closed_references/dna-sequences.fasta"
export QIIME2_SCRIPT="${QIIME2_DIR}qiime2.sh"

export REFERENCE_NAME="SILVA"
export REFERENCE_VERSION="138.1"

export REFERENCE_DIR="/data/reference/${REFERENCE_NAME}/${REFERENCE_NAME}_v${REFERENCE_VERSION}/"
export TAXA_REFERENCE_FEATURES="${REFERENCE_DIR}${REFERENCE_NAME}-${REFERENCE_VERSION}-full-length-seq-taxonomy.qza"
export TAXA_REFERENCE_CLASSIFIER="${REFERENCE_DIR}${REFERENCE_NAME}-${REFERENCE_VERSION}-SSURef-full-length-classifier.qza"
export TAXA_REFERENCE_SEQUENCES="${REFERENCE_DIR}${REFERENCE_NAME}-${REFERENCE_VERSION}-SSURef-Full-Seqs.qza"
export TAXA_REFERENCE_HEADER="${REFERENCE_DIR}${REFERENCE_NAME}_${REFERENCE_VERSION}_Taxonomy_headed.tsv"

export PICRUST2_DIR="${ROOT_DIR}picrust2/"
export PICRUST2_SCRIPT="${PICRUST2_DIR}picrust2.sh"
export RESULT_DIR="${ROOT_DIR}results/"

export OTU_TABLE="${RESULT_DIR}OTUs_with_taxa.tsv"

cd "${ROOT_DIR}" || exit 1

mkdir \
    --mode 0777 \
    --parents \
    --verbose \
    "${QIIME2_DIR}" \
    "${PICRUST2_DIR}" \
    "${SCRIPT_DIR}" \
    "${LOG_DIR}"



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

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/qiime2_picrust2/shell/2_run_qiime2_dada2.sh" \
    -o "${QIIME2_SCRIPT}"
cd "${QIIME2_DIR}" || exit 1

log "Run QIIME2"

export IMG="$(compose_quay_img qiime2 core)"

force_docker_pull "${IMG}"

docker run \
    --env QIIME2_DIR="${QIIME2_DIR}" \
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
    bash "${QIIME2_SCRIPT}" \
|& tee "${LOGS_DIR}$(basename "${QIIME2_SCRIPT}").log"

rm -f "${QIIME2_SCRIPT}"

cd "${ROOT_DIR}" || exit 1



log "Deploy PICRUSt2 script"

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/qiime2_picrust2/shell/3_run_picrust2.sh" \
    -o "${PICRUST2_SCRIPT}"

cd "${PICRUST2_DIR}" || exit 1

log "Run PICRUSt2"

export IMG="$(compose_quay_img biocontainers picrust2)"

force_docker_pull "${IMG}"

docker run \
    --env QIME2_FEATURES_BIOM="${QIME2_FEATURES_BIOM}" \
    --env QIME2_FEATURES_FASTA="${QIME2_FEATURES_FASTA}" \
    --env PICRUST2_DIR="${PICRUST2_DIR}" \
    --net host \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    --workdir="${PICRUST2_DIR}" \
    "${IMG}" \
    bash "${PICRUST2_SCRIPT}" \
|& tee "${LOGS_DIR}$(basename "${PICRUST2_SCRIPT}").log"

rm -f "${PICRUST2_SCRIPT}"

cd "${ROOT_DIR}" || exit 1



log "Copy QIIME2 and PICRUSt2 pipeline output tables"

mkdir -p "${RESULT_DIR}"

find "${ROOT_DIR}" \
    -type f \( \
        -name "path_abun_unstrat_described.tsv" \
        -o -name "pred_metagenome_contrib.legacy.tsv" \
        -o -name "OTUs_with_taxa.tsv" \
    \) -print0 \
| xargs \
    -0 \
    --max-procs "$(nproc)" \
    -I "{}" \
        bash -c '
            export SRC_FILE="{}";
            export DST_FILE="${RESULT_DIR}$(basename "${SRC_FILE}")";
            echo "Copy \"${SRC_FILE}\" to \"${DST_FILE}\"";
            cp \
                "${SRC_FILE}" \
                "${DST_FILE}";
        '

find "${ROOT_DIR}" \
    -type f \
    -name "pred_metagenome_unstrat_described.tsv" \
    -print0 \
| xargs \
    -0 \
    --max-procs "$(nproc)" \
    -I "{}" \
        bash -c '
            export SRC_FILE="{}";
            export DST_FILE="${RESULT_DIR}$(basename "$(dirname "${SRC_FILE}")")_$(basename "${SRC_FILE}")";
            echo "Copy \"${SRC_FILE}\" to \"${DST_FILE}\"";
            cp \
                "${SRC_FILE}" \
                "${DST_FILE}";
        '



# The first line of the raw file is '# Constructed from biom file'
sed -i '1d' "${OTU_TABLE}"

log "Concatenate tables"

export IMG="ivasilyev/curated_projects:latest"

force_docker_pull "${IMG}"

docker run \
    --env OTU_TABLE="${OTU_TABLE}" \
    --env TAXA_REFERENCE_HEADER="${TAXA_REFERENCE_HEADER}" \
    --net host \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data2:/data2 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    "${IMG}" \
        bash -c '
            OUT_FILE="${OTU_TABLE%.*}_annotated.tsv";
            echo Concatenate table \"${OTU_TABLE}\" and \"${TAXA_REFERENCE_HEADER}\" into \"${OUT_FILE}\";
            git pull --quiet && \
            python3 ./meta/scripts/concatenate_tables.py \
                --axis 1 \
                --index "#OTU ID" \
                --input \
                    "${TAXA_REFERENCE_HEADER}" \
                    "${OTU_TABLE}" \
                --output "${OUT_FILE}"
        ' \
|& tee "${LOGS_DIR}concatenate_tables.log"



log "All pipeline runs ended"
