#!/usr/bin/env bash

export ROOT_DIR=""
export RAW_DIR=""

export SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
export SCRIPT_DIR="${ROOT_DIR}scripts/"
export PIPELINE_SCRIPT="${SCRIPT_DIR}1_run_pipeline"

# sudo rm -rf "${ROOT_DIR}"
mkdir \
    --mode 0777 \
    --parents \
    "${ROOT_DIR}" \
    "${SCRIPT_DIR}"

chmod \
    --recursive \
    --verbose \
    0777 \
    "${ROOT_DIR}" \
    "${SCRIPT_DIR}"

cd "${ROOT_DIR}" || exit 1

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/qiime2_picrust2/shell/1_run_pipeline.sh" \
    -o "${PIPELINE_SCRIPT}"

bash "${PIPELINE_SCRIPT}"

rm \
    --force \
    --verbose \
    "$(realpath "${0}")"
