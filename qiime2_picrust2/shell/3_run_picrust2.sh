#!/usr/bin/env bash

export LINE="======================================================================================"

function log {
    # The PICRUSt2 docker image shell does mot have an alias for nanoseconds
    printf "\n${LINE}\n\n[$(date '+%d-%m-%Y %H:%M:%S')][PICRUSt2] $@\n\n${LINE}\n\n"
}

# Required variables begin
export PICRUST2_DIR="$(realpath "${PICRUST2_DIR}")/"
export QIME2_FEATURES_BIOM="$(realpath "${QIME2_FEATURES_BIOM}")"
export QIME2_FEATURES_FASTA="$(realpath "${QIME2_FEATURES_FASTA}")"
# Required variables end

log "Run PICRUSt2 in ${PICRUST2_DIR}"

export LOG_DIR="${PICRUST2_DIR}logs/"
export PIPELINE_DIR="${PICRUST2_DIR}main_pipeline/"
export TABLES_DIR="${PICRUST2_DIR}described_tables/"
export NPROC="$(grep -c '^processor' "/proc/cpuinfo")"

mkdir -p "${PICRUST2_DIR}" "${LOG_DIR}"
cd "${PICRUST2_DIR}" || exit 1
log Ensure that the output directory does not exist
rm -rf "${PIPELINE_DIR}"

log Run the PICRUSt2 pipeline
picrust2_pipeline.py \
    --coverage \
    --hsp_method mp \
    --input "${QIME2_FEATURES_BIOM}" \
    --processes "${NPROC}" \
    --study_fasta "${QIME2_FEATURES_FASTA}" \
    --output "${PIPELINE_DIR}" \
    --stratified \
    --verbose \
    |& tee "${LOG_DIR}picrust2_pipeline.log"

log Run the PICRUSt2 pathway pipeline
pathway_pipeline.py \
    --input "${PIPELINE_DIR}EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
    --intermediate "${PIPELINE_DIR}pathways_out/intermediate" \
    --out_dir "${PIPELINE_DIR}pathways_out" \
    --processes "${NPROC}" \
    --verbose \
    |& tee "${LOG_DIR}pathway_pipeline.log"

log Convert tables
mkdir -p "${TABLES_DIR}EC_metagenome_out/"
convert_table.py \
    "${PIPELINE_DIR}EC_metagenome_out/pred_metagenome_contrib.tsv.gz" \
    --conversion contrib_to_legacy \
    --output "${TABLES_DIR}EC_metagenome_out/pred_metagenome_contrib.legacy.tsv" \
    |& tee "${LOG_DIR}convert_table.log"

log Add KEGG ENZYME descriptions
add_descriptions.py \
    --input "${PIPELINE_DIR}EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
    --map_type EC \
    --output "${TABLES_DIR}EC_metagenome_out/pred_metagenome_unstrat_described.tsv" \
    |& tee "${LOG_DIR}add_descriptions-EC.log"

log Add KEGG ORTHOLOGY descriptions
add_descriptions.py \
    --input "${PIPELINE_DIR}KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
    --map_type KO \
    --output "${TABLES_DIR}KO_metagenome_out/pred_metagenome_unstrat_described.tsv" \
    |& tee "${LOG_DIR}add_descriptions-KO.log"

log Add MetaCyc descriptions
add_descriptions.py  \
    --input "${PIPELINE_DIR}pathways_out/path_abun_unstrat.tsv.gz" \
    --map_type METACYC \
    --output "${TABLES_DIR}pathways_out/path_abun_unstrat_described.tsv" \
    |& tee "${LOG_DIR}add_descriptions-METACYC.log"

log "Completed running PICRUSt2 in ${PICRUST2_DIR}"
chmod -R 777 "$(pwd)"
cd ..
exit 0
