#!/usr/bin/env bash

export LINE="=========================================="
export _LOG_COUNTER=1


function log {
    # The PICRUSt2 docker image shell does mot have an alias for nanoseconds
    printf "\n${LINE}\n\n[$(date '+%d-%m-%Y %H:%M:%S')][PICRUSt2][OP#$(printf "%02d" ${_LOG_COUNTER})] ${*}\n\n${LINE}\n\n"
    _LOG_COUNTER=$((_LOG_COUNTER + 1))
}


# Required variables start
export PICRUST2_DIR="$(realpath "${PICRUST2_DIR}")/"
# Will be symlinked
export PICRUST2_RESULTS_DIR="$(realpath "${PICRUST2_DIR}")"

export QIIME2_FEATURES_BIOM="$(realpath "${QIIME2_FEATURES_BIOM}")"
export QIIME2_FEATURES_FASTA="$(realpath "${QIIME2_FEATURES_FASTA}")"
# Required variables end

log "Run PATHWAYS in '${PICRUST2_DIR}'"

export LOG_DIR="${PICRUST2_DIR}logs/"
export PIPELINE_DIR="${PICRUST2_DIR}main_pipeline/"
export TABLES_DIR="${PICRUST2_DIR}described_tables/"
export NPROC="$(grep -c '^processor' "/proc/cpuinfo")"

mkdir \
    --mode 0777 \
    --parents \
    --verbose \
    "${PICRUST2_DIR}" \
    "${LOG_DIR}"

cd "${PICRUST2_DIR}" || exit 1
log "Ensure that the output directory does not exist"
rm -rf "${PIPELINE_DIR}"



log "Run the PICRUSt2 pipeline"

export EC_METAGENOMES="${PIPELINE_DIR}EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"
export KO_METAGENOMES="${PIPELINE_DIR}KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"
export PATHWAYS="${PIPELINE_DIR}pathways_out/path_abun_unstrat.tsv.gz"

picrust2_pipeline.py \
    --coverage \
    --hsp_method mp \
    --input "${QIIME2_FEATURES_BIOM}" \
    --processes "${NPROC}" \
    --study_fasta "${QIIME2_FEATURES_FASTA}" \
    --output "${PIPELINE_DIR}" \
    --stratified \
    --verbose \
|& tee "${LOG_DIR}picrust2_pipeline.log"



log "Run the PICRUSt2 pathway pipeline"

pathway_pipeline.py \
    --input "${EC_METAGENOMES}" \
    --intermediate "${PIPELINE_DIR}pathways_out/intermediate" \
    --out_dir "${PIPELINE_DIR}pathways_out" \
    --processes "${NPROC}" \
    --verbose \
|& tee "${LOG_DIR}pathway_pipeline.log"



log "Convert tables"

mkdir -p "${TABLES_DIR}"

export LEGASY_TSV="${TABLES_DIR}EC_pred_metagenome_contrib_legacy.tsv"

convert_table.py \
    "${PIPELINE_DIR}EC_metagenome_out/pred_metagenome_contrib.tsv.gz" \
    --conversion contrib_to_legacy \
    --output "${LEGASY_TSV}.gz" \
|& tee "${LOG_DIR}convert_table.log"

gzip "${LEGASY_TSV}.gz" \
> "${LEGASY_TSV}"



log "Add KEGG ENZYME descriptions"

add_descriptions.py \
    --input "${EC_METAGENOMES}" \
    --map_type EC \
    --output "${TABLES_DIR}EC_pred_metagenome_unstrat_described.tsv" \
|& tee "${LOG_DIR}add_descriptions-EC.log"



log "Add KEGG ORTHOLOGY descriptions"

add_descriptions.py \
    --input "${KO_METAGENOMES}" \
    --map_type KO \
    --output "${TABLES_DIR}KO_pred_metagenome_unstrat_described.tsv" \
|& tee "${LOG_DIR}add_descriptions-KO.log"



log "Add MetaCyc descriptions"

add_descriptions.py  \
    --input "${PATHWAYS}" \
    --map_type METACYC \
    --output "${TABLES_DIR}pathways_abun_unstrat_described.tsv" \
|& tee "${LOG_DIR}add_descriptions-METACYC.log"



log "Export denormalized frequencies to use in report"

ln \
    --symbolic \
    --verbose \
    "${TABLES_DIR}" \
    "${PICRUST2_RESULTS_DIR}"



log "Completed running PICRUSt2 in '${PICRUST2_DIR}'"

chmod -R 777 "$(pwd)"

cd ..

rm -f "$(realpath "${0}")"

exit 0
