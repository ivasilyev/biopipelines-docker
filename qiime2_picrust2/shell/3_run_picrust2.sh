#!/usr/bin/env bash

export LINE="=========================================="
export _LOG_COUNTER=1


function log {
    # The PICRUSt2 docker image shell does mot have an alias for nanoseconds
    printf "\n${LINE}\n\n[$(date '+%d-%m-%Y %H:%M:%S')][PICRUSt2][OP#$(printf "%02d" ${_LOG_COUNTER})] ${*}\n\n${LINE}\n\n"
    _LOG_COUNTER=$((_LOG_COUNTER + 1))
}


function log_skip {
    printf "\n${LINE}\n\n[$(date '+%d-%m-%Y %H:%M:%S.%N')][Pipeline][OP#$(printf "%02d" ${_LOG_COUNTER})] Skip\n\n${LINE}\n\n"
}


# Required variables start
export PICRUST2_DIR="$(realpath "${PICRUST2_DIR}")/"

export QIIME2_FEATURES_BIOM="$(realpath "${QIIME2_FEATURES_BIOM}")"
export QIIME2_FEATURES_FASTA="$(realpath "${QIIME2_FEATURES_FASTA}")"
# Required variables end

log "Run PATHWAYS in '${PICRUST2_DIR}'"

export LOG_DIR="${PICRUST2_DIR}logs/"
export MAIN_PIPELINE_DIR="${PICRUST2_DIR}main_pipeline/"
export PATHWAY_PIPELINE_DIR="${MAIN_PIPELINE_DIR}pathways_out/"
export TABLES_DIR="${PICRUST2_DIR}described_tables/"
export NPROC="$(grep -c '^processor' "/proc/cpuinfo")"

mkdir \
    --mode 0777 \
    --parents \
    --verbose \
    "${PICRUST2_DIR}" \
    "${LOG_DIR}"

cd "${PICRUST2_DIR}" || exit 1



log "Run the PICRUSt2 pipeline"

export EC_METAGENOMES="${MAIN_PIPELINE_DIR}EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"
export KO_METAGENOMES="${MAIN_PIPELINE_DIR}KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"
export PATHWAYS="${MAIN_PIPELINE_DIR}pathways_out/path_abun_unstrat.tsv.gz"
if [ ! -d "${MAIN_PIPELINE_DIR}" ]
    then
        log "Ensured that the output directory does not exist: '${MAIN_PIPELINE_DIR}'"

        # PICRUSt2 drops too many sequences if the BIOM and the FASTA are from DADA2.
        # After dropping ASV, very few sequences are left.
        # So it is practically better to use the BIOM and the FASTA obtained from
        # `vsearch cluster-features-closed-reference`, all preceded by
        # basic quality-score-based filtering and followed by chimera filtering
        # and aggressive OTU filtering (the treacherous trio, a.k.a. the Bokulich method)
        picrust2_pipeline.py \
            --coverage \
            --hsp_method mp \
            --input "${QIIME2_FEATURES_BIOM}" \
            --processes "${NPROC}" \
            --study_fasta "${QIIME2_FEATURES_FASTA}" \
            --output "${MAIN_PIPELINE_DIR}" \
            --stratified \
            --verbose \
        |& tee "${LOG_DIR}picrust2_pipeline.log"

        log "Run the PICRUSt2 pathway pipeline"

        pathway_pipeline.py \
            --input "${EC_METAGENOMES}" \
            --intermediate "${PATHWAY_PIPELINE_DIR}intermediate" \
            --out_dir "${PATHWAY_PIPELINE_DIR}" \
            --processes "${NPROC}" \
            --verbose \
        |& tee "${LOG_DIR}pathway_pipeline.log"

    else
        log_skip
    fi



log "Convert tables"

mkdir \
    --mode 0777 \
    --parents \
    --verbose \
    "${TABLES_DIR}"

export LEGACY_TSV="${TABLES_DIR}EC_pred_metagenome_contrib_legacy.tsv"
export LEGASY_GZIP="${LEGACY_TSV}.gz"

convert_table.py \
    "${MAIN_PIPELINE_DIR}EC_metagenome_out/pred_metagenome_contrib.tsv.gz" \
    --conversion contrib_to_legacy \
    --output "${LEGASY_GZIP}" \
|& tee "${LOG_DIR}convert_table.log"

# The output file is too large and does not have clear purpose
# zcat "${LEGASY_GZIP}" > "${LEGACY_TSV}"



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



log "Completed running PICRUSt2 in '${PICRUST2_DIR}'"

chmod -R 777 "${PICRUST2_DIR}"

cd ..

rm -f "$(realpath "${0}")"

exit 0
