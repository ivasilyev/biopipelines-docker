#!/usr/bin/env bash

export LINE="=========================================="
export _LOG_COUNTER=1

function log {
    printf "\n${LINE}\n\n[$(date '+%d-%m-%Y %H:%M:%S.%N')][QIIME2][OP#$(printf "%02d" ${_LOG_COUNTER})] ${*}\n\n${LINE}\n\n"
    _LOG_COUNTER=$((_LOG_COUNTER + 1))
}

function md {
    for i in "${@}"
        do
        mkdir \
            --mode 0777 \
            --parents \
            --verbose \
            "$(dirname "${i}")"
        done
}

# Required variables begin
export QIIME2_DIR="$(realpath "${QIIME2_DIR}")/"
export QIIME2_FEATURES_BIOM="$(realpath "${QIIME2_FEATURES_BIOM}")"
export QIIME2_FEATURES_FASTA="$(realpath "${QIIME2_FEATURES_FASTA}")"

export SAMPLEDATA_CSV="$(realpath "${SAMPLEDATA_CSV}")"
export METADATA_TSV="$(realpath "${METADATA_TSV}")"

export TAXA_REFERENCE_FEATURES="$(realpath "${TAXA_REFERENCE_FEATURES}")"
export TAXA_REFERENCE_CLASSIFIER="$(realpath "${TAXA_REFERENCE_CLASSIFIER}")"
export TAXA_REFERENCE_SEQUENCES="$(realpath "${TAXA_REFERENCE_SEQUENCES}")"
export TAXA_REFERENCE_HEADER="$(realpath "${TAXA_REFERENCE_HEADER}")"
# Required variables end

log "Run QIIME2 in ${QIIME2_DIR}"

export TOOL_NAME="vsearch"
export TOOL_DIR="${QIIME2_DIR}${TOOL_NAME}/"
export LOG_DIR="${TOOL_DIR}logs/"
export CONSENSUS_THRESHOLD=97
export GROUPING_COLUMN_NAME="SubjectID"
export PREV_CONTROL_COLUMN="Subgroup"
export PREV_CONTROL_INDICATOR="ControlNegative"
export TAXA_COLUMN_NAME="taxonomy"
export NPROC="$(grep -c '^processor' "/proc/cpuinfo")"



mkdir -p "${LOG_DIR}"
cd "${QIIME2_DIR}" || exit 1

qiime dev refresh-cache



log "Import and convert pre-demultiplexed paired-end FASTQ files to QIIME2 artifact"

export DEMULTIPLEXED_DIR="${QIIME2_DIR}demultiplexed_reads/"
export DEMULTIPLEXED_READS="${DEMULTIPLEXED_DIR}demultiplexed_PE_reads.qza"

md "${DEMULTIPLEXED_READS}"

if [[ ! -s "${DEMULTIPLEXED_READS}" ]]
    then

    qiime tools import \
        --input-format PairedEndFastqManifestPhred33 \
        --input-path "${SAMPLEDATA_CSV}" \
        --output-path "${DEMULTIPLEXED_READS}" \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
    |& tee "${LOG_DIR}tools import.log"

    log "Summarize sequences"

    qiime demux summarize \
        --i-data "${DEMULTIPLEXED_READS}" \
        --o-visualization "${DEMULTIPLEXED_DIR}demultiplexed_PE_reads.qzv" \
         --verbose \
    |& tee "${LOG_DIR}demux summarize demux_PE_reads.log"

    else
        echo "Skip"
    fi



log "Merge demultiplexed paired-end reads"

export MERGED_READS="${TOOL_DIR}merged_reads/merged_PE_reads.qza"

md "${MERGED_READS}"

qiime vsearch merge-pairs \
    --i-demultiplexed-seqs "${DEMULTIPLEXED_READS}" \
    --o-merged-sequences "${MERGED_READS}" \
    --p-allowmergestagger \
    --p-threads 8 \
    --verbose



log "Filter merged sequences based on Q scores"

export QUALITY_FILTER_DIR="${TOOL_DIR}quality_filter/"
export QUALITY_FILTERED_SEQUENCES="${QUALITY_FILTER_DIR}q_scoring_sequences.qza"

md "${QUALITY_FILTERED_SEQUENCES}"

if [[ ! -s "${QUALITY_FILTERED_SEQUENCES}" ]]
    then

    md "${QUALITY_FILTERED_SEQUENCES}"

    qiime quality-filter q-score \
        --i-demux "${MERGED_READS}" \
        --o-filtered-sequences "${QUALITY_FILTERED_SEQUENCES}" \
        --o-filter-stats "${QUALITY_FILTER_DIR}q_scoring_statistics.qza" \
        --verbose \
    |& tee "${LOG_DIR}quality-filter q-score.log"

    else
        echo "Skip"
    fi



log "Dereplicate sequences"

export DEREPLICATED_DIR="${TOOL_DIR}dereplicated/"
export DEREPLICATED_SEQUENCES="${DEREPLICATED_DIR}dereplicated_sequences.qza"
export DEREPLICATED_FREQUENCIES="${DEREPLICATED_DIR}dereplicated_frequency_table.qza"

md "${DEREPLICATED_SEQUENCES}"

if [[ ! -s "${DEREPLICATED_SEQUENCES}" ]]
    then

    md "${DEREPLICATED_SEQUENCES}"

    qiime vsearch dereplicate-sequences \
        --i-sequences "${QUALITY_FILTERED_SEQUENCES}" \
        --o-dereplicated-sequences "${DEREPLICATED_SEQUENCES}" \
        --o-dereplicated-table "${DEREPLICATED_FREQUENCIES}" \
        --verbose \
    |& tee "${LOG_DIR}vsearch dereplicate-sequences.log"

    else
        echo "Skip"
    fi



export CLUSTERED_DIR="${TOOL_DIR}cluster_features/"
export CLUSTERED_SEQUENCES="${CLUSTERED_DIR}closed_reference_clustered_sequences.qza"
export CLUSTERED_FREQUENCIES="${CLUSTERED_DIR}closed_reference_clustered_table.qza"

if [[ ! -s "${CLUSTERED_SEQUENCES}" ]]
    then

    log "Cluster closed references at ${CONSENSUS_THRESHOLD} percent"

    md "${CLUSTERED_SEQUENCES}"

    qiime vsearch cluster-features-closed-reference \
        --i-reference-sequences "${TAXA_REFERENCE_SEQUENCES}" \
        --i-sequences "${DEREPLICATED_SEQUENCES}" \
        --i-table "${DEREPLICATED_FREQUENCIES}" \
        --o-clustered-sequences "${CLUSTERED_SEQUENCES}" \
        --o-clustered-table "${CLUSTERED_FREQUENCIES}" \
        --o-unmatched-sequences "${CLUSTERED_DIR}closed_reference_unmatched_sequences.qza" \
        --p-perc-identity 0.${CONSENSUS_THRESHOLD} \
        --p-threads "${NPROC}" \
        --verbose \
    |& tee "${LOG_DIR}vsearch cluster-features-closed-reference.log"

    else
        echo "Skip"
    fi



export DECHIMERIZATION_DIR="${TOOL_DIR}uchime-denovo/"
export CHIMERIC_SEQUENCES="${DECHIMERIZATION_DIR}chimeras.qza"
export CHIMERA_FILTERING_SEQUENCES="${DECHIMERIZATION_DIR}nonchimeras.qza"
export DECHIMERIZATION_STATS="${DECHIMERIZATION_DIR}statistics.qza"

log "Run de novo chimera checking"

qiime vsearch uchime-denovo \
    --i-sequences "${CLUSTERED_SEQUENCES}" \
    --i-table "${CLUSTERED_FREQUENCIES}" \
    --o-chimeras "${CHIMERIC_SEQUENCES}" \
    --o-nonchimeras "${CHIMERA_FILTERING_SEQUENCES}" \
    --o-stats "${DECHIMERIZATION_STATS}" \
    --verbose



log "Visualize chimera check summary"

qiime metadata tabulate \
    --m-input-file "${DECHIMERIZATION_STATS}" \
    --o-visualization "${DECHIMERIZATION_DIR}statistics.qzv"



log "Exclude chimeras and borderline chimeras from feature table"

export NON_CHIMERIC_FREQUENCIES="${DECHIMERIZATION_DIR}table_nonchimeric.qza"

qiime feature-table filter-features \
    --i-table "${CLUSTERED_FREQUENCIES}" \
    --m-metadata-file "${CHIMERA_FILTERING_SEQUENCES}" \
    --o-filtered-table "${NON_CHIMERIC_FREQUENCIES}" \
    --p-no-exclude-ids \
    --verbose

qiime feature-table summarize \
    --i-table "${NON_CHIMERIC_FREQUENCIES}" \
    --o-visualization "${DECHIMERIZATION_DIR}table-nonchimeric.qzv"



log "Exclude chimeras and borderline chimeras from feature sequences"

export NON_CHIMERIC_SEQUENCES="${DECHIMERIZATION_DIR}representative_sequences_nonchimeric.qza"

qiime feature-table filter-seqs \
    --i-data "${CLUSTERED_SEQUENCES}" \
    --m-metadata-file "${CHIMERA_FILTERING_SEQUENCES}" \
    --o-filtered-data "${NON_CHIMERIC_SEQUENCES}" \
    --p-no-exclude-ids \
    --verbose



log "Exclude chimeras but retain borderline chimeras from feature table"

export BORDERLINE_CHIMERIC_FREQUENCIES="${DECHIMERIZATION_DIR}table_borderline_chimeric.qza"

qiime feature-table filter-features \
    --i-table "${CLUSTERED_FREQUENCIES}" \
    --m-metadata-file "${CHIMERIC_SEQUENCES}" \
    --o-filtered-table "${BORDERLINE_CHIMERIC_FREQUENCIES}" \
    --p-exclude-ids \
    --verbose

qiime feature-table summarize \
    --i-table "${BORDERLINE_CHIMERIC_FREQUENCIES}" \
    --o-visualization "${DECHIMERIZATION_DIR}table_borderline_chimeric.qzv"



log "Exclude chimeras but retain borderline chimeras from feature sequences"

export BORDERLINE_CHIMERIC_SEQUENCES="${DECHIMERIZATION_DIR}representative_sequences_borderline_chimeric.qza"

qiime feature-table filter-seqs \
    --i-data "${CLUSTERED_SEQUENCES}" \
    --m-metadata-file "${CHIMERIC_SEQUENCES}" \
    --o-filtered-data "${BORDERLINE_CHIMERIC_SEQUENCES}" \
    --p-exclude-ids \
    --verbose



log "Export OTU"

export BIOM_DIR="${TOOL_DIR}bioms/"
export BIOM_RAW="${BIOM_DIR}feature-table.biom"
export BIOM_ANNOTATED="${BIOM_DIR}OTU_with_taxa.biom"

md "${BIOM_RAW}"

if [[ ! -s "${BIOM_RAW}" ]]
    then

    # Output: 'feature-table.biom'
    qiime tools export \
        --input-path "${FREQUENCY_TABLE}" \
        --output-format BIOMV210DirFmt \
        --output-path "${BIOM_DIR}" \
    |& tee "${LOG_DIR}tools export feature-table.biom.log"

    log "Annotate biom with taxonomy data"

    # The directory was already created
    biom add-metadata \
        --sc-separated "taxonomy" \
        --observation-metadata-fp "${TAXA_REFERENCE_HEADER}" \
        --sample-metadata-fp "${METADATA_TSV}" \
        --input-fp "${BIOM_RAW}" \
        --output-fp "${BIOM_ANNOTATED}" \
    |& tee "${LOG_DIR}biom add-metadata.log"

    log "Convert biom to JSON"

    biom convert \
        --to-json \
        --input-fp "${BIOM_ANNOTATED}" \
        --output-fp "${BIOM_DIR}OTU_with_taxa.json" \
    |& tee "${LOG_DIR}biom convert json.log"

    log "Convert biom to TSV"

    biom convert \
        --to-tsv \
        --input-fp "${BIOM_ANNOTATED}" \
        --output-fp "${BIOM_DIR}OTU_with_taxa.tsv" \
        --header-key "taxonomy" \
    |& tee "${LOG_DIR}biom convert taxa tsv.log"

    mv \
        --verbose \
        "${BIOM_ANNOTATED}" \
        "${QIIME2_FEATURES_BIOM}"

    else
        echo "Skip"
    fi



log "Completed running QIIME2 in ${QIIME2_DIR}"

chmod -R 777 "$(pwd)"

cd ..

exit 0
