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

    log "Create summary statistics for raw sequences"

    qiime demux summarize \
        --i-data "${DEMULTIPLEXED_READS}" \
        --o-visualization "${DEMULTIPLEXED_DIR}_demultiplexed_PE_reads.qzv" \
        --verbose \
    |& tee "${LOG_DIR}demux summarize demux_PE_reads.log"

    else
        echo "Skip"
    fi



log "Merge demultiplexed paired-end reads"

export MERGED_READS="${TOOL_DIR}merged_reads/merged_PE_reads.qza"

md "${MERGED_READS}"

if [[ ! -s "${MERGED_READS}" ]]
    then

    # --p-threads INTEGER Range(0, 8, inclusive_end=True)
    # The number of threads to use for computation. Does not scale much past 4 threads.
    qiime vsearch merge-pairs \
        --i-demultiplexed-seqs "${DEMULTIPLEXED_READS}" \
        --o-merged-sequences "${MERGED_READS}" \
        --p-allowmergestagger \
        --p-threads 8 \
        --verbose

    else
        echo "Skip"
    fi



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


log "Cluster closed references at ${CONSENSUS_THRESHOLD} percent and calculate Operational Taxonomic Units, OTU"

export CLUSTERED_DIR="${TOOL_DIR}cluster_features/"
export CLUSTERED_SEQUENCES="${CLUSTERED_DIR}closed_reference_clustered_sequences.qza"
export CLUSTERED_FREQUENCIES="${CLUSTERED_DIR}closed_reference_clustered_table.qza"

if [[ ! -s "${CLUSTERED_SEQUENCES}" && ! -s "${CLUSTERED_FREQUENCIES}" ]]
    then

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



log "Run de novo chimera checking"

export DECHIMERIZATION_DIR="${TOOL_DIR}uchime-denovo/"
export CHIMERIC_SEQUENCES="${DECHIMERIZATION_DIR}chimeras.qza"
export CHIMERA_FILTERING_SEQUENCES="${DECHIMERIZATION_DIR}nonchimeras.qza"
export DECHIMERIZATION_STATS="${DECHIMERIZATION_DIR}dechimerization_statistics.qza"

md "${CHIMERIC_SEQUENCES}"

if [[ ! -s "${CHIMERIC_SEQUENCES}" ]]
    then

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
        --o-visualization "${DECHIMERIZATION_DIR}dechimerization_statistics.qzv"

    else
        echo "Skip"
    fi



log "Exclude chimeras and borderline chimeras from feature table"

export NON_CHIMERIC_FREQUENCIES="${DECHIMERIZATION_DIR}table_nonchimeric.qza"

if [[ ! -s "${NON_CHIMERIC_FREQUENCIES}" ]]
    then

    qiime feature-table filter-features \
        --i-table "${CLUSTERED_FREQUENCIES}" \
        --m-metadata-file "${CHIMERA_FILTERING_SEQUENCES}" \
        --o-filtered-table "${NON_CHIMERIC_FREQUENCIES}" \
        --p-no-exclude-ids \
        --verbose

    qiime feature-table summarize \
        --i-table "${NON_CHIMERIC_FREQUENCIES}" \
        --o-visualization "${DECHIMERIZATION_DIR}table_nonchimeric.qzv"

    else
        echo "Skip"
    fi



log "Exclude chimeras and borderline chimeras from feature sequences"

export NON_CHIMERIC_SEQUENCES="${DECHIMERIZATION_DIR}representative_sequences_nonchimeric.qza"

if [[ ! -s "${NON_CHIMERIC_SEQUENCES}" ]]
    then

    qiime feature-table filter-seqs \
        --i-data "${CLUSTERED_SEQUENCES}" \
        --m-metadata-file "${CHIMERA_FILTERING_SEQUENCES}" \
        --o-filtered-data "${NON_CHIMERIC_SEQUENCES}" \
        --p-no-exclude-ids \
        --verbose

    else
        echo "Skip"
    fi



log "Exclude chimeras but retain borderline chimeras from feature table"

export BORDERLINE_CHIMERIC_FREQUENCIES="${DECHIMERIZATION_DIR}table_borderline_chimeric.qza"
export BORDERLINE_CHIMERIC_SEQUENCES="${DECHIMERIZATION_DIR}representative_sequences_borderline_chimeric.qza"

if [[ ! -s "${BORDERLINE_CHIMERIC_FREQUENCIES}" ]]
    then

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

    qiime feature-table filter-seqs \
        --i-data "${CLUSTERED_SEQUENCES}" \
        --m-metadata-file "${CHIMERIC_SEQUENCES}" \
        --o-filtered-data "${BORDERLINE_CHIMERIC_SEQUENCES}" \
        --p-exclude-ids \
        --verbose

    else
        echo "Skip"
    fi



if [[ -s "${BORDERLINE_CHIMERIC_SEQUENCES}" && -s "${BORDERLINE_CHIMERIC_FREQUENCIES}" ]]
    then

        echo \
            "The dechimerization was successful, use the output" \
            "representative sequences: '${BORDERLINE_CHIMERIC_SEQUENCES}'" \
            "and frequency table: '${BORDERLINE_CHIMERIC_FREQUENCIES}'"

        export REPRESENTATIVE_SEQUENCES="${BORDERLINE_CHIMERIC_SEQUENCES}"

        export FREQUENCY_TABLE="${BORDERLINE_CHIMERIC_FREQUENCIES}"

    else

        echo \
            "The dechimerization was unsuccessful, keep use the input" \
            "representative sequences: '${CLUSTERED_SEQUENCES}'" \
            "and frequency table: '${CLUSTERED_FREQUENCIES}'"

        export REPRESENTATIVE_SEQUENCES="${CLUSTERED_SEQUENCES}"

        export FREQUENCY_TABLE="${CLUSTERED_FREQUENCIES}"
    fi



export DECONTAMINATION_DIR="${TOOL_DIR}decontam/"
export DECONTAMINATION_SCORES="${DECONTAMINATION_DIR}decontamination_scores_by_prevalence.qza"

if [[ ! -s "${DECONTAMINATION_SCORES}" ]]
    then

    log "Trying to decontaminate frequency table"

    md "${DECONTAMINATION_SCORES}"

    qiime quality-control decontam-identify \
        --i-table "${FREQUENCY_TABLE}" \
        --m-metadata-file "${METADATA_TSV}" \
        --o-decontam-scores "${DECONTAMINATION_SCORES}" \
        --p-method prevalence \
        --p-prev-control-column "${PREV_CONTROL_COLUMN}" \
        --p-prev-control-indicator "${PREV_CONTROL_INDICATOR}" \
        --verbose \
    |& tee "${LOG_DIR}quality-control decontam-identify.log"

    # Output: 'stats.tsv'
    qiime tools export \
        --input-path "${DECONTAMINATION_SCORES}" \
        --output-format DecontamScoreDirFmt \
        --output-path "${DECONTAMINATION_DIR}" \
    |& tee "${LOG_DIR}tools export decontam.log"

    mv \
        --verbose \
        "${DECONTAMINATION_DIR}stats.tsv" \
        "${DECONTAMINATION_DIR}decontamination_stats.tsv"

    qiime quality-control decontam-score-viz \
        --i-decontam-scores "${DECONTAMINATION_SCORES}" \
        --i-table "${FREQUENCY_TABLE}" \
        --o-visualization "${DECONTAMINATION_DIR}decontamination_scores_by_prevalence.qzv" \
        --verbose \
    |& tee "${LOG_DIR}quality-control decontam-score-viz.log"

    else
        echo "Skip"
    fi



export DECONTAMINATION_TABLE="${DECONTAMINATION_DIR}decontamination_filtered_frequencies.qza"

if [[ ! -s "${DECONTAMINATION_TABLE}" ]]
    then

    qiime quality-control decontam-remove \
        --i-decontam-scores "${DECONTAMINATION_SCORES}" \
        --i-table "${FREQUENCY_TABLE}" \
        --o-filtered-table "${DECONTAMINATION_TABLE}" \
        --verbose \
    |& tee "${LOG_DIR}quality-control decontam-remove.log"

    else
        echo "Skip"
    fi



if [[ -s "${DECONTAMINATION_TABLE}" ]]
    then
        echo "The decontamination was successful, use the output file: '${DECONTAMINATION_TABLE}'"
        export FREQUENCY_TABLE="${DECONTAMINATION_TABLE}"
    else
        echo "The decontamination was unsuccessful, keep use the input file: '${FREQUENCY_TABLE}'"
    fi



log "Export the denoised sequences to use in PCRUSt2"

export FASTA_DIR="${TOOL_DIR}fasta/"
export FASTA="${FASTA_DIR}dna-sequences.fasta"

md "${FASTA}"

if [[ ! -s "${FASTA}" ]]
    then

    # Output: 'dna-sequences.fasta'
    qiime tools export \
        --input-path "${REPRESENTATIVE_SEQUENCES}" \
        --output-format DNASequencesDirectoryFormat \
        --output-path "${FASTA_DIR}" \
    |& tee "${LOG_DIR}tools export fasta.log"

    ln \
        --symbolic \
        --verbose \
        "${FASTA}" \
        "${QIIME2_FEATURES_FASTA}"

    else
        echo "Skip"
    fi



log "Visualize denormalized OTU"

export TAXONOMY_DIR="${TOOL_DIR}taxonomy/"
export SPECIES_FREQUENCY_TABLE="${TAXONOMY_DIR}species_frequencies.qza"

md "${SPECIES_FREQUENCY_TABLE}"

if [[ ! -s "${SPECIES_FREQUENCY_TABLE}" ]]
    then

    qiime taxa collapse \
        --i-table "${FREQUENCY_TABLE}" \
        --i-taxonomy "${TAXA_REFERENCE_FEATURES}" \
        --p-level 7 \
        --o-collapsed-table "${SPECIES_FREQUENCY_TABLE}" \
        --verbose

    qiime metadata tabulate \
        --m-input-file "${SPECIES_FREQUENCY_TABLE}" \
        --o-visualization "${TAXONOMY_DIR}species_frequencies.qzv" \
        --verbose

    qiime metadata tabulate \
        --m-input-file "${TAXA_REFERENCE_CLASSIFIER}" \
        --o-visualization "${TAXONOMY_DIR}classified_taxonomy.qzv" \
        --verbose \
    |& tee "${LOG_DIR}metadata tabulate classified_taxonomy.log"

    log "Plot metagenome profile bar charts"

    qiime taxa barplot \
        --m-metadata-file "${METADATA_TSV}" \
        --i-table "${FREQUENCY_TABLE}" \
        --i-taxonomy "${TAXA_REFERENCE_CLASSIFIER}" \
        --o-visualization "${TAXONOMY_DIR}taxonomy_barplots.qzv" \
        --verbose \
    |& tee "${LOG_DIR}taxa barplot.log"

    else
        echo "Skip"
    fi



log "Export denormalized OTU"

export BIOM_DIR="${TOOL_DIR}bioms/"
export BIOM_RAW="${BIOM_DIR}feature-table.biom"
export BIOM_DENORMALIZED="${BIOM_DIR}OTU_denormalized.biom"
export BIOM_DENORMALIZED_ANNOTATED="${BIOM_DIR}OTU_denormalized_with_taxa.biom"
export TSV_DENORMALIZED_ANNOTATED="${BIOM_DIR}OTU_denormalized_with_taxa.tsv"

md "${BIOM_DENORMALIZED}"

if [[ ! -s "${BIOM_DENORMALIZED}" ]]
    then

    # Output: 'feature-table.biom'
    qiime tools export \
        --input-path "${FREQUENCY_TABLE}" \
        --output-format BIOMV210DirFmt \
        --output-path "${BIOM_DIR}" \
    |& tee "${LOG_DIR}tools export feature-table.biom.log"

    mv \
        --verbose \
        "${BIOM_RAW}" \
        "${BIOM_DENORMALIZED}"

    log "Annotate denormalized biom file with taxonomy data"

    biom add-metadata \
        --sc-separated "taxonomy" \
        --observation-metadata-fp "${TAXA_REFERENCE_HEADER}" \
        --sample-metadata-fp "${METADATA_TSV}" \
        --input-fp "${BIOM_DENORMALIZED}" \
        --output-fp "${BIOM_DENORMALIZED_ANNOTATED}" \
    |& tee "${LOG_DIR}biom add-metadata.log"

    log "Convert denormalized biom to JSON"

    biom convert \
        --to-json \
        --input-fp "${BIOM_DENORMALIZED_ANNOTATED}" \
        --output-fp "${BIOM_DIR}OTU_with_taxa.json" \
    |& tee "${LOG_DIR}biom convert json.log"

    log "Convert denormalized biom file to TSV"

    biom convert \
        --to-tsv \
        --input-fp "${BIOM_DENORMALIZED_ANNOTATED}" \
        --output-fp "${TSV_DENORMALIZED_ANNOTATED}" \
        --header-key "taxonomy" \
    |& tee "${LOG_DIR}biom convert taxa tsv.log"

    log "Fix denormalized OTU file"

    sed \
        '/^\# .*/d' \
        --in-place \
        "${TSV_DENORMALIZED_ANNOTATED}"

    log "Export the denormalized frequencies to use in PCRUSt2"

    ln \
        --symbolic \
        --verbose \
        "${BIOM_DENORMALIZED_ANNOTATED}" \
        "${QIIME2_FEATURES_BIOM}"

    else
        echo "Skip"
    fi



log "Export normalized OTU"

export NORMALIZED_FREQUENCIES="${BIOM_DIR}clustered_table_normalized.qza"
export BIOM_NORMALIZED="${BIOM_DIR}OTU_normalized.biom"
export BIOM_NORMALIZED_ANNOTATED="${BIOM_DIR}OTU_normalized_with_taxa.biom"
export TSV_NORMALIZED_ANNOTATED="${BIOM_DIR}OTU_normalized_with_taxa.tsv"

if [[ ! -s "${BIOM_NORMALIZED}" ]]
    then

    log "Normalize clustered features"

    qiime feature-table relative-frequency \
        --i-table "${FREQUENCY_TABLE}" \
        --o-relative-frequency-table "${NORMALIZED_FREQUENCIES}" \
        --verbose

    # Output: 'feature-table.biom'
    qiime tools export \
        --input-path "${NORMALIZED_FREQUENCIES}" \
        --output-format BIOMV210DirFmt \
        --output-path "${BIOM_DIR}"

    mv \
        --verbose \
        "${BIOM_RAW}" \
        "${BIOM_NORMALIZED}"

    log "Annotate normalized biom file with taxonomy data"

    biom add-metadata \
        --input-fp "${BIOM_NORMALIZED}" \
        --observation-metadata-fp "${TAXA_REFERENCE_HEADER}" \
        --output-fp "${BIOM_NORMALIZED_ANNOTATED}" \
        --sample-metadata-fp "${METADATA_TSV}" \
        --sc-separated "taxonomy"

    log "Convert normalized biom file to TSV"

    biom convert \
        --header-key "taxonomy" \
        --input-fp "${BIOM_NORMALIZED_ANNOTATED}" \
        --output-fp "${TSV_NORMALIZED_ANNOTATED}" \
        --to-tsv

    log "Fix normalized OTU TSV file"

    sed \
        '/^\# .*/d' \
        --in-place \
        "${TSV_NORMALIZED_ANNOTATED}"

    log "Export normalized frequencies to use in report"

    ln \
        --symbolic \
        --verbose \
        "${TSV_ANNOTATED}" \
        "${QIIME2_OTU_TABLE}"

    else
        echo "Skip"
    fi



log "Completed running QIIME2 in ${QIIME2_DIR}"

chmod -R 777 "${QIIME2_DIR}"

cd ..

rm -f "$(realpath "${0}")"

exit 0
