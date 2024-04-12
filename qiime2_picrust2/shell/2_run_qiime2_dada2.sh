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

export LOG_DIR="${QIIME2_DIR}logs/"
export CONSENSUS_THRESHOLD=97
export GROUPING_COLUMN_NAME="SubjectID"
export NPROC="$(grep -c '^processor' "/proc/cpuinfo")"

mkdir -p "${LOG_DIR}"
cd "${QIIME2_DIR}" || exit 1



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



mkdir -p "${QIIME2_DIR}visualizations/"



export DENOISING_DIR="${QIIME2_DIR}dada/"
export REPRESENTATIVE_SEQUENCES="${DENOISING_DIR}REPRESENTATIVE_SEQUENCES.qza"
export FREQUENCY_TABLE="${DENOISING_DIR}dada2_frequency_table.qza"
export DENOISING_STATS="${DENOISING_DIR}dada2_denoising_statistics.qza"


if [[ ! -s "${REPRESENTATIVE_SEQUENCES}" ]]
    then

    log "DADA2 denoising"

    md "${REPRESENTATIVE_SEQUENCES}"

    qiime dada2 denoise-paired \
        --p-trunc-len-f 225 \
        --p-trunc-len-r 225 \
        --p-n-reads-learn 30000 \
        --p-n-threads "${NPROC}" \
        --i-demultiplexed-seqs "${DEMULTIPLEXED_READS}" \
        --o-representative-sequences "${REPRESENTATIVE_SEQUENCES}" \
        --o-table "${FREQUENCY_TABLE}" \
        --o-denoising-stats "${DENOISING_STATS}" \
        --verbose \
    |& tee "${LOG_DIR}dada2 denoise-paired.log"

    qiime metadata tabulate \
        --m-input-file "${DENOISING_STATS}" \
        --o-visualization "${DENOISING_DIR}dada2_denoising_statistics.qza" \
        --verbose \
    |& tee "${LOG_DIR}metadata tabulate dada2_denoising_statistics.log"

    log "Import metadata"

    export METADATA_QZV="${DENOISING_DIR}tabulated_sample_metadata.qzv"

    qiime metadata tabulate \
        --m-input-file "${METADATA_TSV}" \
        --o-visualization "${METADATA_QZV}" \
        --verbose \
    |& tee "${LOG_DIR}metadata tabulate metadata.log"

    log "Summarize statistics"

    qiime feature-table summarize \
        --i-table "${FREQUENCY_TABLE}"\
        --o-visualization "${DENOISING_DIR}dada2_frequency_table.qzv" \
        --m-sample-metadata-file "${METADATA_TSV}" \
        --verbose \
    |& tee "${LOG_DIR}feature-table summarize.log"

    qiime feature-table tabulate-seqs \
        --i-data "${REPRESENTATIVE_SEQUENCES}" \
        --o-visualization "${DENOISING_DIR}dada2_representative_sequences.qzv" \
        --verbose \
    |& tee "${LOG_DIR}tabulate-seqs.log"

    else
        echo "Skip"
    fi



export TAXONOMY_DIR="${QIIME2_DIR}taxonomy/"
export CLASSIFIED_TAXONOMY="${TAXONOMY_DIR}classified_taxonomy.qza"

if [[ ! -s "${CLASSIFIED_TAXONOMY}" ]]
    then

    log "Assign taxonomy"

    md "${CLASSIFIED_TAXONOMY}"

    # --p-n-jobs, The maximum number of concurrently worker processes. If -1 all CPUs are used. If 1 is given, no parallel computing code is used at all, which is useful for debugging. For n-jobs below -1, (n_cpus + 1 + n-jobs) are used. Thus for n-jobs = -2, all CPUs but one are used.
    qiime feature-classifier classify-sklearn \
        --p-n-jobs "-1" \
        --p-reads-per-batch 10000 \
        --i-classifier "${TAXA_REFERENCE_CLASSIFIER}" \
        --i-reads "${REPRESENTATIVE_SEQUENCES}" \
        --o-classification "${CLASSIFIED_TAXONOMY}" \
        --verbose \
    |& tee "${LOG_DIR}feature-classifier classify-sklearn.log"

    log "Create Amplicon Sequence Variant table"

    qiime metadata tabulate \
        --m-input-file "${CLASSIFIED_TAXONOMY}" \
        --o-visualization "${TAXONOMY_DIR}classified_taxonomy.qzv" \
        --verbose \
        |& tee "${LOG_DIR}metadata tabulate classified_taxonomy.log"

    log "Make prokaryotic profile"

    qiime taxa barplot \
        --m-metadata-file "${METADATA_TSV}" \
        --i-table "${FREQUENCY_TABLE}" \
        --i-taxonomy "${CLASSIFIED_TAXONOMY}" \
        --o-visualization "${TAXONOMY_DIR}taxonomy_barplots.qzv" \
        --verbose \
        |& tee "${LOG_DIR}taxa barplot.log"

    else
        echo "Skip"
    fi



log "Join paired-end reads"

export MERGED_SEQUENCES_DIR="${QIIME2_DIR}merged_reads/"
export MERGED_SEQUENCES="${QIIME2_DIR}merged_reads/merged_sequences.qza"

md "${MERGED_SEQUENCES}"

if [[ ! -s "${MERGED_SEQUENCES}" ]]
    then

    log "Join paired-end reads"

    md "${MERGED_SEQUENCES}"

    # Threads number must be within [0, 8].
    qiime vsearch merge-pairs \
        --i-demultiplexed-seqs "${DEMULTIPLEXED_READS}" \
        --o-merged-sequences "${MERGED_SEQUENCES}" \
        --p-allowmergestagger \
        --p-threads 8 \
        --verbose \
    |& tee "${LOG_DIR}vsearch merge-pairs.log"

    else
        echo "Skip"
    fi



export QUALITY_FILTERED_SEQUENCES="${QIIME2_DIR}q_score_filtered_reads/sequences_filtered_by_q_score.qza"

if [[ ! -s "${QUALITY_FILTERED_SEQUENCES}" ]]
    then

    log "Filter based on Q scores"

    md "${QUALITY_FILTERED_SEQUENCES}"

    qiime quality-filter q-score \
        --i-demux "${MERGED_SEQUENCES}" \
        --o-filtered-sequences "${QUALITY_FILTERED_SEQUENCES}" \
        --o-filter-stats "${QIIME2_DIR}q_score_filtered_reads/filtering_statistics.qza" \
        --verbose \
    |& tee "${LOG_DIR}quality-filter q-score.log"

    else
        echo "Skip"
    fi



export DEREPLICATED_DIR="${QIIME2_DIR}dereplicated/"
export DEREPLICATED_SEQUENCES="${DEREPLICATED_DIR}dereplicated_sequences.qza"
export DEREPLICATED_FREQUENCIES="${DEREPLICATED_DIR}dereplicated_frequency_table.qza"

md "${DEREPLICATED_SEQUENCES}"

if [[ ! -s "${DEREPLICATED_SEQUENCES}" ]]
    then

    log "Dereplicate sequences"

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



export CLUSTERED_DIR="${QIIME2_DIR}closed_references/"
export CLUSTERED_SEQUENCES="${CLUSTERED_DIR}closed_reference_clustered_sequences.qza"
export CLUSTERED_TABLE="${CLUSTERED_DIR}closed_reference_clustered_table.qza"


if [[ ! -s "${CLUSTERED_SEQUENCES}" ]]
    then

    log "Cluster closed references at ${CONSENSUS_THRESHOLD} percent"

    md "${CLUSTERED_SEQUENCES}"

    qiime vsearch cluster-features-closed-reference \
        --i-reference-sequences "${TAXA_REFERENCE_SEQUENCES}" \
        --i-sequences "${DEREPLICATED_SEQUENCES}" \
        --i-table "${DEREPLICATED_FREQUENCIES}" \
        --o-clustered-sequences "${CLUSTERED_SEQUENCES}" \
        --o-clustered-table "${CLUSTERED_TABLE}" \
        --o-unmatched-sequences "${QIIME2_DIR}closed_references/closed_reference_unmatched_sequences.qza" \
        --p-perc-identity 0.${CONSENSUS_THRESHOLD} \
        --p-threads "${NPROC}" \
        --verbose \
    |& tee "${LOG_DIR}vsearch cluster-features-closed-reference.log"

    else
        echo "Skip"
    fi



export DECONTAMINATED_DIR="${QIIME2_DIR}decontam/"
export DECONTAMINATED_SCORES="${DECONTAMINATED_DIR}decontam_scores_by_prevalence.qza"

if [[ ! -s "${DECONTAMINATED_SCORES}" ]]
    then

    log "Trying to decontaminate clustered sequences"

    md "${DECONTAMINATED_SCORES}"

    qiime quality-control decontam-identify \
        --i-table "${CLUSTERED_TABLE}" \
        --m-metadata-file "${METADATA_TSV}" \
        --o-decontam-scores "${DECONTAMINATED_SCORES}" \
        --p-method prevalence \
        --p-prev-control-column "Subgroup" \
        --p-prev-control-indicator "ControlNegative" \
        --verbose \
    |& tee "${LOG_DIR}quality-control decontam-identify.log"

    # Output: 'stats.tsv'
    qiime tools export \
        --input-path "${DECONTAMINATED_SCORES}" \
        --output-format DecontamScoreDirFmt \
        --output-path "${DECONTAMINATED_DIR}" \
    |& tee "${LOG_DIR}tools export decontam.log"

    qiime quality-control decontam-score-viz \
        --i-decontam-scores "${DECONTAMINATED_SCORES}" \
        --i-table "${CLUSTERED_TABLE}" \
        --o-visualization "${DECONTAMINATED_DIR}decontam_scores_by_prevalence.qzv" \
        --verbose \
    |& tee "${LOG_DIR}quality-control decontam-score-viz.log"

    else
        echo "Skip"
    fi



export DECONTAMINATED_TABLE="${DECONTAMINATED_DIR}decontam_closed_reference_clustered_table.qza"

if [[ ! -s "${DECONTAMINATED_TABLE}" ]]
    then

    qiime quality-control decontam-remove \
        --i-decontam-scores "${DECONTAMINATED_SCORES}" \
        --i-table "${CLUSTERED_TABLE}" \
        --o-filtered-table "${DECONTAMINATED_TABLE}" \
        --verbose \
    |& tee "${LOG_DIR}quality-control decontam-remove.log"

    else
        echo "Skip"
    fi

if [[ -s "${DECONTAMINATED_TABLE}" ]]
    then
        echo "The decontamination was successful, use the output file: '${DECONTAMINATED_TABLE}'"
        export CLUSTERED_TABLE="${DECONTAMINATED_TABLE}"
    else
        echo "The decontamination was unsuccessful, keep use the input file: '${CLUSTERED_TABLE}'"
    fi



log "Export the aligned sequences"

# Output: 'dna-sequences.fasta'
qiime tools export \
    --input-path "${CLUSTERED_SEQUENCES}" \
    --output-format DNASequencesDirectoryFormat \
    --output-path "${QIIME2_DIR}closed_references/" \
    |& tee "${LOG_DIR}tools export fasta.log"



log "Export an OTU table"

export BIOM_DIR="${QIIME2_DIR}bioms/"
export BIOM_RAW="${BIOM_DIR}feature-table.biom"

md "${BIOM_RAW}"

if [[ ! -s "${BIOM_RAW}" ]]
    then

    # Output: 'feature-table.biom'
    qiime tools export \
        --input-path "${CLUSTERED_TABLE}" \
        --output-path "${BIOM_DIR}" \
        --output-format BIOMV210DirFmt \
    |& tee "${LOG_DIR}tools export feature-table.biom.log"

    else
        echo "Skip"
    fi




log "Annotate biom with taxonomy data"

# The directory was already created
export BIOM_ANNOTATED="${BIOM_DIR}OTUs_with_taxa.biom"

if [[ ! -s "${BIOM_ANNOTATED}" ]]
    then

    biom add-metadata \
        --sc-separated "taxonomy" \
        --observation-metadata-fp "${TAXA_REFERENCE_HEADER}" \
        --sample-metadata-fp "${METADATA_TSV}" \
        --input-fp "${BIOM_RAW}" \
        --output-fp "${BIOM_ANNOTATED}" \
    |& tee "${LOG_DIR}biom add-metadata.log"

    else
        echo "Skip"
    fi



log "Convert biom to JSON"

biom convert \
    --to-json \
    --input-fp "${BIOM_ANNOTATED}" \
    --output-fp "${BIOM_DIR}OTUs_with_taxa.json" \
    |& tee "${LOG_DIR}biom convert json.log"



log "Convert biom to TSV"

biom convert \
    --to-tsv \
    --input-fp "${BIOM_ANNOTATED}" \
    --output-fp "${BIOM_DIR}OTUs_with_taxa.tsv" \
    --header-key "taxonomy" \
    |& tee "${LOG_DIR}biom convert taxa tsv.log"



log "Perform de novo multiple sequence alignment"

export ALIGNMENTS_RAW="${QIIME2_DIR}alignments/aligned_sequences.qza"

md "${ALIGNMENTS_RAW}"

if [[ ! -s "${ALIGNMENTS_RAW}" ]]
    then

    qiime alignment mafft \
        --p-n-threads "${NPROC}" \
        --i-sequences "${REPRESENTATIVE_SEQUENCES}" \
        --o-alignment "${ALIGNMENTS_RAW}" \
        --verbose \
    |& tee "${LOG_DIR}alignment mafft.log"

    else
        echo "Skip"
    fi



log "Filter the unconserved and highly variable and gapped columns to avoid overestimate distances"

export ALIGNMENTS_MASKED="${QIIME2_DIR}masked_alignments/masked_aligned_sequences.qza"

md "${ALIGNMENTS_MASKED}"

mkdir -p "${QIIME2_DIR}masked_alignments/"

if [[ ! -s "${ALIGNMENTS_MASKED}" ]]
    then

    qiime alignment mask \
        --i-alignment "${ALIGNMENTS_RAW}" \
        --o-masked-alignment "${ALIGNMENTS_MASKED}" \
        --verbose \
    |& tee "${LOG_DIR}alignment mask.log"

    else
        echo "Skip"
    fi



log "Build a phylogenetic ML tree"

export UNROOTED_TREE="${QIIME2_DIR}unrooted_trees/unrooted_tree.qza"

md "${UNROOTED_TREE}"

if [[ ! -s "${UNROOTED_TREE}" ]]
    then

    qiime phylogeny fasttree \
        --p-n-threads "${NPROC}" \
        --i-alignment "${ALIGNMENTS_MASKED}" \
        --o-tree "${UNROOTED_TREE}" \
        --verbose \
    |& tee "${LOG_DIR}phylogeny fasttree.log"

    else
        echo "Skip"
    fi



log "Root the unrooted tree based on the midpoint rooting method"

export ROOTED_TREE="${QIIME2_DIR}rooted_trees/rooted_tree.qza"

md "${ROOTED_TREE}"

if [[ ! -s "${ROOTED_TREE}" ]]
    then

    qiime phylogeny midpoint-root \
        --i-tree "${UNROOTED_TREE}" \
        --o-rooted-tree "${ROOTED_TREE}" \
        --verbose \
    |& tee "${LOG_DIR}phylogeny midpoint-root.log"

    else
        echo "Skip"
    fi




log "Export frequency BIOM"

export FEATURE_BIOM="${QIIME2_DIR}dada2/feature-table.biom"

md "${FEATURE_BIOM}"

if [[ ! -s "${FEATURE_BIOM}" ]]
    then

    # Output: 'feature-table.biom'
    qiime tools export \
        --input-path "${FREQUENCY_TABLE}" \
        --output-format BIOMV210DirFmt \
        --output-path "${QIIME2_DIR}dada2/"
    else
        echo "Skip"
    fi



log "Convert frequency BIOM to table"

export FEATURE_TABLE="${QIIME2_DIR}dada2/feature-table.tsv"

if [[ ! -s "${FEATURE_TABLE}" ]]
    then

    biom convert \
        --to-tsv \
        --input-fp "${FEATURE_BIOM}" \
        --output-fp "${FEATURE_TABLE}" \
    |& tee "${LOG_DIR}biom convert dada2 tsv.log"

    else
        echo "Skip"
    fi



log "Analyze the core diversity using the phylogenetic pipeline"

export CORE_METRICS_DIR="${QIIME2_DIR}phylogenetic_core_metrics/"

# '--output-dir' must not exist!
rm -rf "${CORE_METRICS_DIR}"

qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "${ROOTED_TREE}" \
    --i-table "${FREQUENCY_TABLE}" \
    --m-metadata-file "${METADATA_TSV}" \
    --output-dir "${CORE_METRICS_DIR}" \
    --p-n-jobs-or-threads "${NPROC}" \
    --p-sampling-depth 10 \
    --verbose \
|& tee "${LOG_DIR}diversity core-metrics-phylogenetic.log"

export FAITH_VECTOR="${CORE_METRICS_DIR}faith_pd_vector.qza"

qiime metadata tabulate \
    --m-input-file "${FAITH_VECTOR}" \
    --o-visualization "${CORE_METRICS_DIR}faith-pd-group-significance.qzv" \
    --verbose

# Output: alpha-diversity.tsv
qiime tools export \
    --input-path "${FAITH_VECTOR}" \
    --output-format AlphaDiversityDirectoryFormat \
    --output-path "${CORE_METRICS_DIR}" \
|& tee "${LOG_DIR}tools export faith_pd_vector.log"



log "Visualize alpha diversity"

qiime diversity alpha-group-significance \
    --i-alpha-diversity "${FAITH_VECTOR}" \
    --m-metadata-file "${METADATA_TSV}" \
    --o-visualization "${CORE_METRICS_DIR}alpha_faith_pd_group_significance.qzv" \
    --verbose \
    |& tee "${LOG_DIR}diversity alpha-group-significance faith_pd_vector.log"

qiime diversity alpha-group-significance \
    --i-alpha-diversity "${CORE_METRICS_DIR}evenness_vector.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --o-visualization "${CORE_METRICS_DIR}alpha_evenness_group_significance.qzv" \
    --verbose \
    |& tee "${LOG_DIR}diversity alpha-group-significance evenness_vector.log"

# The first 2 lines are '# Constructed from biom file' and header
export DENOISED_SAMPLES=$(( $(wc -l "${FEATURE_TABLE}" | awk '{ print $1 }') / 2 ))
export ALPHA_RAREFACTION="${CORE_METRICS_DIR}alpha_rarefaction.qzv"

if [[ ! -s "${DENOISED_SAMPLES}" ]]
    then

    qiime diversity alpha-rarefaction \
        --m-metadata-file "${METADATA_TSV}" \
        --i-table "${FREQUENCY_TABLE}" \
        --i-phylogeny "${ROOTED_TREE}" \
        --o-visualization "${ALPHA_RAREFACTION}" \
        --p-max-depth ${DENOISED_SAMPLES} \
        --verbose \
    |& tee "${LOG_DIR}diversity alpha-rarefaction.log"

    else
        echo "Skip"
    fi



log "Visualize beta diversity"

export UNIFRAC_MATRIX="${CORE_METRICS_DIR}unweighted_unifrac_distance_matrix.qza"

if [[ ! -s "${UNIFRAC_MATRIX}" ]]
    then

    qiime diversity beta-group-significance \
        --i-distance-matrix "${UNIFRAC_MATRIX}" \
        --m-metadata-file "${METADATA_TSV}" \
        --m-metadata-column "SampleSource" \
        --o-visualization "${CORE_METRICS_DIR}beta_unweighted_unifrac_SampleSource_significance.qzv" \
        --p-pairwise \
        --verbose \
    |& tee "${LOG_DIR}diversity beta-group-significance SampleSource.log"

    qiime diversity beta-group-significance \
        --i-distance-matrix "${UNIFRAC_MATRIX}" \
        --m-metadata-file "${METADATA_TSV}" \
        --m-metadata-column "${GROUPING_COLUMN_NAME}" \
        --o-visualization "${CORE_METRICS_DIR}beta_unweighted_unifrac_${GROUPING_COLUMN_NAME}_significance.qzv" \
        --p-pairwise \
        --verbose \
    |& tee "${LOG_DIR}diversity beta-group-significance ${GROUPING_COLUMN_NAME}.log"

    qiime emperor plot \
        --i-pcoa "${CORE_METRICS_DIR}unweighted_unifrac_pcoa_results.qza" \
        --m-metadata-file "${METADATA_TSV}" \
        --o-visualization "${CORE_METRICS_DIR}unweighted-unifrac-emperor.qzv" \
        --verbose \
    |& tee "${LOG_DIR}emperor plot.log"

    else
        echo "Skip"
    fi

log "Completed running QIIME2 in ${QIIME2_DIR}"

chmod -R 777 "$(pwd)"

cd ..

exit 0
