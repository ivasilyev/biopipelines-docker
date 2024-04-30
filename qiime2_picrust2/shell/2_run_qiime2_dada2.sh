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



log "Denoise demultiplexed sequences with DADA2"

export DADA2_DIR="${QIIME2_DIR}dada2/"
export DADA2_DENOISING_DIR="${DADA2_DIR}denoising/"
export DADA2_REPRESENTATIVE_SEQUENCES="${DADA2_DENOISING_DIR}representative_sequences.qza"
export DADA2_FREQUENCY_TABLE="${DADA2_DENOISING_DIR}frequency_table.qza"
export DADA2_DENOISING_STATS="${DADA2_DENOISING_DIR}denoising_statistics.qza"

if [[ ! -s "${DADA2_REPRESENTATIVE_SEQUENCES}" ]]
    then

    md "${DADA2_REPRESENTATIVE_SEQUENCES}"

    # Q-score based filtering is built in to DADA2,
    # so doing this quality-filter step prior to denoising with DADA2 is unnecessary.
    # It also will denoise the forward and reverse reads independently and then join (merge) them
    qiime dada2 denoise-paired \
        --p-trunc-len-f 225 \
        --p-trunc-len-r 225 \
        --p-n-reads-learn 30000 \
        --p-n-threads "${NPROC}" \
        --i-demultiplexed-seqs "${DEMULTIPLEXED_READS}" \
        --o-representative-sequences "${DADA2_REPRESENTATIVE_SEQUENCES}" \
        --o-table "${DADA2_FREQUENCY_TABLE}" \
        --o-denoising-stats "${DADA2_DENOISING_STATS}" \
        --verbose \
    |& tee "${LOG_DIR}dada2 denoise-paired.log"

    else
        echo "Skip"
    fi



export DADA2_DECONTAMINATED_DIR="${DADA2_DIR}decontam/"
export DADA2_DECONTAMINATED_SCORES="${DADA2_DECONTAMINATED_DIR}decontam_scores_by_prevalence.qza"

if [[ ! -s "${DADA2_DECONTAMINATED_SCORES}" ]]
    then

    log "Trying to decontaminate frequency table"

    md "${DADA2_DECONTAMINATED_SCORES}"

    qiime quality-control decontam-identify \
        --i-table "${DADA2_FREQUENCY_TABLE}" \
        --m-metadata-file "${METADATA_TSV}" \
        --o-decontam-scores "${DADA2_DECONTAMINATED_SCORES}" \
        --p-method prevalence \
        --p-prev-control-column "${PREV_CONTROL_COLUMN}" \
        --p-prev-control-indicator "${PREV_CONTROL_INDICATOR}" \
        --verbose \
    |& tee "${LOG_DIR}quality-control decontam-identify.log"

    # Output: 'stats.tsv'
    qiime tools export \
        --input-path "${DADA2_DECONTAMINATED_SCORES}" \
        --output-format DecontamScoreDirFmt \
        --output-path "${DADA2_DECONTAMINATED_DIR}" \
    |& tee "${LOG_DIR}tools export decontam.log"

    qiime quality-control decontam-score-viz \
        --i-decontam-scores "${DADA2_DECONTAMINATED_SCORES}" \
        --i-table "${DADA2_FREQUENCY_TABLE}" \
        --o-visualization "${DADA2_DECONTAMINATED_DIR}decontam_scores_by_prevalence.qzv" \
        --verbose \
    |& tee "${LOG_DIR}quality-control decontam-score-viz.log"

    else
        echo "Skip"
    fi



export DADA2_DECONTAMINATED_TABLE="${DADA2_DECONTAMINATED_DIR}decontam_filtered_frequencies.qza"

if [[ ! -s "${DADA2_DECONTAMINATED_TABLE}" ]]
    then

    qiime quality-control decontam-remove \
        --i-decontam-scores "${DADA2_DECONTAMINATED_SCORES}" \
        --i-table "${DADA2_FREQUENCY_TABLE}" \
        --o-filtered-table "${DADA2_DECONTAMINATED_TABLE}" \
        --verbose \
    |& tee "${LOG_DIR}quality-control decontam-remove.log"

    else
        echo "Skip"
    fi



if [[ -s "${DADA2_DECONTAMINATED_TABLE}" ]]
    then
        echo "The decontamination was successful, use the output file: '${DADA2_DECONTAMINATED_TABLE}'"
        export DADA2_FREQUENCY_TABLE="${DADA2_DECONTAMINATED_TABLE}"
    else
        echo "The decontamination was unsuccessful, keep use the input file: '${DADA2_FREQUENCY_TABLE}'"
    fi



log "Generate a tabular view of DADA2 denoising metadata"

qiime metadata tabulate \
    --m-input-file "${DADA2_DENOISING_STATS}" \
    --o-visualization "${DADA2_DENOISING_DIR}denoising_statistics.qza" \
    --verbose \
|& tee "${LOG_DIR}metadata tabulate dada2_denoising_statistics.log"



log "Summarize statistics"

qiime feature-table summarize \
    --i-table "${DADA2_FREQUENCY_TABLE}"\
    --o-visualization "${DADA2_DENOISING_DIR}frequency_table.qzv" \
    --m-sample-metadata-file "${METADATA_TSV}" \
    --verbose \
|& tee "${LOG_DIR}feature-table summarize.log"



log "Generate tabular view of feature identifier to sequence mapping"

qiime feature-table tabulate-seqs \
    --i-data "${DADA2_REPRESENTATIVE_SEQUENCES}" \
    --o-visualization "${DADA2_DENOISING_DIR}representative_sequences.qzv" \
    --verbose \
|& tee "${LOG_DIR}tabulate-seqs.log"



log "Assign taxonomy as Amplicon Sequence Variants, ASV"

export DADA2_TAXONOMY_DIR="${DADA2_DIR}taxonomy/"
export DADA2_CLASSIFIED_TAXONOMY="${DADA2_TAXONOMY_DIR}classified_taxonomy.qza"

if [[ ! -s "${DADA2_CLASSIFIED_TAXONOMY}" ]]
    then

    md "${DADA2_CLASSIFIED_TAXONOMY}"

    # --p-n-jobs, The maximum number of concurrently worker processes. If -1 all CPUs are used. If 1 is given, no parallel computing code is used at all, which is useful for debugging. For n-jobs below -1, (n_cpus + 1 + n-jobs) are used. Thus for n-jobs = -2, all CPUs but one are used.
    qiime feature-classifier classify-sklearn \
        --p-n-jobs "-1" \
        --p-reads-per-batch 10000 \
        --i-classifier "${TAXA_REFERENCE_CLASSIFIER}" \
        --i-reads "${DADA2_REPRESENTATIVE_SEQUENCES}" \
        --o-classification "${DADA2_CLASSIFIED_TAXONOMY}" \
        --verbose \
    |& tee "${LOG_DIR}feature-classifier classify-sklearn.log"

    log "Create ASV table"

    qiime metadata tabulate \
        --m-input-file "${DADA2_CLASSIFIED_TAXONOMY}" \
        --o-visualization "${DADA2_TAXONOMY_DIR}classified_taxonomy.qzv" \
        --verbose \
    |& tee "${LOG_DIR}metadata tabulate classified_taxonomy.log"

    log "Make prokaryotic profile"

    qiime taxa barplot \
        --m-metadata-file "${METADATA_TSV}" \
        --i-table "${DADA2_FREQUENCY_TABLE}" \
        --i-taxonomy "${DADA2_CLASSIFIED_TAXONOMY}" \
        --o-visualization "${DADA2_TAXONOMY_DIR}taxonomy_barplots.qzv" \
        --verbose \
    |& tee "${LOG_DIR}taxa barplot.log"

    else
        echo "Skip"
    fi



log "Perform de novo multiple sequence alignment"

export ALIGNMENTS_RAW="${QIIME2_DIR}alignments/aligned_sequences.qza"

md "${ALIGNMENTS_RAW}"

if [[ ! -s "${ALIGNMENTS_RAW}" ]]
    then

    qiime alignment mafft \
        --p-n-threads "${NPROC}" \
        --i-sequences "${DADA2_REPRESENTATIVE_SEQUENCES}" \
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



log "Analyze the core diversity using the phylogenetic pipeline"

export CORE_METRICS_DIR="${QIIME2_DIR}phylogenetic_core_metrics/"

# '--output-dir' must not exist!
rm -rf "${CORE_METRICS_DIR}"

qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "${ROOTED_TREE}" \
    --i-table "${DADA2_FREQUENCY_TABLE}" \
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
        --i-table "${DADA2_FREQUENCY_TABLE}" \
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
