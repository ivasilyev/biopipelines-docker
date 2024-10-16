#!/usr/bin/env bash

export LINE="=========================================="
export _LOG_COUNTER=1

function log {
    printf "\n${LINE}\n\n[$(date '+%d-%m-%Y %H:%M:%S.%N')][QIIME2][OP#$(printf "%02d" ${_LOG_COUNTER})] ${*}\n\n${LINE}\n\n"
    _LOG_COUNTER=$((_LOG_COUNTER + 1))
}


function log_skip {
    printf "\n${LINE}\n\n[$(date '+%d-%m-%Y %H:%M:%S.%N')][Pipeline][OP#$(printf "%02d" ${_LOG_COUNTER})] Skip\n\n${LINE}\n\n"
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
export SAMPLEDATA_CSV="$(realpath "${SAMPLEDATA_CSV}")"
export METADATA_TSV="$(realpath "${METADATA_TSV}")"
export TAXA_REFERENCE_FEATURES="$(realpath "${TAXA_REFERENCE_FEATURES}")"
export TAXA_REFERENCE_CLASSIFIER="$(realpath "${TAXA_REFERENCE_CLASSIFIER}")"
export TAXA_REFERENCE_SEQUENCES="$(realpath "${TAXA_REFERENCE_SEQUENCES}")"
export TAXA_REFERENCE_HEADER="$(realpath "${TAXA_REFERENCE_HEADER}")"
# Required variables end
export TOOL_NAME="dada2"
export TOOL_DIR="${QIIME2_DIR}${TOOL_NAME}/"
export LOG_DIR="${TOOL_DIR}logs/"
export CONSENSUS_THRESHOLD=97
export GROUPING_COLUMN_NAME="SubjectID"
export CONTROL_COLUMN_NAME="Subgroup"
export CONTROL_INDICATOR_VALUE="ControlNegative"
export TAXA_COLUMN_NAME="taxonomy"
export NPROC="$(grep -c '^processor' "/proc/cpuinfo")"

mkdir -p "${LOG_DIR}"
log "Run QIIME2 in ${QIIME2_DIR}"
cd "${QIIME2_DIR}" || exit 1
qiime dev refresh-cache

log "Import and convert pre-demultiplexed paired-end FASTQ files to QIIME2 artifact"
export DEMULTIPLEXED_DIR="${QIIME2_DIR}demultiplexed_reads/"
export DEMULTIPLEXED_READS="${DEMULTIPLEXED_DIR}demultiplexed_PE_reads.qza"
md "${DEMULTIPLEXED_READS}"
if [[ ! -s "${DEMULTIPLEXED_READS}" ]]
    then
    echo "qiime tools import"
    qiime tools import \
        --input-format PairedEndFastqManifestPhred33 \
        --input-path "${SAMPLEDATA_CSV}" \
        --output-path "${DEMULTIPLEXED_READS}" \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
    |& tee "${LOG_DIR}tools import.log"
    log "Create summary statistics for raw sequences"
    echo "qiime demux summarize"
    qiime demux summarize \
        --i-data "${DEMULTIPLEXED_READS}" \
        --o-visualization "${DEMULTIPLEXED_DIR}_demultiplexed_PE_reads.qzv" \
         --verbose \
    |& tee "${LOG_DIR}demux summarize demux_PE_reads.log"
    else
        log_skip
    fi

log "Denoise demultiplexed sequences with DADA2 and calculate Amplicon Sequence Variants, ASV"

export DENOISING_DIR="${TOOL_DIR}denoising/"
export REPRESENTATIVE_SEQUENCES="${DENOISING_DIR}representative_sequences.qza"
export FREQUENCY_TABLE="${DENOISING_DIR}frequency_table.qza"
export DENOISING_STATS="${DENOISING_DIR}denoising_statistics.qza"

if [[ ! -s "${REPRESENTATIVE_SEQUENCES}" ]]
    then

    md "${REPRESENTATIVE_SEQUENCES}"

    # Q-score based filtering is built in to DADA2,
    # so doing this quality-filter step prior to denoising with DADA2 is unnecessary.
    # It also will denoise the forward and reverse reads independently and then join (merge) them
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

    else
        log_skip
    fi



log "Generate tabular view of feature identifier to sequence mapping"

qiime feature-table tabulate-seqs \
    --i-data "${REPRESENTATIVE_SEQUENCES}" \
    --o-visualization "${DENOISING_DIR}representative_sequences.qzv" \
    --verbose \
|& tee "${LOG_DIR}tabulate-seqs.log"



log "Generate a tabular view of DADA2 denoising metadata"

qiime metadata tabulate \
    --m-input-file "${DENOISING_STATS}" \
    --o-visualization "${DENOISING_DIR}denoising_statistics.qzv" \
    --verbose \
|& tee "${LOG_DIR}metadata tabulate dada2_denoising_statistics.log"



log "Try to decontaminate frequency table"
export DECONTAMINATION_DIR="${TOOL_DIR}decontam/"
export DECONTAMINATION_SCORES="${DECONTAMINATION_DIR}decontamination_scores_by_prevalence.qza"
if [[ ! -s "${DECONTAMINATION_SCORES}" ]]
    then
    md "${DECONTAMINATION_SCORES}"
    echo "qiime quality-control decontam-identify"
    qiime quality-control decontam-identify \
        --i-table "${FREQUENCY_TABLE}" \
        --m-metadata-file "${METADATA_TSV}" \
        --o-decontam-scores "${DECONTAMINATION_SCORES}" \
        --p-method prevalence \
        --p-prev-control-column "${CONTROL_COLUMN_NAME}" \
        --p-prev-control-indicator "${CONTROL_INDICATOR_VALUE}" \
        --verbose \
    |& tee "${LOG_DIR}quality-control decontam-identify.log"
    # Output: 'stats.tsv'
    echo "qiime tools export"
    qiime tools export \
        --input-path "${DECONTAMINATION_SCORES}" \
        --output-format DecontamScoreDirFmt \
        --output-path "${DECONTAMINATION_DIR}" \
    |& tee "${LOG_DIR}tools export decontam.log"
    mv \
        --verbose \
        "${DECONTAMINATION_DIR}stats.tsv" \
        "${DECONTAMINATION_DIR}decontamination_stats.tsv"
    echo "qiime quality-control decontam-score-viz"
    qiime quality-control decontam-score-viz \
        --i-decontam-scores "${DECONTAMINATION_SCORES}" \
        --i-table "${FREQUENCY_TABLE}" \
        --o-visualization "${DECONTAMINATION_DIR}decontamination_scores_by_prevalence.qzv" \
        --verbose \
    |& tee "${LOG_DIR}quality-control decontam-score-viz.log"
    else
        log_skip
    fi

export DECONTAMINATION_TABLE="${DECONTAMINATION_DIR}decontamination_filtered_frequencies.qza"
if [[ ! -s "${DECONTAMINATION_TABLE}" ]]
    then
    echo "qiime quality-control decontam-remove"
    qiime quality-control decontam-remove \
        --i-decontam-scores "${DECONTAMINATION_SCORES}" \
        --i-table "${FREQUENCY_TABLE}" \
        --o-filtered-table "${DECONTAMINATION_TABLE}" \
        --verbose \
    |& tee "${LOG_DIR}quality-control decontam-remove.log"
    else
        log_skip
    fi

if [[ -s "${DECONTAMINATION_TABLE}" ]]
    then
        echo "The decontamination was successful, use the output file: '${DECONTAMINATION_TABLE}'"
        export FREQUENCY_TABLE="${DECONTAMINATION_TABLE}"
    else
        echo "The decontamination was unsuccessful, keep use the input file: '${FREQUENCY_TABLE}'"
    fi

log "Assign taxonomy into ASV"
export TAXONOMY_DIR="${TOOL_DIR}taxonomy/"
export CLASSIFIED_TAXONOMY="${TAXONOMY_DIR}classified_taxonomy.qza"
export OTU_ASV_MAPPER="${TAXONOMY_DIR}ASV_confidences.tsv"
if [[ ! -s "${CLASSIFIED_TAXONOMY}" ]]
    then
    md "${CLASSIFIED_TAXONOMY}"
    # ASV creation is ruining OTU references.
    # Thus, a separated mapper must be created each time.
    # --p-n-jobs, The maximum number of concurrently worker processes.
    # If -1 all CPUs are used.
    # If 1 is given, no parallel computing code is used at all, which is useful for debugging.
    # For n-jobs below -1, (n_cpus + 1 + n-jobs) are used.
    # Thus for n-jobs = -2, all CPUs but one are used.
    echo "qiime feature-classifier classify-sklearn"
    qiime feature-classifier classify-sklearn \
        --p-n-jobs "-1" \
        --p-reads-per-batch 10000 \
        --i-classifier "${TAXA_REFERENCE_CLASSIFIER}" \
        --i-reads "${REPRESENTATIVE_SEQUENCES}" \
        --o-classification "${CLASSIFIED_TAXONOMY}" \
        --verbose \
    |& tee "${LOG_DIR}feature-classifier classify-sklearn.log"
    log "Create ASV-OTU mapper"
    # Output file: 'taxonomy.tsv'
    echo "qiime tools export"
    qiime tools export \
        --input-path "${CLASSIFIED_TAXONOMY}" \
        --output-format TSVTaxonomyDirectoryFormat \
        --output-path "${TAXONOMY_DIR}"
    log "Fix ASV-OTU mapper header"
    sed \
        's|Feature ID\tTaxon\tConfidence|#OTU ID\ttaxonomy\tconfidence|' \
        "${TAXONOMY_DIR}taxonomy.tsv" \
    > "${OTU_ASV_MAPPER}"
    echo "qiime metadata tabulate"
    qiime metadata tabulate \
        --m-input-file "${CLASSIFIED_TAXONOMY}" \
        --o-visualization "${TAXONOMY_DIR}classified_taxonomy.qzv" \
        --verbose \
    |& tee "${LOG_DIR}metadata tabulate classified_taxonomy.log"
    log "Plot metagenome profile bar charts"
    echo "qiime taxa barplot"
    qiime taxa barplot \
        --m-metadata-file "${METADATA_TSV}" \
        --i-table "${FREQUENCY_TABLE}" \
        --i-taxonomy "${CLASSIFIED_TAXONOMY}" \
        --o-visualization "${TAXONOMY_DIR}taxonomy_barplots.qzv" \
        --verbose \
    |& tee "${LOG_DIR}taxa barplot.log"
    else
        log_skip
    fi

log "Visualize denormalized ASV"
export SPECIES_FREQUENCY_TABLE="${TAXONOMY_DIR}species_frequencies.qza"
if [[ ! -s "${SPECIES_FREQUENCY_TABLE}" ]]
    then
    echo "qiime taxa collapse"
    qiime taxa collapse \
        --i-table "${FREQUENCY_TABLE}" \
        --i-taxonomy "${CLASSIFIED_TAXONOMY}" \
        --p-level 7 \
        --o-collapsed-table "${SPECIES_FREQUENCY_TABLE}" \
        --verbose
    echo "qiime metadata tabulate"
    qiime metadata tabulate \
        --m-input-file "${SPECIES_FREQUENCY_TABLE}" \
        --o-visualization "${TAXONOMY_DIR}species_frequencies.qzv" \
        --verbose
    else
        log_skip
    fi

log "Export denormalized ASV"
export BIOM_DIR="${TOOL_DIR}bioms/"
export BIOM_RAW="${BIOM_DIR}feature-table.biom"
export TSV_RAW="${BIOM_DIR}ASV.tsv"
export BIOM_ANNOTATED="${BIOM_DIR}ASV_with_taxa.biom"
md "${BIOM_RAW}"
if [[ ! -s "${BIOM_RAW}" ]]
    then
    # Output: 'feature-table.biom'
    echo "qiime tools export"
    qiime tools export \
        --input-path "${FREQUENCY_TABLE}" \
        --output-format BIOMV210DirFmt \
        --output-path "${BIOM_DIR}" \
    |& tee "${LOG_DIR}tools export feature-table.biom.log"
    log "Convert denormalized biom file to TSV"
    # The annotated biom file cannot include the `confidence` column from mapper
    # into output TSV file where it is required due to ridiculous format constraints,
    # so it is better to merge the raw TSV data with the mapper elsewhere
    echo "biom convert"
    biom convert \
        --header-key "taxonomy" \
        --input-fp "${BIOM_RAW}" \
        --output-fp "${TSV_RAW}" \
        --to-tsv \
    |& tee "${LOG_DIR}biom convert taxa tsv.log"
    log "Fix denormalized ASV file"
    sed \
        '/^\# .*/d' \
        --in-place \
        "${TSV_RAW}"
    log "Annotate denormalized biom with taxonomy data"
    # The directory was already created
    echo "biom add-metadata"
    biom add-metadata \
        --input-fp "${BIOM_RAW}" \
        --float-fields "confidence" \
        --observation-metadata-fp "${OTU_ASV_MAPPER}" \
        --output-fp "${BIOM_ANNOTATED}" \
        --sample-metadata-fp "${METADATA_TSV}" \
        --sc-separated "taxonomy" \
    |& tee "${LOG_DIR}biom add-metadata.log"
    log "Convert denormalized biom to JSON"
    echo "biom convert"
    biom convert \
        --input-fp "${BIOM_ANNOTATED}" \
        --output-fp "${BIOM_DIR}ASV_with_taxa.json" \
        --to-json \
    |& tee "${LOG_DIR}biom convert json.log"
    else
        log_skip
    fi

log "Perform de novo multiple sequence alignment"
export ALIGNMENTS_RAW="${TOOL_DIR}alignments/aligned_sequences.qza"
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
        log_skip
    fi



log "Filter the unconserved and highly variable and gapped columns to avoid overestimate distances"

export ALIGNMENTS_MASKED="${TOOL_DIR}masked_alignments/masked_aligned_sequences.qza"

md "${ALIGNMENTS_MASKED}"

if [[ ! -s "${ALIGNMENTS_MASKED}" ]]
    then

    qiime alignment mask \
        --i-alignment "${ALIGNMENTS_RAW}" \
        --o-masked-alignment "${ALIGNMENTS_MASKED}" \
        --verbose \
    |& tee "${LOG_DIR}alignment mask.log"

    else
        log_skip
    fi



log "Build a phylogenetic ML tree"

export UNROOTED_TREE="${TOOL_DIR}unrooted_trees/unrooted_tree.qza"

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
        log_skip
    fi



log "Root the unrooted tree based on the midpoint rooting method"

export ROOTED_TREE="${TOOL_DIR}rooted_trees/rooted_tree.qza"

md "${ROOTED_TREE}"

if [[ ! -s "${ROOTED_TREE}" ]]
    then

    qiime phylogeny midpoint-root \
        --i-tree "${UNROOTED_TREE}" \
        --o-rooted-tree "${ROOTED_TREE}" \
        --verbose \
    |& tee "${LOG_DIR}phylogeny midpoint-root.log"

    else
        log_skip
    fi



log "Analyze the core diversity using the phylogenetic pipeline"

export CORE_METRICS_DIR="${TOOL_DIR}phylogenetic_core_metrics/"

# '--output-dir' must not exist!
rm \
    --force \
    --recursive \
    --verbose \
    "${CORE_METRICS_DIR}"

# The description for output files:
# https://docs.qiime2.org/jupyterbooks/cancer-microbiome-intervention-tutorial/030-tutorial-downstream/050-core-metrics.html#core-phylogenetic-diversity-metrics
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "${ROOTED_TREE}" \
    --i-table "${FREQUENCY_TABLE}" \
    --m-metadata-file "${METADATA_TSV}" \
    --output-dir "${CORE_METRICS_DIR}" \
    --p-n-jobs-or-threads "${NPROC}" \
    --p-sampling-depth 10 \
    --verbose \
|& tee "${LOG_DIR}diversity core-metrics-phylogenetic.log"

log "Visualize alpha diversity"
find . \
    -type f \
    -name "*_vector.qza" \
    -print0 \
| xargs \
    -0 \
    --max-procs "$(nproc)" \
    -I "{}" \
        bash -c '
            FILE="{}";
            ALPHA_METRIC_DIR_NAME="${FILE%.*}/";
            BASE_NAME="$(basename "${ALPHA_METRIC_DIR_NAME}")";
            mkdir -p "${ALPHA_METRIC_DIR_NAME}";
            echo "Visualize alpha diversity data for \`${FILE}\`";
            echo "qiime diversity alpha-group-significance";
            qiime diversity alpha-group-significance \
                --i-alpha-diversity "${FILE}" \
                --m-metadata-file "${METADATA_TSV}" \
                --o-visualization "${ALPHA_METRIC_DIR_NAME}${BASE_NAME}_significance.qzv" \
                --verbose;
            # Output: alpha-diversity.tsv
            echo "qiime diversity alpha-group-significance";
            qiime tools export \
                --input-path "${FILE}" \
                --output-format AlphaDiversityDirectoryFormat \
                --output-path "${ALPHA_METRIC_DIR_NAME}";
            mv \
                --verbose \
                "${ALPHA_METRIC_DIR_NAME}alpha-diversity.tsv" \
                "${ALPHA_METRIC_DIR_NAME}${BASE_NAME}_alpha_diversity.tsv";
        '

log "Summarize frequency statistics"
export SAMPLE_FREQUENCY_DETAILS_DIR="${TOOL_DIR}sample_frequency_details/"
export SUMMARY_STATISTICS_QZV="${SAMPLE_FREQUENCY_DETAILS_DIR}frequency_table.qzv"
md "${SUMMARY_STATISTICS_QZV}"
qiime feature-table summarize \
    --i-table "${FREQUENCY_TABLE}" \
    --o-visualization "${SUMMARY_STATISTICS_QZV}" \
    --m-sample-metadata-file "${METADATA_TSV}" \
    --verbose \
|& tee "${LOG_DIR}feature-table summarize.log"

log "Export frequencies per sample"
export SAMPLE_FREQUENCY_DETAILS_CSV="${SAMPLE_FREQUENCY_DETAILS_DIR}sample-frequency-detail.csv"
# Output: directory with the file 'sample-frequency-detail.csv'
qiime tools export \
    --input-path "${SUMMARY_STATISTICS_QZV}" \
    --output-path "${SAMPLE_FREQUENCY_DETAILS_DIR}"
export SAMPLE_FREQUENCY_VALUES="${SAMPLE_FREQUENCY_DETAILS_DIR}values.txt"

log "Sort frequencies per sample"
awk \
    -F "," \
    '{print $NF}' \
    "${SAMPLE_FREQUENCY_DETAILS_CSV}" \
| sort --general-numeric-sort \
| sed 's|\..*$||g' \
> "${SAMPLE_FREQUENCY_VALUES}"
# Frequencies per sample
export MIN_FPS="$(head -n 1 "${SAMPLE_FREQUENCY_VALUES}")"
export MAX_FPS="$(tail -n 1 "${SAMPLE_FREQUENCY_VALUES}")"

log "Generate interactive alpha rarefaction curves"
# The first 2 lines are '# Constructed from biom file' and header
export ALPHA_RAREFACTION_DIR="${TOOL_DIR}alpha_rarefaction"
export ALPHA_RAREFACTION_FILE="${CORE_METRICS_DIR}alpha_rarefaction.qzv"
export ALPHA_RAREFACTION_LOG="${LOG_DIR}diversity alpha-rarefaction.log"
if [[ ! -s "${ALPHA_RAREFACTION_FILE}" ]]
    then
    md "${ALPHA_RAREFACTION_FILE}"
    echo "Use maximal frequency per sample: ${MAX_FPS}"
    echo "qiime diversity alpha-rarefaction"
    qiime diversity alpha-rarefaction \
        --m-metadata-file "${METADATA_TSV}" \
        --i-table "${FREQUENCY_TABLE}" \
        --i-phylogeny "${ROOTED_TREE}" \
        --o-visualization "${ALPHA_RAREFACTION_FILE}" \
        --p-max-depth "${MAX_FPS}" \
        --verbose \
    |& tee "${ALPHA_RAREFACTION_LOG}"
    if [[ ! -s "${ALPHA_RAREFACTION_FILE}" ]]
        then
        log "Use log from unsuccessful alpha rarefaction as it was from duck typing"
        MAX_FPS="$(
            grep \
                --perl-regexp \
                --only-matching \
                '(?<= is greater than the maximum sample total frequency of the feature_table \().*(?=\)\.)' \
                "${ALPHA_RAREFACTION_LOG}" \
            | head -n 1
        )"
        echo "Use maximal frequency per sample: ${MAX_FPS}"
        echo "qiime diversity alpha-rarefaction"
        qiime diversity alpha-rarefaction \
            --m-metadata-file "${METADATA_TSV}" \
            --i-table "${FREQUENCY_TABLE}" \
            --i-phylogeny "${ROOTED_TREE}" \
            --o-visualization "${ALPHA_RAREFACTION_FILE}" \
            --p-max-depth "${MAX_FPS}" \
            --verbose \
        |& tee -a "${ALPHA_RAREFACTION_LOG}"
        else
            echo "Alpha rarefaction wss successful from the first attempt"
        fi
    else
        log_skip
    fi

log "Visualize beta diversity"
find . \
    -type f \
    -name "*_distance_matrix.qza" \
    -print0 \
| xargs \
    -0 \
    --max-procs "$(nproc)" \
    -I "{}" \
        bash -c '
            FILE="{}";
            BETA_METRICS_DIR_NAME="${FILE%.*}/";
            BASE_NAME="$(basename "${BETA_METRICS_DIR_NAME}")";
            mkdir -p "${BETA_METRICS_DIR_NAME}";
            echo "Visualize beta diversity data for \`${FILE}\`";
            echo "qiime diversity beta-group-significance";
            qiime diversity beta-group-significance \
                --i-distance-matrix "${FILE}" \
                --m-metadata-file "${METADATA_TSV}" \
                --m-metadata-column "${CONTROL_COLUMN_NAME}" \
                --o-visualization "${BETA_METRICS_DIR_NAME}${BASE_NAME}_${CONTROL_COLUMN_NAME}_significance.qzv" \
                --p-pairwise \
                --verbose;
            # Output: distance-matrix.tsv
            echo "qiime tools export";
            qiime tools export \
                --input-path "${FILE}" \
                --output-format DistanceMatrixDirectoryFormat \
                --output-path "${BETA_METRICS_DIR_NAME}";
            mv \
                --verbose \
                "${BETA_METRICS_DIR_NAME}distance-matrix.tsv" \
                "${BETA_METRICS_DIR_NAME}${BASE_NAME}_beta_diversity.tsv";
        '

log "Calculate ANCOM, ANalysis of Composition Of Microbiomes"
export ANCOM_DIR="${TOOL_DIR}ancom/"
export ANCOM_PSEUDOCOUNT_FREQUENCY_TABLE="${ANCOM_DIR}pseudocounted_frequencies.qza"
md "${ANCOM_PSEUDOCOUNT_FREQUENCY_TABLE}"
# qiime feature-table filter-features
echo "qiime composition add-pseudocount"
qiime composition add-pseudocount \
    --i-table "${SPECIES_FREQUENCY_TABLE}" \
    --o-composition-table "${ANCOM_PSEUDOCOUNT_FREQUENCY_TABLE}" \
    --verbose
echo "qiime composition ancom"
qiime composition ancom \
    --i-table "${ANCOM_PSEUDOCOUNT_FREQUENCY_TABLE}" \
    --m-metadata-file "${METADATA_TSV}" \
    --m-metadata-column "${CONTROL_COLUMN_NAME}" \
    --o-visualization "${ANCOM_DIR}ancom_by_${CONTROL_COLUMN_NAME}.qzv" \
    --verbose

log "Completed running QIIME2 in ${QIIME2_DIR}"
chmod -R 777 "${QIIME2_DIR}"
cd ..
rm -f "$(realpath "${0}")"
exit 0
