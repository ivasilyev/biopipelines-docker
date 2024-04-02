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

export DEMULTIPLEXED_READS="${QIIME2_DIR}demultiplexed_reads/demultiplexed_PE_reads.qza"

md "${DEMULTIPLEXED_READS}"

qiime tools import \
    --input-format PairedEndFastqManifestPhred33 \
    --input-path "${SAMPLEDATA_CSV}" \
    --output-path "${DEMULTIPLEXED_READS}" \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    |& tee "${LOG_DIR}tools import.log"



log "Summarize sequences"

mkdir -p "${QIIME2_DIR}visualizations/"

qiime demux summarize \
    --i-data "${DEMULTIPLEXED_READS}" \
    --o-visualization "${QIIME2_DIR}visualizations/demultiplexed_PE_reads.qzv" \
     --verbose \
    |& tee "${LOG_DIR}demux summarize demux_PE_reads.log"



log "DADA2 denoising"

export REPRESENTATIVE_SEQUENCES="${QIIME2_DIR}dada2/REPRESENTATIVE_SEQUENCES.qza"
export FREQUENCY_TABLE="${QIIME2_DIR}dada2/dada2_frequency_table.qza"
export DENOISING_STATS="${QIIME2_DIR}dada2/dada2_denoising_statistics.qza"

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
    --o-visualization "${QIIME2_DIR}visualizations/dada2_denoising_statistics.qza" \
    --verbose \
    |& tee "${LOG_DIR}metadata tabulate dada2_denoising_statistics.log"



log "Import metadata"

export METADATA_QZV="${QIIME2_DIR}visualizations/tabulated_sample_metadata.qzv"

qiime metadata tabulate \
    --m-input-file "${METADATA_TSV}" \
    --o-visualization "${METADATA_QZV}" \
    --verbose \
    |& tee "${LOG_DIR}metadata tabulate metadata.log"



log "Summarize statistics"

qiime feature-table summarize \
    --i-table "${FREQUENCY_TABLE}"\
    --o-visualization "${QIIME2_DIR}visualizations/dada2_frequency_table.qzv" \
    --m-sample-metadata-file "${METADATA_TSV}" \
    --verbose \
    |& tee "${LOG_DIR}feature-table summarize.log"

qiime feature-table tabulate-seqs \
    --i-data "${REPRESENTATIVE_SEQUENCES}" \
    --o-visualization "${QIIME2_DIR}visualizations/dada2_representative_sequences.qzv" \
    --verbose \
    |& tee "${LOG_DIR}tabulate-seqs.log"



log "Assign taxonomy"

export CLASSIFIED_TAXONOMY="${QIIME2_DIR}taxonomy/classified_taxonomy.qza"

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
    --o-visualization "${QIIME2_DIR}visualizations/classified_taxonomy.qzv" \
    --verbose \
    |& tee "${LOG_DIR}metadata tabulate classified_taxonomy.log"



log "Make prokaryotic profile"

qiime taxa barplot \
    --m-metadata-file "${METADATA_TSV}" \
    --i-table "${FREQUENCY_TABLE}" \
    --i-taxonomy "${CLASSIFIED_TAXONOMY}" \
    --o-visualization "${QIIME2_DIR}visualizations/taxonomy_barplots.qzv" \
    --verbose \
    |& tee "${LOG_DIR}taxa barplot.log"



log "Join paired-end reads"

export MERGED_SEQUENCES_DIR="${QIIME2_DIR}merged_reads/"
export MERGED_SEQUENCES="${QIIME2_DIR}merged_reads/merged_sequences.qza"

md "${MERGED_SEQUENCES}"

# Threads number must be within [0, 8].
qiime vsearch merge-pairs \
    --i-demultiplexed-seqs "${DEMULTIPLEXED_READS}" \
    --o-merged-sequences "${MERGED_SEQUENCES}" \
    --p-allowmergestagger \
    --p-threads 8 \
    --verbose \
    |& tee "${LOG_DIR}vsearch merge-pairs.log"



log "Filter based on Q scores"

export QUALITY_FILTERED_SEQUENCES="${QIIME2_DIR}q_score_filtered_reads/sequences_filtered_by_q_score.qza"

md "${QUALITY_FILTERED_SEQUENCES}"

qiime quality-filter q-score \
    --i-demux "${MERGED_SEQUENCES}" \
    --o-filtered-sequences "${QUALITY_FILTERED_SEQUENCES}" \
    --o-filter-stats "${QIIME2_DIR}q_score_filtered_reads/filtering_statistics.qza" \
    --verbose \
    |& tee "${LOG_DIR}quality-filter q-score.log"



log "Dereplicate sequences"

export DEREPLICATED_DIR="${QIIME2_DIR}dereplicated/"
export DEREPLICATED_SEQUENCES="${DEREPLICATED_DIR}dereplicated_sequences.qza"
export DEREPLICATED_FREQUENCIES="${DEREPLICATED_DIR}dereplicated_frequency_table.qza"

md "${DEREPLICATED_SEQUENCES}"

qiime vsearch dereplicate-sequences \
    --i-sequences "${QUALITY_FILTERED_SEQUENCES}" \
    --o-dereplicated-table "${DEREPLICATED_FREQUENCIES}" \
    --o-dereplicated-sequences "${DEREPLICATED_FREQUENCIES}" \
    --verbose \
    |& tee "${LOG_DIR}vsearch dereplicate-sequences.log"



log "Cluster closed references at ${CONSENSUS_THRESHOLD} percent"

export CLUSTERED_DIR="${QIIME2_DIR}closed_references/"
export CLUSTERED_SEQUENCES="${CLUSTERED_DIR}closed_reference_clustered_sequences.qza"
export CLUSTERED_TABLE="${CLUSTERED_DIR}closed_reference_clustered_table.qza"

md "${CLUSTERED_SEQUENCES}"

qiime vsearch cluster-features-closed-reference \
    --p-threads "${NPROC}" \
    --i-reference-sequences "${TAXA_REFERENCE_SEQUENCES}" \
    --i-table "${DEREPLICATED_FREQUENCIES}" \
    --i-sequences "${DEREPLICATED_FREQUENCIES}" \
    --o-clustered-table "${CLUSTERED_TABLE}" \
    --o-clustered-sequences "${CLUSTERED_SEQUENCES}" \
    --o-unmatched-sequences "${QIIME2_DIR}closed_references/closed_reference_unmatched_sequences.qza" \
    --p-perc-identity 0.${CONSENSUS_THRESHOLD} \
    --verbose \
    |& tee "${LOG_DIR}vsearch cluster-features-closed-reference.log"



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

# Output: 'feature-table.biom'
qiime tools export \
    --input-path "${CLUSTERED_TABLE}" \
    --output-path "${BIOM_DIR}" \
    --output-format BIOMV210DirFmt \
    |& tee "${LOG_DIR}tools export feature-table.biom.log"



log "Annotate biom with taxonomy data"

# The directory was already created
export BIOM_ANNOTATED="${BIOM_DIR}OTUs_with_taxa.biom"

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

qiime alignment mafft \
    --p-n-threads "${NPROC}" \
    --i-sequences "${REPRESENTATIVE_SEQUENCES}" \
    --o-alignment "${ALIGNMENTS_RAW}" \
    --verbose \
    |& tee "${LOG_DIR}alignment mafft.log"



log "Filter the unconserved and highly variable and gapped columns to avoid overestimate distances"

export ALIGNMENTS_MASKED="${QIIME2_DIR}masked_alignments/masked_aligned_sequences.qza"

md "${ALIGNMENTS_MASKED}"

mkdir -p "${QIIME2_DIR}masked_alignments/"

qiime alignment mask \
    --i-alignment "${ALIGNMENTS_RAW}" \
    --o-masked-alignment "${ALIGNMENTS_MASKED}" \
    --verbose \
    |& tee "${LOG_DIR}alignment mask.log"



log "Build a phylogenetic ML tree"

export UNROOTED_TREE="${QIIME2_DIR}unrooted_trees/unrooted_tree.qza"

md "${UNROOTED_TREE}"

qiime phylogeny fasttree \
    --p-n-threads "${NPROC}" \
    --i-alignment "${ALIGNMENTS_MASKED}" \
    --o-tree "${UNROOTED_TREE}" \
    --verbose \
    |& tee "${LOG_DIR}phylogeny fasttree.log"



log "Root the unrooted tree based on the midpoint rooting method"

export ROOTED_TREE="${QIIME2_DIR}rooted_trees/rooted_tree.qza"

md "${ROOTED_TREE}"

qiime phylogeny midpoint-root \
    --i-tree "${UNROOTED_TREE}" \
    --o-rooted-tree "${ROOTED_TREE}" \
    --verbose \
    |& tee "${LOG_DIR}phylogeny midpoint-root.log"



log "Export frequency BIOM"

export FEATURE_BIOM="${QIIME2_DIR}dada2/feature-table.biom" 

md "${FEATURE_BIOM}"

# Output: 'feature-table.biom'
qiime tools export \
    --input-path "${FREQUENCY_TABLE}" \
    --output-format BIOMV210DirFmt \
    --output-path "${QIIME2_DIR}dada2/"



log "Convert frequency BIOM to table"

export FEATURE_TABLE="${QIIME2_DIR}dada2/feature-table.tsv"

biom convert \
    --to-tsv \
    --input-fp "${FEATURE_BIOM}" \
    --output-fp "${FEATURE_TABLE}" \
    |& tee "${LOG_DIR}biom convert dada2 tsv.log"



log "Analyze the core diversity using the phylogenetic pipeline"

# '--output-dir' must not exist!
rm -rf "${QIIME2_DIR}phylogenetic_core_metrics/"

qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "${ROOTED_TREE}" \
    --i-table "${FREQUENCY_TABLE}" \
    --m-metadata-file "${METADATA_TSV}" \
    --output-dir "${QIIME2_DIR}phylogenetic_core_metrics/" \
    --p-n-jobs-or-threads "${NPROC}" \
    --p-sampling-depth 10 \
    --verbose \
    |& tee "${LOG_DIR}diversity core-metrics-phylogenetic.log"

export FAITH_VECTOR="${QIIME2_DIR}phylogenetic_core_metrics/faith_pd_vector.qza"

qiime metadata tabulate \
    --m-input-file "${FAITH_VECTOR}" \
    --o-visualization "${QIIME2_DIR}visualizations/faith-pd-group-significance.qzv" \
    --verbose

# Output: alpha-diversity.tsv
qiime tools export \
    --input-path "${FAITH_VECTOR}" \
    --output-format AlphaDiversityDirectoryFormat \
    --output-path "${QIIME2_DIR}phylogenetic_core_metrics/" \
    |& tee "${LOG_DIR}tools export faith_pd_vector.log"



log "Visualize alpha diversity"

qiime diversity alpha-group-significance \
    --i-alpha-diversity "${FAITH_VECTOR}" \
    --m-metadata-file "${METADATA_TSV}" \
    --o-visualization "${QIIME2_DIR}visualizations/alpha_faith_pd_group_significance.qzv" \
    --verbose \
    |& tee "${LOG_DIR}diversity alpha-group-significance faith_pd_vector.log"

qiime diversity alpha-group-significance \
    --i-alpha-diversity "${QIIME2_DIR}phylogenetic_core_metrics/evenness_vector.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --o-visualization "${QIIME2_DIR}visualizations/alpha_evenness_group_significance.qzv" \
    --verbose \
    |& tee "${LOG_DIR}diversity alpha-group-significance evenness_vector.log"

# The first 2 lines are '# Constructed from biom file' and header
export DENOISED_SAMPLES=$(( $(wc -l "${FEATURE_TABLE}" | awk '{ print $1 }') / 2 ))

qiime diversity alpha-rarefaction \
    --m-metadata-file "${METADATA_TSV}" \
    --i-table "${FREQUENCY_TABLE}" \
    --i-phylogeny "${ROOTED_TREE}" \
    --o-visualization "${QIIME2_DIR}visualizations/alpha_rarefaction.qzv" \
    --p-max-depth ${DENOISED_SAMPLES} \
    --verbose \
    |& tee "${LOG_DIR}diversity alpha-rarefaction.log"



log "Visualize beta diversity"

export UNIFRAC_MATRIX="${QIIME2_DIR}phylogenetic_core_metrics/unweighted_unifrac_distance_matrix.qza"

qiime diversity beta-group-significance \
    --i-distance-matrix "${UNIFRAC_MATRIX}" \
    --m-metadata-file "${METADATA_TSV}" \
    --m-metadata-column "SampleSource" \
    --o-visualization "${QIIME2_DIR}visualizations/beta_unweighted_unifrac_SampleSource_significance.qzv" \
    --p-pairwise \
    --verbose \
    |& tee "${LOG_DIR}diversity beta-group-significance SampleSource.log"

qiime diversity beta-group-significance \
    --i-distance-matrix "${UNIFRAC_MATRIX}" \
    --m-metadata-file "${METADATA_TSV}" \
    --m-metadata-column "${GROUPING_COLUMN_NAME}" \
    --o-visualization "${QIIME2_DIR}visualizations/beta_unweighted_unifrac_${GROUPING_COLUMN_NAME}_significance.qzv" \
    --p-pairwise \
    --verbose \
    |& tee "${LOG_DIR}diversity beta-group-significance ${GROUPING_COLUMN_NAME}.log"

qiime emperor plot \
    --i-pcoa "${QIIME2_DIR}phylogenetic_core_metrics/unweighted_unifrac_pcoa_results.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --o-visualization "${QIIME2_DIR}visualizations/unweighted-unifrac-emperor.qzv" \
    --verbose \
    |& tee "${LOG_DIR}emperor plot.log"

log "Completed running QIIME2 in ${QIIME2_DIR}"

chmod -R 777 "$(pwd)"

cd ..

exit 0
