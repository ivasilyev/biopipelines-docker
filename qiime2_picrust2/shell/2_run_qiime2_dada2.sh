#!/usr/bin/env bash

export LINE="======================================================================================"

function log {
    printf "\n${LINE}\n\n[$(date '+%d-%m-%Y %H:%M:%S.%N')][QIIME2] $@\n\n${LINE}\n\n"
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



log Import and convert pre-demultiplexed paired-end FASTQ files to QIIME2 artifact

mkdir -p "${QIIME2_DIR}demultiplexed_reads/"

qiime tools import --input-format PairedEndFastqManifestPhred33 \
    --input-path "${SAMPLEDATA_CSV}" \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --output-path "${QIIME2_DIR}demultiplexed_reads/demultiplexed_PE_reads.qza" \
    |& tee "${LOG_DIR}tools import.log"



log Summarize sequences

mkdir -p "${QIIME2_DIR}visualizations/"

qiime demux summarize \
    --i-data "${QIIME2_DIR}demultiplexed_reads/demultiplexed_PE_reads.qza" \
    --o-visualization "${QIIME2_DIR}visualizations/demultiplexed_PE_reads.qzv" \
     --verbose \
    |& tee "${LOG_DIR}demux summarize demux_PE_reads.log"



log DADA2 denoising

mkdir -p "${QIIME2_DIR}dada2/"

qiime dada2 denoise-paired \
    --p-trunc-len-f 225 \
    --p-trunc-len-r 225 \
    --p-n-reads-learn 30000 \
    --p-n-threads "${NPROC}" \
    --i-demultiplexed-seqs "${QIIME2_DIR}demultiplexed_reads/demultiplexed_PE_reads.qza" \
    --o-representative-sequences "${QIIME2_DIR}dada2/dada2_representative_sequences.qza" \
    --o-table "${QIIME2_DIR}dada2/dada2_frequency_table.qza" \
    --o-denoising-stats "${QIIME2_DIR}dada2/dada2_denoising_statistics.qza" \
    --verbose \
    |& tee "${LOG_DIR}dada2 denoise-paired.log"

qiime metadata tabulate \
    --m-input-file "${QIIME2_DIR}dada2/dada2_denoising_statistics.qza" \
    --o-visualization "${QIIME2_DIR}visualizations/dada2_denoising_statistics.qza" \
    --verbose \
    |& tee "${LOG_DIR}metadata tabulate dada2_denoising_statistics.log"



log Import metadata

export METADATA_QZV="${QIIME2_DIR}visualizations/tabulated_sample_metadata.qzv"

qiime metadata tabulate \
    --m-input-file "${METADATA_TSV}" \
    --o-visualization "${METADATA_QZV}" \
    --verbose \
    |& tee "${LOG_DIR}metadata tabulate metadata.log"



log Summarize statistics

qiime feature-table summarize \
    --i-table "${QIIME2_DIR}dada2/dada2_frequency_table.qza"\
    --o-visualization "${QIIME2_DIR}visualizations/dada2_frequency_table.qzv" \
    --m-sample-metadata-file "${METADATA_TSV}" \
    --verbose \
    |& tee "${LOG_DIR}feature-table summarize.log"

qiime feature-table tabulate-seqs \
    --i-data "${QIIME2_DIR}dada2/dada2_representative_sequences.qza" \
    --o-visualization "${QIIME2_DIR}visualizations/dada2_representative_sequences.qzv" \
    --verbose \
    |& tee "${LOG_DIR}tabulate-seqs.log"



log Assign taxonomy

mkdir -p "${QIIME2_DIR}taxonomy/"

# --p-n-jobs, The maximum number of concurrently worker processes. If -1 all CPUs are used. If 1 is given, no parallel computing code is used at all, which is useful for debugging. For n-jobs below -1, (n_cpus + 1 + n-jobs) are used. Thus for n-jobs = -2, all CPUs but one are used.
qiime feature-classifier classify-sklearn \
    --p-n-jobs "-1" \
    --p-reads-per-batch 10000 \
    --i-classifier "${TAXA_REFERENCE_CLASSIFIER}" \
    --i-reads "${QIIME2_DIR}dada2/dada2_representative_sequences.qza" \
    --o-classification "${QIIME2_DIR}taxonomy/classified_taxonomy.qza" \
    --verbose \
    |& tee "${LOG_DIR}feature-classifier classify-sklearn.log"



log Create Amplicon Sequence Variant table

qiime metadata tabulate \
    --m-input-file "${QIIME2_DIR}taxonomy/classified_taxonomy.qza" \
    --o-visualization "${QIIME2_DIR}visualizations/classified_taxonomy.qzv" \
    --verbose \
    |& tee "${LOG_DIR}metadata tabulate classified_taxonomy.log"



log Make prokaryotic profile

qiime taxa barplot \
    --m-metadata-file "${METADATA_TSV}" \
    --i-table "${QIIME2_DIR}dada2/dada2_frequency_table.qza" \
    --i-taxonomy "${QIIME2_DIR}taxonomy/classified_taxonomy.qza" \
    --o-visualization "${QIIME2_DIR}visualizations/taxonomy_barplots.qzv" \
    --verbose \
    |& tee "${LOG_DIR}taxa barplot.log"



log Join paired-end reads

mkdir -p "${QIIME2_DIR}joined_reads/"

# Threads number must be within [0, 8].
qiime vsearch join-pairs \
    --p-threads 8 \
    --i-demultiplexed-seqs "${QIIME2_DIR}demultiplexed_reads/demultiplexed_PE_reads.qza" \
    --p-allowmergestagger \
    --o-joined-sequences "${QIIME2_DIR}joined_reads/joined_PE_reads.qza" \
    --verbose \
    |& tee "${LOG_DIR}vsearch join-pairs.log"



log Filter based on Q scores

mkdir -p "${QIIME2_DIR}q_score_filtered_reads/"

qiime quality-filter q-score \
    --i-demux "${QIIME2_DIR}joined_reads/joined_PE_reads.qza" \
    --o-filtered-sequences "${QIIME2_DIR}q_score_filtered_reads/sequences_filtered_by_q_score.qza" \
    --o-filter-stats "${QIIME2_DIR}q_score_filtered_reads/filtering_statistics.qza" \
    --verbose \
    |& tee "${LOG_DIR}quality-filter q-score.log"



log Dereplicate sequences

mkdir -p "${QIIME2_DIR}dereplicated/"

qiime vsearch dereplicate-sequences \
    --i-sequences "${QIIME2_DIR}q_score_filtered_reads/sequences_filtered_by_q_score.qza" \
    --o-dereplicated-table "${QIIME2_DIR}dereplicated/dereplicated_frequency_table.qza" \
    --o-dereplicated-sequences "${QIIME2_DIR}dereplicated/dereplicated_sequences.qza" \
    --verbose \
    |& tee "${LOG_DIR}vsearch dereplicate-sequences.log"



log Cluster closed references at ${CONSENSUS_THRESHOLD}%

mkdir -p "${QIIME2_DIR}closed_references/"

qiime vsearch cluster-features-closed-reference \
    --p-threads "${NPROC}" \
    --i-reference-sequences "${TAXA_REFERENCE_SEQUENCES}" \
    --i-table "${QIIME2_DIR}dereplicated/dereplicated_frequency_table.qza" \
    --i-sequences "${QIIME2_DIR}dereplicated/dereplicated_sequences.qza" \
    --o-clustered-table "${QIIME2_DIR}closed_references/closed_reference_clustered_table.qza" \
    --o-clustered-sequences "${QIIME2_DIR}closed_references/closed_reference_clustered_sequences.qza" \
    --o-unmatched-sequences "${QIIME2_DIR}closed_references/closed_reference_unmatched_sequences.qza" \
    --p-perc-identity 0.${CONSENSUS_THRESHOLD} \
    --verbose \
    |& tee "${LOG_DIR}vsearch cluster-features-closed-reference.log"



log Export the aligned sequences

# Output: 'dna-sequences.fasta'
qiime tools export \
    --input-path "${QIIME2_DIR}closed_references/closed_reference_clustered_sequences.qza" \
    --output-format DNASequencesDirectoryFormat \
    --output-path "${QIIME2_DIR}closed_references/" \
    |& tee "${LOG_DIR}tools export fasta.log"



log Export an OTU table

mkdir -p "${QIIME2_DIR}bioms/"

# Output: 'feature-table.biom'
qiime tools export \
    --input-path "${QIIME2_DIR}closed_references/closed_reference_clustered_table.qza" \
    --output-path "${QIIME2_DIR}bioms/" \
    --output-format BIOMV210DirFmt \
    |& tee "${LOG_DIR}tools export feature-table.biom.log"



log Annotate biom with taxonomy data

biom add-metadata \
    --sc-separated "taxonomy" \
    --observation-metadata-fp "${TAXA_REFERENCE_HEADER}" \
    --sample-metadata-fp "${METADATA_TSV}" \
    --input-fp "${QIIME2_DIR}bioms/feature-table.biom" \
    --output-fp "${QIIME2_DIR}bioms/OTUs_with_taxa.biom" \
    |& tee "${LOG_DIR}biom add-metadata.log"



log Convert biom to JSON

biom convert \
    --to-json \
    --input-fp "${QIIME2_DIR}bioms/OTUs_with_taxa.biom" \
    --output-fp "${QIIME2_DIR}bioms/OTUs_with_taxa.json" \
    |& tee "${LOG_DIR}biom convert json.log"



log Convert biom to TSV

biom convert \
    --to-tsv \
    --input-fp "${QIIME2_DIR}bioms/OTUs_with_taxa.biom" \
    --output-fp "${QIIME2_DIR}bioms/OTUs_with_taxa.tsv" \
    --header-key "taxonomy" \
    |& tee "${LOG_DIR}biom convert taxa tsv.log"



log Perform de novo multiple sequence alignment

mkdir -p "${QIIME2_DIR}alignments/"

qiime alignment mafft \
    --p-n-threads "${NPROC}" \
    --i-sequences "${QIIME2_DIR}dada2/dada2_representative_sequences.qza" \
    --o-alignment "${QIIME2_DIR}alignments/aligned_sequences.qza" \
    --verbose \
    |& tee "${LOG_DIR}alignment mafft.log"



log Filter the unconserved and highly variable and gapped columns to avoid overestimate distances

mkdir -p "${QIIME2_DIR}masked_alignments/"

qiime alignment mask \
    --i-alignment "${QIIME2_DIR}alignments/aligned_sequences.qza" \
    --o-masked-alignment "${QIIME2_DIR}masked_alignments/masked_aligned_sequences.qza" \
    --verbose \
    |& tee "${LOG_DIR}alignment mask.log"



log Build a phylogenetic ML tree

mkdir -p "${QIIME2_DIR}unrooted_trees/"

qiime phylogeny fasttree \
    --p-n-threads "${NPROC}" \
    --i-alignment "${QIIME2_DIR}masked_alignments/masked_aligned_sequences.qza" \
    --o-tree "${QIIME2_DIR}unrooted_trees/unrooted_tree.qza" \
    --verbose \
    |& tee "${LOG_DIR}phylogeny fasttree.log"



log Root the unrooted tree based on the midpoint rooting method

mkdir -p "${QIIME2_DIR}rooted_trees/"

qiime phylogeny midpoint-root \
    --i-tree "${QIIME2_DIR}unrooted_trees/unrooted_tree.qza" \
    --o-rooted-tree "${QIIME2_DIR}rooted_trees/rooted_tree.qza" \
    --verbose \
    |& tee "${LOG_DIR}phylogeny midpoint-root.log"



# Output: 'feature-table.biom'
qiime tools export \
    --input-path "${QIIME2_DIR}dada2/dada2_frequency_table.qza" \
    --output-format BIOMV210DirFmt \
    --output-path "${QIIME2_DIR}dada2/"



biom convert \
    --to-tsv \
    --input-fp "${QIIME2_DIR}dada2/feature-table.biom" \
    --output-fp "${QIIME2_DIR}dada2/feature-table.tsv" \
    |& tee "${LOG_DIR}biom convert dada2 tsv.log"



log Analyze the core diversity using the phylogenetic pipeline

# '--output-dir' must not exist!
rm -rf "${QIIME2_DIR}phylogenetic_core_metrics/"

qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "${QIIME2_DIR}rooted_trees/rooted_tree.qza" \
    --i-table "${QIIME2_DIR}dada2/dada2_frequency_table.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --output-dir "${QIIME2_DIR}phylogenetic_core_metrics/" \
    --p-n-jobs-or-threads "${NPROC}" \
    --p-sampling-depth 10 \
    --verbose \
    |& tee "${LOG_DIR}diversity core-metrics-phylogenetic.log"

qiime metadata tabulate \
    --m-input-file "${QIIME2_DIR}phylogenetic_core_metrics/faith_pd_vector.qza" \
    --o-visualization "${QIIME2_DIR}visualizations/faith-pd-group-significance.qzv" \
    --verbose

# Output: alpha-diversity.tsv
qiime tools export \
    --input-path "${QIIME2_DIR}phylogenetic_core_metrics/faith_pd_vector.qza" \
    --output-format AlphaDiversityDirectoryFormat \
    --output-path "${QIIME2_DIR}phylogenetic_core_metrics/" \
    |& tee "${LOG_DIR}tools export faith_pd_vector.log"



log Visualize alpha diversity

qiime diversity alpha-group-significance \
    --i-alpha-diversity "${QIIME2_DIR}phylogenetic_core_metrics/faith_pd_vector.qza" \
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
export DENOISED_SAMPLES=$(( $(wc -l "${QIIME2_DIR}dada2/feature-table.tsv" | awk '{ print $1 }') - 2 ))

qiime diversity alpha-rarefaction \
    --m-metadata-file "${METADATA_TSV}" \
    --i-table "${QIIME2_DIR}dada2/dada2_frequency_table.qza" \
    --i-phylogeny "${QIIME2_DIR}rooted_trees/rooted_tree.qza" \
    --o-visualization "${QIIME2_DIR}visualizations/alpha_rarefaction.qzv" \
    --p-max-depth ${DENOISED_SAMPLES} \
    --verbose \
    |& tee "${LOG_DIR}diversity alpha-rarefaction.log"



log Visualize beta diversity

qiime diversity beta-group-significance \
    --i-distance-matrix "${QIIME2_DIR}phylogenetic_core_metrics/unweighted_unifrac_distance_matrix.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --m-metadata-column "SampleSource" \
    --o-visualization "${QIIME2_DIR}visualizations/beta_unweighted_unifrac_SampleSource_significance.qzv" \
    --p-pairwise \
    --verbose \
    |& tee "${LOG_DIR}diversity beta-group-significance SampleSource.log"

qiime diversity beta-group-significance \
    --i-distance-matrix "${QIIME2_DIR}phylogenetic_core_metrics/unweighted_unifrac_distance_matrix.qza" \
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



log Test differential abundance with ANCOM

mkdir -p "${QIIME2_DIR}differential_abundances/"

qiime feature-table filter-samples \
    --i-table "${QIIME2_DIR}dada2/dada2_frequency_table.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --o-filtered-table "${QIIME2_DIR}differential_abundances/source_features_table.qza" \
    --verbose \
    |& tee "${LOG_DIR}feature-table filter-samples.log"

qiime composition add-pseudocount \
    --i-table "${QIIME2_DIR}differential_abundances/source_features_table.qza" \
    --o-composition-table "${QIIME2_DIR}differential_abundances/pseudocounted_features_table.qza" \
    --verbose \
    |& tee "${LOG_DIR}composition add-pseudocount.log"

qiime composition ancom \
    --i-table "${QIIME2_DIR}differential_abundances/pseudocounted_features_table.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --m-metadata-column "${GROUPING_COLUMN_NAME}" \
    --o-visualization "${QIIME2_DIR}visualizations/ancom_${GROUPING_COLUMN_NAME}.qzv" \
    --verbose \
    |& tee "${LOG_DIR}composition ancom.log"



export COLLAPSE_LEVEL=6

log Test differential abundance with ANCOM for collapsed features

mkdir -p "${QIIME2_DIR}collapsed_differential_abundances/"

qiime taxa collapse \
    --i-table "${QIIME2_DIR}differential_abundances/source_features_table.qza" \
    --i-taxonomy "${TAXA_REFERENCE_FEATURES}" \
    --p-level ${COLLAPSE_LEVEL} \
    --o-collapsed-table "${QIIME2_DIR}collapsed_differential_abundances/collapsed_source_features_table.qza" \
    --verbose \
    |& tee "${LOG_DIR}taxa collapse.log"

qiime composition add-pseudocount \
    --i-table "${QIIME2_DIR}collapsed_differential_abundances/collapsed_source_features_table.qza" \
    --o-composition-table "${QIIME2_DIR}collapsed_differential_abundances/collapsed_pseudocounted_features_table.qza" \
    --verbose \
    |& tee "${LOG_DIR}composition add-pseudocount collapsed.log"

qiime composition ancom \
    --i-table "${QIIME2_DIR}collapsed_differential_abundances/collapsed_pseudocounted_features_table.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --m-metadata-column "${GROUPING_COLUMN_NAME}" \
    --o-visualization "${QIIME2_DIR}visualizations/collapsed_ancom_${GROUPING_COLUMN_NAME}.qzv" \
    |& tee "${LOG_DIR}composition ancom collapsed.log"



log "Completed running QIIME2 in ${QIIME2_DIR}"

chmod -R 777 "$(pwd)"

cd ..

exit 0
