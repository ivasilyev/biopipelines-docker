
conda install \
    q2-picrust2 \
    --yes \
    -c conda-forge \
    -c bioconda \
    -c gavinmdouglas


export PICRUST2_DIR="${TOOL_DIR}picrust2_pic_sepp/"
export PICRUST2_KO_FREQUENCY_TABLE="${PICRUST2_DIR}ko_metagenome.qza"

# Output files: 'ec_metagenome.qza', 'ko_metagenome.qza', 'pathway_abundance.qza'
qiime picrust2 full-pipeline \
    --i-seq "${REPRESENTATIVE_SEQUENCES}" \
    --i-table "${FREQUENCY_TABLE}" \
    --output-dir "${PICRUST2_DIR}" \
    --p-hsp-method pic \
    --p-placement-tool sepp \
    --p-threads $(nproc) \
    --verbose
