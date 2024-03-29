FROM continuumio/miniconda3:4.6.14
# BASE IMAGE ACCESSING COMMAND:
# export IMG="continuumio/miniconda3:4.6.14" && docker pull "${IMG}" && docker run --rm --net=host -it "${IMG}"

ENV DEBIAN_FRONTEND noninteractive

# The software version check is hardcoded within `mgefinder`
RUN echo "Install Conda packages" && \
    conda install --quiet --yes --channel bioconda \
    bowtie2=2.3.5.0 \
    bwa=0.7.17 \
    cd-hit=4.8.1 \
    emboss=6.6.0.0 \
    htslib \
    pysam \
    python=3.6 \
    samtools=1.9 \
    snakemake=3.13.3 \
    tbb=2020.2 \
    && \
    conda clean --all --yes

RUN echo "Install Python packages" && \
    pip install mgefinder && \
    pip cache purge

ENV CONVERTER="${HOME}/scripts/convert_sequence_files.py"

RUN echo "Get the sequence conversion script" && \
    mkdir -p "$(dirname "${CONVERTER}")" && \
    curl -fsSL \
        "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/meta/scripts/convert_sequence_files.py" \
        -o "${CONVERTER}" && \
    chmod a+x "${CONVERTER}"

RUN mgefinder --help

# MANUAL BUILD COMMAND:
# export REPO="mgefinder" && export TAG="ivasilyev/${REPO}:latest" && docker build --network=host --tag "${REPO}" . && docker tag "${REPO}" "${TAG}" && docker push "${TAG}"

# PULL COMMAND:
# export IMG="ivasilyev/mgefinder:latest" && docker pull "${IMG}" && docker run --rm --net=host -it "${IMG}"
