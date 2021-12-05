FROM continuumio/miniconda3:4.6.14
# BASE IMAGE ACCESSING COMMAND:
# export IMG="continuumio/miniconda3:4.6.14" && docker pull "${IMG}" && docker run --rm --net=host -it "${IMG}"

# The software version check is hardcoded within `mgefinder`
RUN conda install --quiet --yes --channel bioconda \
    bowtie2=2.3.5.0 \
    cd-hit=4.8.1 \
    emboss=6.6.0.0 \
    htslib \
    pysam \
    python=3.6 \
    samtools=1.9 \
    snakemake=3.13.3 \
    tbb=2020.2

RUN pip install --quiet mgefinder

RUN mgefinder --help

# MANUAL BUILD COMMAND:
# export IMG="mgefinder" && export TAG="ivasilyev/${IMG}:latest" && docker build --network=host --tag "${IMG}" . && docker tag "${IMG}" "${TAG}" && docker push "${TAG}"
