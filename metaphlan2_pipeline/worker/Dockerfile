################## BASE IMAGE ######################
FROM ivasilyev/bwt_filtering_pipeline_worker

################## METADATA ######################
LABEL base.image="ivasilyev/bwt_filtering_pipeline_worker"
LABEL version="1.3"
LABEL software="metaphlan2"
LABEL description="Docker image containing all required packages to run filtering Bowtie-based pipeline"
LABEL website="https://github.com/ivasilyev/bowtie-tools"
LABEL documentation="https://github.com/ivasilyev/ThreeBees"
LABEL license="GPL-2.0"
LABEL tags="Genomics"

################## MAINTAINER ######################
MAINTAINER Ilya Vasilyev <u0412u0418u042e@gmail.com>

################## INSTALLATION ######################
# Base image CLI: docker run --rm -v /data:/data -it ivasilyev/bwt_filtering_pipeline_worker bash

# Clone the MetaPhlAn 2 repository and add it to PATH:
RUN cd ~/scripts && \
    hg clone https://bitbucket.org/biobakery/metaphlan2 && \
    cd metaphlan2 && \
    echo "export mpa_dir=${PWD}" | tee -a ~/.bashrc && \
    echo "export PATH=$PATH:${PWD}" | tee -a ~/.bashrc
