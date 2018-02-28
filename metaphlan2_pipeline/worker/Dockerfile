################## BASE IMAGE ######################
FROM ivasilyev/bwt_filtering_pipeline_worker

################## METADATA ######################
LABEL base.image="ivasilyev/bwt_filtering_pipeline_worker"
LABEL version="1.3"
LABEL software="metaphlan2"
LABEL description="Docker image containing all required packages to run MetaPhlAn2"
LABEL website="https://github.com/ivasilyev/bowtie-tools"
LABEL documentation="https://github.com/ivasilyev/ThreeBees"
LABEL license="GPL-2.0"
LABEL tags="Genomics"

################## MAINTAINER ######################
MAINTAINER Ilya Vasilyev <u0412u0418u042e@gmail.com>

################## INSTALLATION ######################
# Base image CLI: docker run --rm -v /data:/data -it ivasilyev/bwt_filtering_pipeline_worker bash

# Clone the MetaPhlAn2 repository and add it to PATH:
RUN cd ~/scripts && \
    hg clone https://bitbucket.org/biobakery/metaphlan2 && \
    cd metaphlan2 && \
    echo "export mpa_dir=${PWD}" | tee -a ~/.bashrc && \
    echo "export PATH=$PATH:${PWD}" | tee -a ~/.bashrc && \
    chmod -R 755 ${PWD}

# Get the MetaPhlAn2 custom scripts
RUN cd ~/scripts && \
    git clone https://github.com/ivasilyev/biopipelines-docker.git && \
    mv biopipelines-docker/metaphlan2_pipeline/worker/pipeline_wrapper.py .

# Get the MetaPhlAn2 reference sequences database - not implemented due to large size (~ 250 Mb)
#RUN mkdir /reference/metaphlan2 && \
#    cd /reference/metaphlan2 && \
#    wget --quiet https://bitbucket.org/biobakery/metaphlan2/downloads/mpa_v20_m200.tar && \
#    tar xopf mpa_v20_m200.tar && \
#    bzip2 -d mpa_v20_m200.fna.bz2 && \
#    python3 ~/scripts/cook_the_reference.py -i mpa_v20_m200.fna -o /reference/metaphlan2/index && \
#    rm -f mpa_v20_m200.tar

# Directories to mount
VOLUME ["/data", "/data1", "/data2", "/config", "/reference"]

# Overwrite this with 'CMD []' in a dependent Dockerfile
# CMD ["/bin/bash"]
CMD ["/bin/bash"]

# Setup the default directory
WORKDIR /home/docker
