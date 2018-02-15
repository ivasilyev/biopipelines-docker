################## BASE IMAGE ######################
FROM redis

################## METADATA ######################
LABEL base.image="ubuntu:17.10"
LABEL version="1.2"
LABEL software="bwt_filtering_pipeline_docker"
LABEL description="Docker image containing all required packages to run filtering Bowtie-based pipeline"
LABEL website="https://github.com/ivasilyev/bowtie-tools"
LABEL documentation="https://github.com/ivasilyev/ThreeBees"
LABEL license="GPL-2.0"
LABEL tags="Genomics"

################## MAINTAINER ######################
MAINTAINER Ilya Vasilyev <u0412u0418u042e@gmail.com>

################## INSTALLATION ######################
# Base image CLI: docker run --rm -it redis bash
# Enable apt-get noninteractive install
ENV DEBIAN_FRONTEND noninteractive

# Software update
RUN apt-get clean all && \
    apt-get -y --force-yes update && \
    apt-get -y --force-yes upgrade

# Software installation
RUN apt-get -y --force-yes install \
    git \
    fuse \
    python3 \
    python3-pip \
    dpkg \
    curl

# Cleanup
RUN apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Tini - A tiny but valid init for containers
RUN TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

RUN pip3 install redis

# Create user docker with password docker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse docker && \
    echo `echo "root\nroot\n" | passwd root` && \
    echo `echo "docker\ndocker\n" | passwd docker` && \
    mkdir -p /data && \
    chown docker:docker /data

# Change user (CLI: su - docker)
USER docker

# Update Path variables
ENV PATH=$PATH:/home/docker/bin
ENV HOME=/home/docker

# Get the pipeline scripts
RUN mkdir ~/scripts ~/bin && \
    cd ~/scripts && \
    git clone --recurse-submodules https://github.com/ivasilyev/bwt_filtering_pipeline_docker.git && \
    cd ~/scripts/bwt_filtering_pipeline_docker && \
    find ${PWD}/ -type d -exec chmod 777 {} \; && \
    find ${PWD}/ -type f -exec chmod 644 {} \; && \
    ln -s ${PWD}/master/sampledata2jsons_queue.py /data/

# Directories to mount
VOLUME ["/data"]

# Overwrite this with 'CMD []' in a dependent Dockerfile
# CMD ["/bin/bash"]
CMD ["python3", "sampledata2jsons_queue.py"]

# Setup the default directory
WORKDIR /data

# Go to Dockerfile dir & run:
# docker build -t bwt_filtering_pipeline_master .
# docker tag bwt_filtering_pipeline_master ivasilyev/bwt_filtering_pipeline_master:latest && docker push ivasilyev/bwt_filtering_pipeline_master:latest

