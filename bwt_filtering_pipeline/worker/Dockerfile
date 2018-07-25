################## BASE IMAGE ######################
FROM ubuntu:17.10

################## METADATA ######################
LABEL base.image="ubuntu:18.04"
LABEL pipeline="bwt_filtering"
LABEL name="bwt_filtering_pipeline_worker"
LABEL description="Docker image containing all required packages to run filtering Bowtie-based pipeline"
LABEL website="https://github.com/ivasilyev/bowtie-tools"
LABEL documentation="https://github.com/ivasilyev/ThreeBees"
LABEL license="GPL-2.0"
LABEL tags="Genomics"
LABEL maintainer="Ilya Vasilyev <u0412u0418u042e@gmail.com>"

################## INSTALLATION ######################
# Base image CLI: docker run --rm --net=host -it ubuntu:18.04 bash

# Enable apt-get noninteractive install
ENV DEBIAN_FRONTEND noninteractive

RUN echo "Software update" && \
    apt-get -y update && \
    apt-get -y install \
            software-properties-common \
            apt-utils && \
    add-apt-repository universe && \
    add-apt-repository multiverse

RUN echo "Install packages" && \
    apt-get clean all && \
    apt-get -y update && \
    apt-get -y upgrade && \
    apt-get -y install \
               autotools-dev   \
               automake        \
               autoconf        \
               build-essential \
               cmake           \
               flex            \
               curl            \
               grep            \
               sed             \
               dpkg            \
               fuse            \
               git             \
               wget            \
               zip             \
               nano            \
               openjdk-8-jre   \
               pkg-config      \
               python          \
               python-dev      \
               python-pip      \
               python3         \
               python3-pip     \
               bzip2           \
               ca-certificates \
               lib32ncurses5-dev \
               libatlas-base-dev \
               libblas-dev     \
               libboost-all-dev \
               libbz2-dev      \
               libcurl4-openssl-dev \
               libev-dev       \
               libglib2.0-0    \
               libhdf5-dev     \
               liblapack-dev   \
               liblzma-dev     \
               libsm6          \
               libssl-dev      \
               libtbb2         \
               libtbb-dev      \
               libxext6        \
               libxrender1     \
               mercurial       \
               subversion      \
               zlib1g-dev

RUN echo "Install Tini - A tiny but valid init for containers" && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb

RUN echo "Cleanup" && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN echo "Deploy Bowtie" && \
    cd /tmp && \
    export BOWTIE_VERSION="1.2.2" && \
    curl -fsSL https://github.com/BenLangmead/bowtie/releases/download/v1.2.2_p1/bowtie-1.2.2-linux-x86_64.zip -o bowtie.zip && \
    unzip bowtie.zip && \
    cd bowtie-${BOWTIE_VERSION}-linux-x86_64 && \
    cp bowtie bowtie-align-l bowtie-align-s bowtie-build bowtie-build-l bowtie-build-s bowtie-inspect bowtie-inspect-l bowtie-inspect-s /usr/local/bin/ && \
    cd /usr/local/bin/ && \
    chmod ugo+x bowtie bowtie-align-l bowtie-align-s bowtie-build bowtie-build-l bowtie-build-s bowtie-inspect bowtie-inspect-l bowtie-inspect-s && \
    cd /tmp && \
    rm -rf bowtie-${BOWTIE_VERSION}-linux-x86_64 bowtie.zip

RUN echo "Deploy Bowtie2" && \
    cd /tmp && \
    export BOWTIE2_VERSION="2.3.4.2-alpha.2" && \
    curl -fsSL https://github.com/BenLangmead/bowtie2/releases/download/v2.3.4.2-alpha.2/bowtie2-2.3.4.2-alpha.2-linux-x86_64.zip -o bowtie2.zip && \
    unzip bowtie2.zip && \
    cd bowtie2-${BOWTIE2_VERSION}-linux-x86_64 && \
    cp bowtie2 bowtie2-align-l bowtie2-align-s bowtie2-build bowtie2-build-l bowtie2-build-s bowtie2-inspect bowtie2-inspect-l bowtie2-inspect-s /usr/local/bin/ && \
    cd /usr/local/bin/ && \
    chmod ugo+x bowtie2 bowtie2-align-l bowtie2-align-s bowtie2-build bowtie2-build-l bowtie2-build-s bowtie2-inspect bowtie2-inspect-l bowtie2-inspect-s && \
    cd /tmp && \
    rm -rf bowtie2-${BOWTIE2_VERSION}-linux-x86_64 bowtie2.zip

RUN echo "Install SAMtools" && \
    mkdir -p /tmp/samtools && \
    cd /tmp/samtools && \
    export SAMTOOLS_VERSION="1.9" && \
    printf "\nCompile HTSlib\n\n" && \
    curl -fsSL https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 -o htslib.tar.bz2 && \
    tar xvjf htslib.tar.bz2 && \
    cd htslib-${SAMTOOLS_VERSION} && \
    autoheader && \
    autoconf && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    printf "\nCompile SAMtools\n\n" && \
    curl -fsSL https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 -o samtools.tar.bz2 && \
    tar xvjf samtools.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure && \
    make && \
    make install && \
    cd /tmp && \
    rm -rf samtools

RUN echo "Compile BEDtools2"
    mkdir -p /tmp/bedtools2 && \
    cd /tmp/bedtools2 && \
    export BEDTOOLS_VERSION="2.27.1" && \
    curl -fsSL https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz -o bedtools2.tar.gz && \
    tar xf bedtools2.tar.gz && \
    cd bedtools2 && \
    make && \
    make install && \
    cd /tmp && \
    rm -rf bedtools2

# Install Python packages
RUN pip3 install numpy pandas redis

# Create user docker with password docker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse docker && \
    echo `echo "docker\ndocker\n" | passwd docker` && \
    mkdir /data /config ${REF_DIR} && \
    chown docker:docker /data /config ${REF_DIR} && \
    chmod -R 755 /data /config ${REF_DIR}

# Change user (CLI: su - docker)
USER docker

# Update Path variables
ENV PATH=$PATH:/home/docker/bin
ENV HOME=/home/docker

# Get the pipeline scripts
RUN mkdir ~/scripts && \
    cd ~/scripts && \
    git clone --recurse-submodules https://github.com/ivasilyev/bowtie-tools.git && \
    git clone https://github.com/ivasilyev/biopipelines-docker.git && \
    mv bowtie-tools/nBee.py bowtie-tools/cook_the_reference.py bowtie-tools/sam2coverage.py bowtie-tools/verify_coverages.py biopipelines-docker/bwt_filtering_pipeline/worker/worker.py . && \
    rm -rf bowtie-tools biopipelines-docker && \
    chmod -R +x ~/scripts && \
    echo "export PATH=$PATH:${PWD}" | tee -a ~/.bashrc

# Directories to mount
VOLUME ["/data", "/config", "/reference"]

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["python3", "/home/docker/scripts/pipeline_wrapper.py"]

# Setup the default directory
WORKDIR /home/docker

# MANUAL BUILD COMMAND:
# docker build -t bwt_filtering_pipeline_worker . && docker tag bwt_filtering_pipeline_worker ivasilyev/bwt_filtering_pipeline_worker:latest && docker push ivasilyev/bwt_filtering_pipeline_worker:latest
