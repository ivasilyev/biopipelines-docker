FROM ubuntu:18.04
# Base image CLI: docker run --rm --net=host -it ubuntu:18.04 bash
# In the case of network failure, try:
# echo deb http://mirror.truenetwork.ru/ubuntu/ bionic-backports main restricted universe multiverse | tee -a /etc/apt/sources.list

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
               zip             \
               nano            \
               pkg-config      \
               python3-dev     \
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
               zlib1g-dev

RUN echo "Create directories" && \
    mkdir -p /opt /tmp /data /config

RUN echo "Install Tini - A tiny but valid init for containers" && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb

RUN echo "Cleanup" && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN TOOL="bowtie" && \
    echo "Download ${TOOL}" && \
    REPO="BenLangmead/bowtie" && \
    PAGE=$(curl -s https://api.github.com/repos/${REPO}/releases/latest) && \
    VER=$(echo "${PAGE}" | grep "tag_name" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    echo "Version ${VER} found" && \
    URL=$(echo "${PAGE}" | grep "browser_download_url" | grep "linux" | cut -d ":" -f 2,3 | tr -d '\"') && \
    ZIP=/tmp/${TOOL}.zip && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} downloaded" && \
    unzip ${ZIP} -d /opt/ && \
    rm -f ${ZIP} && \
    mv $(ls -d /opt/${TOOL}*) /opt/${TOOL} && \
    chmod a+x -R /opt/${TOOL} && \
    echo "${TOOL} installed"

RUN TOOL="bowtie2" && \
    echo "Download ${TOOL}" && \
    REPO="BenLangmead/bowtie2" && \
    PAGE=$(curl -s https://api.github.com/repos/${REPO}/releases/latest) && \
    VER=$(echo "${PAGE}" | grep "tag_name" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    echo "Version ${VER} found" && \
    URL=$(echo "${PAGE}" | grep "browser_download_url" | grep "linux" | grep -v "sra" | cut -d ":" -f 2,3 | tr -d '\"') && \
    ZIP=/tmp/${TOOL}.zip && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} downloaded" && \
    unzip ${ZIP} -d /opt/ && \
    rm -f ${ZIP} && \
    mv $(ls -d /opt/${TOOL}*) /opt/${TOOL} && \
    chmod a+x -R /opt/${TOOL} && \
    echo "${TOOL} installed"

RUN TOOL="samtools" && \
    echo "Download ${TOOL}" && \
    REPO="samtools/samtools" && \
    PAGE=$(curl -s https://api.github.com/repos/${REPO}/releases/latest) && \
    VER=$(echo "${PAGE}" | grep "tag_name" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    echo "Version ${VER} found" && \
    URL=$(echo "${PAGE}" | grep "browser_download_url" | grep "tar.bz2" | cut -d ":" -f 2,3 | tr -d '\"') && \
    ZIP=/tmp/${TOOL}.tar.bz2 && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} downloaded" && \
    tar xvjf ${ZIP} -C /opt/ && \
    rm -f ${ZIP} && \
    mv $(ls -d /opt/${TOOL}*) /opt/${TOOL} && \
    chmod a+x -R /opt/${TOOL} && \
    cd /opt/${TOOL} && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure && \
    make && \
    make install && \
    cd / && \
    rm -rf /opt/${TOOL} && \
    echo "${TOOL} installed"

RUN TOOL="bedtools2" && \
    echo "Download ${TOOL}" && \
    REPO="arq5x/bedtools2" && \
    PAGE=$(curl -s https://api.github.com/repos/${REPO}/releases/latest) && \
    VER=$(echo "${PAGE}" | grep "tag_name" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    echo "Version ${VER} found" && \
    URL=$(echo "${PAGE}" | grep "browser_download_url" | grep "tar.gz" | cut -d ":" -f 2,3 | tr -d '\"') && \
    ZIP=/tmp/${TOOL}.tar.gz && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} downloaded" && \
    tar xf ${ZIP} -C /opt/ && \
    rm -f ${ZIP} && \
    chmod a+x -R /opt/${TOOL} && \
    cd /opt/${TOOL} && \
    make && \
    make install && \
    cd / && \
    rm -rf /opt/${TOOL} && \
    echo "${TOOL} installed"

RUN echo "Install packages for Python 3" && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    ln -s /usr/bin/pip3 /usr/bin/pip && \
    pip3 install numpy pandas redis

# Create user docker with password docker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse docker && \
    echo `echo "docker\ndocker\n" | passwd docker` && \
    chown docker:docker /data /config  && \
    chmod -R 755 /data /config

# Change user (CLI: su - docker)
USER docker

# Get the pipeline scripts
RUN mkdir -p /home/docker/scripts && \
    cd /home/docker/scripts && \
    git clone --recurse-submodules https://github.com/ivasilyev/bowtie-tools.git && \
    git clone https://github.com/ivasilyev/biopipelines-docker.git && \
    mv bowtie-tools/nBee.py bowtie-tools/cook_the_reference.py bowtie-tools/sam2coverage.py bowtie-tools/verify_coverages.py bowtie-tools/modules biopipelines-docker/bwt_filtering_pipeline/worker/worker.py . && \
    rm -rf bowtie-tools biopipelines-docker && \
    chmod a+x -R /home/docker/scripts

# Update Path variables
ENV PATH=$PATH:/opt/bowtie:/opt/bowtie2:/home/docker/bin:/home/docker/scripts
ENV HOME=/home/docker

# Directories to mount
VOLUME ["/data", "/config", "/reference"]

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["python3", "/home/docker/scripts/pipeline_wrapper.py"]

# Setup the default directory
WORKDIR /home/docker

# MANUAL BUILD COMMAND:
# docker build -t bwt_filtering_pipeline_worker . && docker tag bwt_filtering_pipeline_worker ivasilyev/bwt_filtering_pipeline_worker:latest && docker push ivasilyev/bwt_filtering_pipeline_worker:latest
