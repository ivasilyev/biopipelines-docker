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

RUN echo "Install packages for Python 3" && \
    ln -sfn /usr/bin/python3 /usr/bin/python && \
    ln -sfn /usr/bin/pip3 /usr/bin/pip && \
    pip3 install redis

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

RUN TOOL="bcl2fastq" && \
    echo "Download ${TOOL}" && \
    PAGE=$(curl -s https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html) && \
    VER=$(echo "${PAGE}" | grep "Linux tarball" | grep -Po 'v[^ a-z]+') && \
    echo "Version ${VER} found" && \
    URL=$(echo "${PAGE}" | sed -n 's/.*href="\([^"]*\).*/\1/p' | grep "tar.zip") && \
    ZIP1=/tmp/${TOOL}.tar.zip && \
    curl -fsSL ${URL} -o ${ZIP1} && \
    unzip ${ZIP1} -d /tmp/ && \
    rm -f ${ZIP1} && \
    ZIP2=$(ls -d /tmp/${TOOL}*) && \
    tar xvf ${ZIP2} -C /opt/ && \
    rm -f ${ZIP2} && \
    chmod a+x -R /opt/${TOOL} && \
    cd /opt/${TOOL} && \
    ./src/configure --prefix=/opt/${TOOL} && \
    make && \
    make install && \
    cd / && \
    rm -rf /opt/${TOOL} && \
    echo "${TOOL} installed"

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
ENV PATH=$PATH:/opt/bowtie:/opt/bowtie2:/home/docker/bin:/home/docker/scripts
ENV HOME=/home/docker

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["bash"]

# Setup the default directory
WORKDIR /home/docker
