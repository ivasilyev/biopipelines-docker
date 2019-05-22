FROM ubuntu:18.04
# Base image CLI: docker run --rm --net=host -it ubuntu:18.04 bash
# In the case of network failure, try:
# echo deb http://mirror.truenetwork.ru/ubuntu/ bionic-backports main restricted universe multiverse | tee -a /etc/apt/sources.list

ENV DEBIAN_FRONTEND noninteractive

RUN echo "Software update" && \
    apt-get -y update && \
    apt-get -y install software-properties-common apt-utils && \
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
    mkdir -p /opt /tmp

RUN TOOL="tini" && \
    echo "Download ${TOOL}" && \
    REPO="krallin/tini" && \
    PAGE=$(curl -s https://api.github.com/repos/${REPO}/releases/latest) && \
    VER=$(echo "${PAGE}" | grep "tag_name" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    echo "Version ${VER} found" && \
    URL=$(echo "${PAGE}" | grep "browser_download_url" | grep -Po "\"http.+\.deb\"" | sed 's/\"//g' | sort | tail -n 1) && \
    ZIP=/tmp/${TOOL}.deb && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} downloaded" && \
    dpkg -i ${ZIP} && \
    rm -f ${ZIP} && \
    echo "${TOOL} installed"

# .NET is required by Illumina Canvas, details: https://github.com/Illumina/canvas/issues/99
RUN echo "Install M$ .NET Core SDK x64" && \
    apt-get install -y apt-transport-https && \
    ZIP=/tmp/packages-microsoft-prod.deb && \
    curl -fsSL https://packages.microsoft.com/config/ubuntu/18.04/packages-microsoft-prod.deb -o ${ZIP} && \
    dpkg -i ${ZIP} && \
    rm -f ${ZIP} && \
    apt-get update -y && \
    apt-get install -y dotnet-sdk-2.1

RUN echo "Install Mono" && \
    apt-get install -y gnupg ca-certificates && \
    apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF && \
    echo "deb https://download.mono-project.com/repo/ubuntu stable-bionic main" | tee /etc/apt/sources.list.d/mono-official-stable.list && \
    apt-get update -y && \
    apt-get install -y mono-devel

RUN echo "Cleanup" && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

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

RUN TOOL="bcftools" && \
    echo "Download ${TOOL}" && \
    REPO="samtools/bcftools" && \
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

RUN TOOL="vcftools" && \
    echo "Download ${TOOL}" && \
    REPO="vcftools/vcftools" && \
    PAGE=$(curl -s https://api.github.com/repos/${REPO}/releases/latest) && \
    VER=$(echo "${PAGE}" | grep "tag_name" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    echo "Version ${VER} found" && \
    URL=$(echo "${PAGE}" | grep "browser_download_url" | grep "tar.gz" | cut -d ":" -f 2,3 | tr -d '\"') && \
    ZIP=/tmp/${TOOL}.tar.gz && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} downloaded" && \
    tar xf ${ZIP} -C /opt/ && \
    rm -f ${ZIP} && \
    mv $(ls -d /opt/${TOOL}*) /opt/${TOOL} && \
    chmod a+x -R /opt/${TOOL} && \
    cd /opt/${TOOL} && \
    ./configure && \
    make && \
    make install && \
    cd / && \
    rm -rf /opt/${TOOL} && \
    echo "${TOOL} installed"

RUN TOOL="Canvas" && \
    echo "Download ${TOOL}" && \
    REPO="Illumina/canvas" && \
    PAGE=$(curl -s https://api.github.com/repos/${REPO}/releases/latest) && \
    # VER=$(echo "${PAGE}" | grep "tag_name" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    VER="1.35.1.1316+master" && \
    echo "Version ${VER} found" && \
    # URL=$(echo "${PAGE}" | grep "browser_download_url" | grep "tar.gz" | cut -d ":" -f 2,3 | tr -d '\"') && \
    URL="https://github.com/Illumina/canvas/releases/download/1.35.1.1316%2Bmaster/Canvas-1.35.1.1316.master_x64.tar.gz" && \
    ZIP=/tmp/${TOOL}.tar.gz && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} downloaded" && \
    tar xf ${ZIP} -C /opt/ && \
    # The archive also contains the redundant empty folder
    rm -rf ${ZIP} /opt/${VER} && \
    mv $(ls -d /opt/${TOOL}*) /opt/${TOOL} && \
    chmod a+x -R /opt/${TOOL}

RUN echo "Install packages for Python 3" && \
    ln -sfn /usr/bin/python3 /usr/bin/python && \
    ln -sfn /usr/bin/pip3 /usr/bin/pip && \
    pip3 install redis

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
ENV PATH=$PATH:/opt/Canvas
ENV HOME=/home/docker

# Directories to mount
VOLUME ["/data", "/config", "/reference"]

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["bash"]

# Setup the default directory
WORKDIR /home/docker
