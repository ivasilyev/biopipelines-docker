FROM ubuntu:19.04
# Base image CLI:
# IMG=ubuntu:19.04 && docker pull ${IMG} && docker run --rm -it --net=host ${IMG} bash

ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN echo "Software update" && \
    apt-get update -y && \
    apt-get install -y \
            software-properties-common \
            apt-utils && \
    add-apt-repository universe && \
    add-apt-repository multiverse

RUN echo "Software install" && \
    apt-get clean all && \
    apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y  \
            build-essential \
            curl \
            fuse \
            nano \
            python3-dev \
            python3-pip \
            unzip

RUN TOOL="tini" && \
    echo "Download ${TOOL}" && \
    REPO="krallin/tini" && \
    PAGE=$(curl -fsSL https://api.github.com/repos/${REPO}/releases/latest) && \
    VER=$(echo "${PAGE}" | grep "tag_name" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    echo "Version ${VER} was found" && \
    URL=$(echo "${PAGE}" | grep "browser_download_url" | grep -Po "\"http.+\.deb\"" | sed 's/\"//g' | sort | tail -n 1) && \
    ZIP=/tmp/${TOOL}.deb && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} was downloaded" && \
    dpkg -i ${ZIP} && \
    rm -f ${ZIP} && \
    echo "${TOOL} was installed"

RUN TOOL="miniconda" && \
    echo "Download ${TOOL}" && \
    _SH=/tmp/${TOOL}.sh && \
    curl -fsSL "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" -o ${_SH} && \
    echo 'export PATH=${PATH}:/opt/conda/bin' > /etc/profile.d/conda.sh && \
    bash ${_SH} -b -p "/opt/conda" && \
    rm -f ${_SH} && \
    chmod 777 -R /opt/conda && \
    echo "${TOOL} was installed"

RUN echo "Cleanup" && \
    apt-get autoremove --purge -y && \
    apt-get clean all && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV PATH=${PATH}:/opt/conda/bin

RUN conda config --add channels bioconda && \
    conda upgrade -y conda

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["bash"]

# Setup the default directory
WORKDIR /
