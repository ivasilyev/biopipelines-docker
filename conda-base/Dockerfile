FROM ubuntu:19.04
# Base image CLI:
# IMG=ubuntu:19.04 && docker pull ${IMG} && docker run --rm -it --net=host ${IMG} bash

ENV DEBIAN_FRONTEND noninteractive

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
            git \
            curl \
            openjdk-8-jre \
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
    echo "${TOOL} was installed"

RUN echo "Cleanup" && \
    apt-get autoremove --purge && \
    apt-get clean all && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Create user docker with password docker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups fuse docker && \
    usermod -aG wheel docker && \
    echo -e "docker\ndocker" | passwd docker

# Change user (CLI: su - docker)
USER docker

# Update Path variables
ENV HOME=/home/docker

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["bash"]

# Setup the default directory
WORKDIR /home/docker
