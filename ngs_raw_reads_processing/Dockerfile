FROM centos:7.6.1810
# Base image CLI: docker run --rm --net=host -it centos:7.6.1810 bash

RUN echo "Software update" && \
    yum update -y && \
    yum install -y epel-release && \
    yum install -y https://centos7.iuscommunity.org/ius-release.rpm && \
    yum upgrade -y

RUN echo "Install packages" && \
    yum install -y \
        curl \
        java-1.8.0-openjdk-devel \
        git \
        nano \
        python36u-devel \
        python36u-pip \
        unzip \
        yum-utils && \
    yum groupinstall -y "Development tools"

RUN echo "Install packages for Python 3" && \
    PY_VER="3.6" && \
    ln -s  "/usr/bin/python${PY_VER}" "/usr/bin/python3" && \
    ln -s "/usr/bin/pip${PY_VER}" "/usr/bin/pip3" && \
    pip3 install --upgrade pip && \
    pip3 install redis

RUN echo "Create directories" && \
    mkdir -p /opt /tmp

RUN TOOL="tini" && \
    echo "Download ${TOOL}" && \
    REPO="krallin/tini" && \
    PAGE=$(curl -s https://api.github.com/repos/${REPO}/releases/latest) && \
    VER=$(echo "${PAGE}" | grep "tag_name" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    echo "Version ${VER} found" && \
    URL=$(echo "${PAGE}" | grep "browser_download_url" | grep -Po "\"http.+\.rpm\"" | sed 's/\"//g' | sort | tail -n 1) && \
    ZIP=/tmp/${TOOL}.rpm && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} downloaded" && \
    yum -y --nogpgcheck localinstall ${ZIP} && \
    rm -f ${ZIP} && \
    echo "${TOOL} installed"

RUN TOOL="bcl2fastq" && \
    echo "Download ${TOOL}" && \
    PAGE=$(curl -s https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html) && \
    VER=$(echo "${PAGE}" | grep "Linux tarball" | grep -Po 'v[^ a-z]+') && \
    echo "Version ${VER} found" && \
    URL="https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-20-0-linux-x86-64.zip" && \
    ZIP1=/tmp/${TOOL}.zip && \
    curl -fsSL ${URL} -o ${ZIP1} && \
    echo "${TOOL} downloaded" && \
    unzip ${ZIP1} -d /tmp/ && \
    rm -f ${ZIP1} && \
    ZIP2=$(ls -d /tmp/${TOOL}*) && \
    yum -y --nogpgcheck localinstall ${ZIP2} && \
    rm -f ${ZIP2} && \
    echo "${TOOL} installed"

RUN echo "Cleanup" && \
    yum clean all && \
    rm -rf /var/cache/yum/* /tmp/* /var/tmp/*

RUN TOOL="solid-denovo-assembly" && \
    echo "Download ${TOOL}" && \
    URL="https://www.biostars.org/static/downloads/solid/solid-denovo-assembly.tar.gz" && \
    ZIP=/tmp/${TOOL}.tar.gz && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} downloaded" && \
    tar xf ${ZIP} -C /opt/ && \
    rm -f ${ZIP} && \
    cd /opt/${TOOL} && \
    chmod a+x -R /opt/${TOOL} && \
    make && \
    cd / && \
    # TODO Add binaries to $PATH
    echo "${TOOL} is in \"/opt/${TOOL}\""

RUN TOOL="XSQ_Tools" && \
    echo "Download ${TOOL}" && \
    URL="https://www.biostars.org/static/downloads/solid/xsq.tools.tar.gz" && \
    ZIP=/tmp/${TOOL}.tar.gz && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} downloaded" && \
    tar xf ${ZIP} -C /opt/ && \
    rm -f ${ZIP} && \
    chmod a+x -R /opt/${TOOL} && \
    cd / && \
    # TODO Add binaries to $PATH
    echo "${TOOL} is in \"/opt/${TOOL}\""

# Create user docker with password docker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups fuse docker && \
    usermod -aG wheel docker && \
    echo `echo "docker\ndocker\n" | passwd docker`

# Change user (CLI: su - docker)
USER docker

# Update Path variables
ENV PATH=$PATH:/opt/bowtie:/opt/bowtie2:/home/docker/bin:/home/docker/scripts
ENV HOME=/home/docker

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["bash"]

# Setup the default directory
WORKDIR /home/docker
