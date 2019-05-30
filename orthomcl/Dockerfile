FROM ubuntu:19.04
# Base image CLI:
# IMG=ubuntu:19.04 && docker pull ${IMG} && docker run --rm --net=host -it ${IMG} bash

ENV DEBIAN_FRONTEND noninteractive

RUN echo "Software update" && \
    apt-get -y update && \
    apt-get -y install \
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
            mysql-server \
            libmysqlclient-dev \
            perl \
            python3 \
            python3-pip

RUN echo "Cleanup" && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN TOOL="tini" && \
    echo "Download ${TOOL}" && \
    REPO="krallin/tini" && \
    PAGE=$(curl -s https://api.github.com/repos/${REPO}/releases/latest) && \
    VER=$(echo "${PAGE}" | grep "tag_name" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    echo "Version ${VER} was found" && \
    URL=$(echo "${PAGE}" | grep "browser_download_url" | grep -Po "\"http.+\.deb\"" | sed 's/\"//g' | sort | tail -n 1) && \
    ZIP=/tmp/${TOOL}.deb && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} was downloaded" && \
    dpkg -i ${ZIP} && \
    rm -f ${ZIP} && \
    echo "${TOOL} was installed"

RUN TOOL="ncbi-blast" && \
    echo "Download ${TOOL}" && \
    URL1="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/" && \
    PAGE=$(curl -fsSL ${URL1}) && \
    URL2=${URL1}$(echo "${PAGE}" | grep "src.tar.gz$" | awk '{print $NF}') && \
    ZIP=/tmp/${TOOL}.tar.gz && \
    curl -fsSL ${URL2} -o ${ZIP} && \
    echo "${TOOL} was downloaded" && \
    tar xf ${ZIP} -C /opt/ && \
    rm -f ${ZIP} && \
    mv $(ls -d /opt/${TOOL}*) /opt/${TOOL} && \
    chmod a+x -R /opt/${TOOL} && \
    cd /opt/${TOOL}/c++ && \
    ./configure && \
    cd /opt/${TOOL}/c++/ReleaseMT/build && \
    make all_r && \
    cp -r /opt/${TOOL}/c++/ReleaseMT/bin/* /opt/${TOOL}/ && \
    cd / && \
    rm -rf /opt/${TOOL}/c++ && \
    echo "${TOOL} was installed"

RUN TOOL="mcl" && \
    echo "Download ${TOOL}" && \
    URL="https://micans.org/mcl/src/mcl-latest.tar.gz" && \
    ZIP=/tmp/${TOOL}.tar.gz && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} was downloaded" && \
    tar xf ${ZIP} -C /opt/ && \
    rm -f ${ZIP} && \
    mv $(ls -d /opt/${TOOL}*) /opt/${TOOL} && \
    chmod a+x -R /opt/${TOOL} && \
    cd /opt/${TOOL} && \
    ./configure --enable-blast --prefix=/opt/${TOOL} && \
    make && \
    make install \ &&
    cd / && \
    rm -rf /opt/${TOOL} && \
    echo "${TOOL} was installed"

RUN TOOL="diamond" && \
    echo "Download ${TOOL}" && \
    REPO="bbuchfink/diamond" && \
    PAGE=$(curl -s https://api.github.com/repos/${REPO}/releases/latest) && \
    VER=$(echo "${PAGE}" | grep "tag_name" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    echo "Version ${VER} was found" && \
    URL=$(echo "${PAGE}" | grep "browser_download_url" | grep "linux64.tar.gz" | awk '{print $NF}' | tr -d '\", ') && \
    ZIP=/tmp/${TOOL}.tar.gz && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} was downloaded" && \
    mkdir -p /opt/${TOOL} && \
    tar xf ${ZIP} -C /opt/${TOOL} && \
    rm -f ${ZIP} && \
    chmod a+x -R /opt/${TOOL} && \
    ln -s "/opt/${TOOL}/fastqc" /usr/local/bin/diamond && \
    echo "${TOOL} was placed into \"/opt/${TOOL}\""


# Create user docker with password docker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups fuse docker && \
    usermod -aG wheel docker && \
    echo -e "docker\ndocker" | passwd docker

# Change user (CLI: su - docker)
USER docker

# Update Path variables
ENV HOME=/home/docker
ENV PATH=$PATH:/opt/ncbi-blast

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["bash"]

# Setup the default directory
WORKDIR /home/docker





RUN echo "Install NCBI BLAST" && \
    mkdir -p /opt /tmp && \
    curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz -o /tmp/ncbi-blast.tar.gz && \
    tar -xzvf /tmp/ncbi-blast.tar.gz -C /opt/ && \
    rm -f /tmp/ncbi-blast.tar.gz && \
    mv $(ls -d /opt/* | grep ncbi-blast) /opt/ncbi-blast

RUN echo "Install MCL" && \
    curl -fsSL https://micans.org/mcl/src/mcl-latest.tar.gz -o /tmp/mcl.tar.gz && \
    tar -xzvf /tmp/mcl.tar.gz -C /tmp/ && \
    rm -f /tmp/mcl.tar.gz && \
    mv $(ls -d /tmp/* | grep mcl) /tmp/mcl && \
    cd /tmp/mcl && \
    ./configure --enable-blast --prefix=/tmp/mcl \ &&
    make install \ &&
    cd .. \ &&
    rm -rf /tmp/mcl

RUN echo "Install OrthoMCL"

# Install some utilities
RUN yum install -y wget tar perl perl-DBI perl-DBD-MySQL gcc make mysql

# Make a software folder
RUN mkdir /software

# Install blast
RUN cd /software \
    && wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/blast-2.2.26-x64-linux.tar.gz \
    && tar -xzvf blast-2.2.26-x64-linux.tar.gz \
    && rm blast-2.2.26-x64-linux.tar.gz

# Install MCL
RUN cd /software \
    && wget http://micans.org/mcl/src/mcl-14-137.tar.gz \
    && tar -xzvf mcl-14-137.tar.gz \
    && rm mcl-14-137.tar.gz \
    && cd /software/mcl-14-137 \
    && ./configure --enable-blast --prefix=/software/mcl \
    && make install \
    && cd /software \
    && rm -rf /software/mcl-14-137

# Install OrthoMCL
RUN cd /software \
    && wget http://www.orthomcl.org/common/downloads/software/v2.0/orthomclSoftware-v2.0.9.tar.gz \
    && tar -xzvf orthomclSoftware-v2.0.9.tar.gz \
    && rm orthomclSoftware-v2.0.9.tar.gz

# Setup the paths
ENV PATH $PATH:/software/blast-2.2.26/bin
ENV PATH $PATH:/software/mcl/bin
ENV PATH $PATH:/software/orthomclSoftware-v2.0.9/bin
