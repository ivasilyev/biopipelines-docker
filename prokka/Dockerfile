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
            git \
            locales \
            parallel \
            python3-dev \
            python3-pip \
            tar \
            unzip \
            libpango1.0-dev

RUN echo "Cleanup" && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN echo "Reconfigure locale" && \
    locale-gen --purge en_US.UTF-8 && \
    DEBIAN_FRONTEND="noninteractive" dpkg-reconfigure locales && \
    update-locale LANG=en_US.UTF-8 LANGUAGE=en_US.UTF-8 LC_ALL=en_US.UTF-8

RUN TOOL="barrnap" && \
    echo "Download ${TOOL}" && \
    REPO="tseemann/barrnap" && \
    PAGE=$(curl -s https://api.github.com/repos/${REPO}/releases/latest) && \
    VER=$(echo "${PAGE}" | grep "tag_name" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    echo "Version ${VER} found" && \
    URL=$(echo "${PAGE}" | grep "zipball_url" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    ZIP=/tmp/${TOOL}.zip && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} was downloaded" && \
    unzip ${ZIP} -d /opt/ && \
    rm -f ${ZIP} && \
    mv $(ls -d /opt/*${TOOL}*) /opt/${TOOL} && \
    chmod a+x -R /opt/${TOOL} && \
    ln -s "/opt/${TOOL}/bin/barrnap" /usr/local/bin/barrnap && \
    echo "${TOOL} was placed into \"/opt/${TOOL}\""

RUN TOOL="infernal" && \
    echo "Download ${TOOL}" && \
    URL1="http://eddylab.org/infernal" && \
    PAGE=$(curl -fsSL ${URL1}) && \
    URL2=$(echo "${PAGE}" | grep -E 'linux.*\.tar\.gz' | cut -d '"' -f 2) && \
    ZIP=/tmp/${TOOL}.tar.gz && \
    curl -fsSL ${URL2} -o ${ZIP} && \
    echo "${TOOL} was downloaded" && \
    tar xf ${ZIP} -C /opt/ && \
    rm -f ${ZIP} && \
    mv $(ls -d /opt/*${TOOL}*) /opt/${TOOL} && \
    chmod a+x -R /opt/${TOOL} && \
    cd /opt/${TOOL} && \
    ./configure && \
    make && \
    make install && \
    cd easel && \
    make install && \
    cd / && \
    rm -rf /opt/${TOOL} && \
    echo "${TOOL} was installed"

RUN TOOL="genometools" && \
    echo "Download ${TOOL}" && \
    URL1="http://genometools.org/pub/" && \
    PAGE=$(curl -fsSL ${URL1}) && \
    URL2=${URL1}$(echo "${PAGE}" | grep -E 'genometools-[0-9].*\.tar\.gz' | \
                 grep -o ">genometools.*\.tar\.gz<" | tr -d '<>' | sort | tail -n 1) && \
    ZIP=/tmp/${TOOL}.tar.gz && \
    curl -fsSL ${URL2} -o ${ZIP} && \
    echo "${TOOL} was downloaded" && \
    tar xf ${ZIP} -C /opt/ && \
    rm -f ${ZIP} && \
    mv $(ls -d /opt/*${TOOL}*) /opt/${TOOL} && \
    chmod a+x -R /opt/${TOOL} && \
    cd /opt/${TOOL} && \
    make && \
    make install && \
    cd gtpython && \
    python setup.py install && \
    cd / && \
    rm -rf /opt/${TOOL} && \
    echo "${TOOL} was installed"

# make: *** [Makefile:811: obj/src/extended/gff3_escaping.o] Error 1 ;_;

# Create user docker with password docker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse docker && \
    echo `echo "docker\ndocker\n" | passwd docker`

# Change user (CLI: su - docker)
USER docker

# Update Path variables
ENV HOME=/home/docker
ENV PATH=$PATH:/opt/infernal/binaries

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["bash"]

# Setup the default directory
WORKDIR /home/docker
