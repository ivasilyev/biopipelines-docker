FROM ubuntu:18.04
# CLI:
# export IMG=ubuntu:18.04 && docker pull ${IMG} && docker run --rm --net=host -it ${IMG} bash

RUN echo Install packages && \
    apt-get -qq update && apt-get -y upgrade && \
	apt-get install -y apt-utils python3-pip python3-dev build-essential curl bzip2 libz-dev

RUN echo "Install packages for Python 3" && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    ln -s /usr/bin/pip3 /usr/bin/pip && \
    pip3 install regex cutadapt

RUN echo "Cleanup" && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN TOOL="SPAdes" && \
    echo "Download ${TOOL}" && \
    REPO="ablab/spades" && \
    PAGE=$(curl -s https://api.github.com/repos/${REPO}/releases/latest) && \
    VER=$(echo "${PAGE}" | grep "tag_name" | cut -d ":" -f 2,3 | tr -d '\", ') && \
    echo "Version ${VER} found" && \
    URL=$(echo "${PAGE}" | grep "browser_download_url" | grep "Linux" | cut -d ":" -f 2,3 | tr -d '\" ') && \
    mkdir -p /opt /tmp && \
    ZIP=/tmp/${TOOL}.tar.gz && \
    curl -fsSL ${URL} -o ${ZIP} && \
    echo "${TOOL} downloaded" && \
    tar -xvzf ${ZIP} -C /opt/ && \
    rm -f ${ZIP} && \
    mv $(ls -d /opt/${TOOL}*) /opt/${TOOL} && \
    chmod a+x -R /opt/${TOOL} && \
    echo "${TOOL} installed"

RUN echo "Install Tini - A tiny but valid init for containers" && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb

# Create user docker with password docker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse docker && \
    echo `echo "docker\ndocker\n" | passwd docker` && \
    mkdir -p /data /config /home/docker/bin && \
    chown docker:docker /data /config /home/docker/bin && \
    chmod -R 755 /data /config /home/docker/bin

# Change user (CLI: su - docker)
USER docker

# Update PATH
ENV PATH=/opt/SPAdes/bin/:${PATH}
ENV HOME=/home/docker
WORKDIR /home/docker
