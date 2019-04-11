FROM ubuntu:18.04
# CLI:
# export IMG=ubuntu:18.04 && docker pull ${IMG} && docker run --rm --net=host -it ${IMG} bash

RUN echo Install packages && \
    apt-get -qq update && apt-get -y upgrade && \
	apt-get install -y python3-pip python3-dev build-essential curl bzip2 libz-dev

RUN pip3 install regex cutadapt

RUN echo Cleanup && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN echo Download SPAdes && \
    REPO="ablab/spades" && \
    VER=$(curl -s https://api.github.com/repos/${REPO}/releases/latest | \
    grep "tag_name" | \
    cut -d ":" -f 2,3 | \
    tr -d '\", ') && \
    \
    echo Version ${VER} found && \
    URL=$(curl -s https://api.github.com/repos/${REPO}/releases/latest | \
    grep "browser_download_url" | \
    grep "Linux" | \
    cut -d ":" -f 2,3 | \
    tr -d '\"') && \
    \
    mkdir -p /opt /tmp && \
    echo Download started && \
    curl -fsSL ${URL} -o /tmp/spades.tar.gz && \
    tar -xvzf /tmp/spades.tar.gz -C /opt/ && \
    rm -f /tmp/spades.tar.gz && \

RUN echo Install Tini - A tiny but valid init for containers && \
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
ENV PATH="$(ls -d /opt/* | grep SPAdes)/bin/:${PATH}"
ENV HOME=/home/docker
WORKDIR /home/docker
