# Base image CLI:
# IMG=qiime2/core && docker pull ${IMG} && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ${IMG} bash
FROM qiime2/core

ENV DEBIAN_FRONTEND noninteractive

RUN echo "Software update" && \
    apt-get -y update && \
    apt-get -y upgrade && \
    apt-get -y autoremove

RUN echo "Software install" && \
    apt-get -y install software-properties-common apt-utils dpkg gnupg2 curl

RUN echo "Install Tini - A tiny but valid init for containers" && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb

RUN echo "Cleanup" && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN echo "Create user" && \
    groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse docker && \
    echo `echo "docker\ndocker\n" | passwd docker` && \
    mkdir -p /data /config /home/docker/bin && \
    chown docker:docker /data /config /home/docker/bin && \
    chmod -R 755 /data /config /home/docker/bin

# Change user (CLI: su - docker)
USER docker

# Update environment variables
ENV PATH=$PATH:/home/docker/bin
ENV HOME=/home/docker
WORKDIR /home/docker

CMD ["/bin/bash"]
