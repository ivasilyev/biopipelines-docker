# Source: https://github.com/broadinstitute/gatk/blob/master/Dockerfile

# Using OpenJDK 8
FROM broadinstitute/gatk:latest

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update && \
    apt-get -y upgrade && \
    apt-get -y autoremove

# Install Tini - A tiny but valid init for containers
RUN TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb

# Cleanup
RUN apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Create user docker with password docker
RUN groupadd fuse && \
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
