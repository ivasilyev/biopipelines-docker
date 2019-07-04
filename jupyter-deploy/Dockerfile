# Base image CLI:
# IMG=ubuntu:18.04 && docker pull ${IMG} && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ${IMG} bash
FROM ubuntu:18.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update && \
    apt-get -y upgrade && \
    apt-get -y autoremove && \
    apt-get -y install \
            # Base packages
            ca-certificates \
            curl \
            git \
            nano \
            # Python
            python3-dev \
            python3-pip \
            # LXML dependenciea
            libxml2-dev \
            libxslt1-dev \
            zlib1g-dev \
            # Required by matplotlib
            python3-tk \
            # NodeJS backend
            nodejs \
            npm

# Install Tini - A tiny but valid init for containers
RUN TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb

# Cleanup
RUN apt-get -y autoremove && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Python packages
RUN pip3 install --upgrade pip
RUN pip3 install setuptools \
    # Base data analysis
    pandas statsmodels scipy xlrd \
    # Sequence data processing
    biopython \
    # Web scraping
    requests bs4 lxml jinja2 pyyaml \
    # Visualization
    matplotlib seaborn \
    # Web interface
    ipython jupyter jupyterlab \
    # ML-related
    scikit-learn tensorflow

# Web port forwarding
EXPOSE 80 443 31522

# Create user docker with password docker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse docker && \
    echo `echo "docker\ndocker\n" | passwd docker` && \
    mkdir -p /data /config /home/docker/bin && \
    chown docker:docker /data /config /home/docker/bin && \
    chmod -R 755 /data /config /home/docker/bin

# Change user (CLI: su - docker)
USER docker

# Get the automation script for Jupyter Lab instances
RUN mkdir -p /home/docker/scripts && \
    curl -fsSL "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/jupyter-deploy/jupyter-deploy.sh" -o /home/docker/scripts/jupyter-deploy.sh

# Update environment variables
ENV PATH=$PATH:/home/docker/bin
ENV HOME=/home/docker
WORKDIR /home/docker

CMD ["/bin/bash"]

# MANUAL BUILD COMMAND:
# export DOCKER_IMAGE_NAME=jupyter-deploy && docker build -t ${DOCKER_IMAGE_NAME} . && docker tag ${DOCKER_IMAGE_NAME} ivasilyev/${DOCKER_IMAGE_NAME}:latest && docker push ivasilyev/${DOCKER_IMAGE_NAME}:latest
