FROM ubuntu:20.04

# BASE IMAGE ACCESSING COMMAND:
# export IMG="ubuntu:20.04" && docker pull "${IMG}" && docker run --rm --net=host -it "${IMG}"

ENV DEBIAN_FRONTEND noninteractive

RUN echo "Update packages" && \
    apt-get -y update && \
    apt-get -y upgrade && \
    apt-get -y autoremove && \
    echo "Install base packages" && \
    apt-get -y install \
            apt-utils \
            ca-certificates \
            software-properties-common \
            && \
    echo "Install main packages" && \
    apt-get -y install \
            # Base packages
            curl \
            git \
            nano \
            # Unicode stuff
            catdoc \
            convmv \
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
            npm \
            # ODBC header files
            unixodbc-dev \
            && \
    echo "Cleanup" && \
    apt-get -y autoremove && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN echo "Install Tini - A tiny but valid init for containers" && \
    TINI_VERSION=`curl -s "https://github.com/krallin/tini/releases/latest" | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    DL_FILE="/tmp/tini.deb" && \
    mkdir -p "(basename "${DL_FILE}")" && \
    curl -fsSL "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" -o "${DL_FILE}" && \
    dpkg -i "${DL_FILE}" && \
    rm -f "${DL_FILE}"

RUN echo "Install primary Python packages" && \
    cd /tmp/ && \
    pip3 install --upgrade --no-cache-dir pip setuptools wheel && \
    curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/jupyter-deploy/requirements.txt" && \
    curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/jupyter-deploy/requirements-linux.txt" && \
    pip3 install --no-cache-dir --requirement "requirements.txt" --requirement "requirements-linux.txt" && \
    echo "Clean Python package cache" && \
    pip cache purge && \
    rm -rf /tmp/*

# Web port forwarding
EXPOSE 80 443 31522

RUN echo "Create user docker with password docker" && \
    groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse docker && \
    echo `echo "docker\ndocker\n" | passwd docker` && \
    mkdir -p /config /home/docker/bin && \
    chown docker:docker /config /home/docker/bin && \
    chmod -R a+rw /config /home/docker/bin

# Change user (CLI: su - docker)
USER docker

RUN echo "Get the automation script for Jupyter Lab instances" && \
    mkdir -p /home/docker/scripts && \
    curl -fsSL "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/jupyter-deploy/jupyter-deploy.sh" -o /home/docker/scripts/jupyter-deploy.sh

# Update environment variables
ENV PATH=$PATH:/home/docker/bin
ENV HOME=/home/docker
WORKDIR /home/docker

CMD ["/bin/bash"]

# MANUAL BUILD COMMAND:
# export REPO="jupyter-deploy" && export TAG="ivasilyev/${REPO}:latest" && docker build --network=host --tag "${REPO}" . && docker tag "${REPO}" "${TAG}" && docker push "${TAG}" && rm -f Dockerfile

# PULL COMMAND:
# export IMG="ivasilyev/jupyter-deploy:latest" && docker pull "${IMG}" && docker run --rm --net=host -it "${IMG}"
