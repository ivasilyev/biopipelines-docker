# Base image CLI:
# IMG=ubuntu:18.04 && docker pull ${IMG} && docker run --rm -it --net=host ${IMG} bash
FROM ubuntu:18.04

# Enable apt-get noninteractive install
ENV DEBIAN_FRONTEND noninteractive

RUN echo "Software update" && \
    apt-get -y update && \
    apt-get -y install \
            software-properties-common \
            apt-utils && \
    add-apt-repository universe && \
    add-apt-repository multiverse

RUN apt-get clean all && \
    apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y  \
        autotools-dev   \
        automake        \
        cmake           \
        curl            \
        grep            \
        sed             \
        dpkg            \
        fuse            \
        git             \
        wget            \
        zip             \
        openjdk-8-jre   \
        build-essential \
        pkg-config      \
        python3         \
        python3-pip     \
        bzip2           \
        ca-certificates \
        libglib2.0-0    \
        libxext6        \
        libsm6          \
        libxrender1     \
        git             \
        mercurial       \
        subversion      \
        zlib1g-dev      \
        fonts-dejavu    \
        gfortran        \
        gcc

# Install Tini - A tiny but valid init for containers
RUN TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -fsSL "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" -o tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

# Cleanup
RUN apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Python packages
RUN pip3 install --upgrade pip
pip3 install setuptools
# Base data analysis
pip3 install pandas statsmodels scipy xlrd
# Web scraping
pip3 install  requests bs4 lxml jinja2 pyyaml
# Visualization
pip3 install matplotlib seaborn
# External access
pip3 install ipython jupyter jupyterlab
# ML-related
pip3 install scikit-learn tensorflow

# Install MiniConda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add user docker with password docker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse docker && \
    echo `echo "docker\ndocker\n" | passwd docker`

# Give write permissions to conda folder
RUN chmod 777 -R /opt/conda/

# Change user
USER docker

ENV PATH=$PATH:/opt/conda/bin
ENV PATH=$PATH:/home/docker/bin
ENV HOME=/home/docker

RUN mkdir /home/docker/bin

RUN conda config --add channels r bioconda && \
    conda upgrade -y conda


# Install additional Conda packages
RUN for c in "r-essentials" \
            "r-base=3.4.1" \
            "r-irkernel=0.8*" \
            "r-plyr=1.8*" \
            "r-devtools=1.13*" \
            "r-tidyverse=1.1*" \
            "r-shiny=1.0*" \
            "r-rmarkdown=1.8*" \
            "r-forecast=8.2*" \
            "r-rsqlite=2.0*" \
            "r-reshape2=1.4*" \
            "r-nycflights13=0.2*" \
            "r-caret=6.0*" \
            "r-rcurl=1.95*" \
            "r-crayon=1.3*" \
            "r-randomforest=4.6*" \
            "r-htmltools=0.3*" \
            "r-sparklyr=0.7*" \
            "r-htmlwidgets=1.0*" \
            "r-hexbin=1.27*"
    do conda install --quiet -y ${c} && conda upgrade -y -c ${c}
    done
    conda clean -tipsy && \
    fix-permissions /opt/conda/

RUN conda upgrade conda

RUN echo "Prepare R Essentials Environment In Jupyter" && \
    conda create -n my-r-env -c r r-essentials && \
    conda skeleton cran ldavis conda build r-ldavis

CMD ["/bin/bash"]

WORKDIR /data