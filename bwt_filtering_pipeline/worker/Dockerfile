################## BASE IMAGE ######################
FROM ubuntu:17.10

################## METADATA ######################
LABEL base.image="ubuntu:17.10"
LABEL version="1.4"
LABEL software="bwt_filtering_pipeline"
LABEL description="Docker image containing all required packages to run filtering Bowtie-based pipeline"
LABEL website="https://github.com/ivasilyev/bowtie-tools"
LABEL documentation="https://github.com/ivasilyev/ThreeBees"
LABEL license="GPL-2.0"
LABEL tags="Genomics"

################## MAINTAINER ######################
MAINTAINER Ilya Vasilyev <u0412u0418u042e@gmail.com>

################## INSTALLATION ######################
# Base image CLI: docker run --rm -it ubuntu:17.10 bash
# Declare variables
ENV MINICONDA_URL="https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh"
ENV REF_DIR="/reference"

# Enable apt-get noninteractive install
ENV DEBIAN_FRONTEND noninteractive

# Enable repositories
RUN mv /etc/apt/sources.list /etc/apt/sources.list.bak && \
    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse\n\
                      deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse\n\
                      deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-backports main restricted universe multiverse\n\
                      deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
    cat /etc/apt/sources.list.bak >> /etc/apt/sources.list && \
    cat /etc/apt/sources.list

# Software update
RUN apt-get clean all && \
    apt-get -y update && \
    apt-get -y upgrade

# Software installation
RUN apt-get -y install \
    autotools-dev   \
    automake        \
    build-essential \
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
    pkg-config      \
    python          \
	python-dev      \
    python-pip      \
    python3         \
    python3-pip     \
    bzip2           \
    ca-certificates \
    libglib2.0-0    \
    libtbb2         \
    libtbb-dev      \
    libxext6        \
    libsm6          \
    libxrender1     \
    libev-dev       \
    libblas-dev     \
    liblapack-dev   \
    libatlas-base-dev \
    libhdf5-dev     \
    mercurial       \
    subversion      \
    zlib1g-dev

# Install Tini - A tiny but valid init for containers
RUN TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

# Cleanup
RUN apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Conda - a cross-platform, language-agnostic binary package manager
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet ${MINICONDA_URL} -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Install Python packages
RUN pip3 install numpy \
                 scipy \
                 h5py \
                 biom-format \
                 pandas \
                 redis && \
    pip install numpy \
                 scipy \
                 h5py \
                 biom-format \
                 pandas \
                 redis


# Create user docker with password docker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse docker && \
    echo `echo "docker\ndocker\n" | passwd docker` && \
    mkdir /data /data1 /data2 /config ${REF_DIR} && \
    chown docker:docker /data /data1 /data2 /config ${REF_DIR}
    
# Give ultimate access to the Conda directory
RUN chmod 777 -R /opt/conda

# Change user (CLI: su - docker)
USER docker

# Update Path variables
ENV PATH=$PATH:/opt/conda/bin
ENV PATH=$PATH:/home/docker/bin
ENV HOME=/home/docker

# Setup Conda
RUN conda config --add channels r && \
    conda config --add channels bioconda && \
    conda config --add channels intel && \
    conda upgrade -y conda

# Install additional Conda packages 
RUN for s in "intel tbb" "bioconda bowtie=1.1.2" "bioconda bowtie2" "bioconda samtools" "bioconda bedtools"; do conda install -y -c ${s} && conda upgrade -y -c ${s}; done

# Get the pipeline scripts
RUN mkdir ~/scripts && \
    cd ~/scripts && \
    git clone --recurse-submodules https://github.com/ivasilyev/bowtie-tools.git && \
    git clone https://github.com/ivasilyev/biopipelines-docker.git && \
    mv bowtie-tools/nBee.py bowtie-tools/cook_the_reference.py bowtie-tools/sam2coverage.py biopipelines-docker/bwt_filtering_pipeline/worker/pipeline_wrapper.py biopipelines-docker/bwt_filtering_pipeline/worker/queue_handler.py . && \
    rm -rf bowtie-tools biopipelines-docker && \
    echo "export PATH=$PATH:${PWD}" | tee -a ~/.bashrc

# Human genome reference indexes, use if required
ENV HG19_COLORSPACE="ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19_c.ebwt.zip"
ENV HG19_BWT2I="ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip"
ENV HG37_COLORSPACE="ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/h_sapiens_37_asm_c.ebwt.zip"
ENV HG37_BWT2I="ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz"

# Download references - not implemented due to large size (~ 12 Gb)
# The hg19 index masks are 'hg19_c' and 'hg19'
# RUN mkdir ${REF_DIR}/hg19 && \
#    cd ${REF_DIR}/hg19 && \
#    wget --quiet ${HG19_COLORSPACE} && \
#    unzip $(ls ./*.zip) && \
#    rm $(ls ./*.zip) && \
#    wget --quiet ${HG19_BWT2I} && \
#    unzip $(ls ./*.zip) && \
#    rm $(ls ./*.zip) && \
#    printf "${REF_DIR}/hg19/hg19.fasta\t${REF_DIR}/hg19/hg19_c\t${REF_DIR}/hg19/hg19\tnone\tnone\tnone\n" > ${REF_DIR}/hg19/hg19.refdata

# The hg37 index masks are 'h_sapiens_37_asm_c' and 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'
# RUN mkdir ${REF_DIR}/hg37 && \
#    cd ${REF_DIR}/hg37 && \
#    wget --quiet ${HG37_COLORSPACE} && \
#    unzip $(ls ./*.zip) && \
#    rm $(ls ./*.zip) && \
#    wget --quiet ${HG37_BWT2I} && \
#    tar -xvzf $(ls ./*.tar.gz) && \
#    rm $(ls ./*.tar.gz) && \
#    printf "${REF_DIR}/hg37/hg37.fasta\t${REF_DIR}/hg37/h_sapiens_37_asm_c\t${REF_DIR}/hg37/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index\tnone\tnone\tnone\n" > ${REF_DIR}/hg37/hg37.refdata

# Directories to mount
VOLUME ["/data", "/data1", "/data2", "/config", "/reference"]

# Overwrite this with 'CMD []' in a dependent Dockerfile
# CMD ["/bin/bash"]
CMD ["python3", "pipeline_wrapper.py"]

# Setup the default directory
WORKDIR /home/docker

# Go to Dockerfile dir & run:
# docker build -t bwt_filtering_pipeline_worker .
# docker tag bwt_filtering_pipeline_worker ivasilyev/bwt_filtering_pipeline_worker:latest && docker push ivasilyev/bwt_filtering_pipeline_worker:latest
