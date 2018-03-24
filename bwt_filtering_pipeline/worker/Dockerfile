################## BASE IMAGE ######################
FROM ubuntu:17.10

################## METADATA ######################
LABEL base.image="ubuntu:17.10"
LABEL pipeline="bwt_filtering"
LABEL name="bwt_filtering_pipeline_worker"
LABEL description="Docker image containing all required packages to run filtering Bowtie-based pipeline"
LABEL website="https://github.com/ivasilyev/bowtie-tools"
LABEL documentation="https://github.com/ivasilyev/ThreeBees"
LABEL license="GPL-2.0"
LABEL tags="Genomics"
LABEL maintainer="Ilya Vasilyev <u0412u0418u042e@gmail.com>"

################## INSTALLATION ######################
# Base image CLI: docker run --rm -it ubuntu:17.10 bash
# Declare variables
ENV REF_DIR="/reference"

# Enable apt-get noninteractive install
ENV DEBIAN_FRONTEND noninteractive

# Enable repositories - disabling test
#RUN mv /etc/apt/sources.list /etc/apt/sources.list.bak && \
#    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse\n\
#                      deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse\n\
#                      deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-backports main restricted universe multiverse\n\
#                      deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
#    cat /etc/apt/sources.list.bak >> /etc/apt/sources.list && \
#    cat /etc/apt/sources.list

RUN apt-get -y update && \
    apt-get -y install \
            software-properties-common \
            apt-utils && \
    add-apt-repository universe && \
    add-apt-repository multiverse

# Software update
RUN apt-get clean all && \
    apt-get -y update && \
    apt-get -y upgrade && \
    apt-get -y install \
               autotools-dev   \
               automake        \
               autoconf        \
               build-essential \
               cmake           \
               flex            \
               curl            \
               grep            \
               sed             \
               dpkg            \
               fuse            \
               git             \
               wget            \
               zip             \
               nano            \
               openjdk-8-jre   \
               pkg-config      \
               python          \
               python-dev      \
               python-pip      \
               python3         \
               python3-pip     \
               bzip2           \
               ca-certificates \
               lib32ncurses5-dev \
               libatlas-base-dev \
               libblas-dev     \
               libboost-all-dev \
               libbz2-dev      \
               libcurl4-openssl-dev \
               libev-dev       \
               libglib2.0-0    \
               libhdf5-dev     \
               liblapack-dev   \
               liblzma-dev     \
               libsm6          \
               libssl-dev      \
               libtbb2         \
               libtbb-dev      \
               libxext6        \
               libxrender1     \
               mercurial       \
               subversion      \
               zlib1g-dev

# Install Tini - A tiny but valid init for containers
RUN TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb

# Cleanup
RUN apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Compile Bowtie - crashes
#RUN mkdir -p /tmp/bowtie && \
#    cd /tmp/bowtie && \
#    export BOWTIE_VERSION=1.2.2_p1 && \
#    curl -fsSL https://github.com/BenLangmead/bowtie/archive/v${BOWTIE_VERSION}.tar.gz -o bowtie.tar.gz && \
#    tar xf bowtie.tar.gz && \
#    cd bowtie-${BOWTIE_VERSION} && \
#    make

# Deploy bowtie
RUN cd /tmp && \
    wget -q https://github.com/BenLangmead/bowtie/releases/download/v1.2.2_p1/bowtie-1.2.2-linux-x86_64.zip && \
    unzip bowtie-1.2.2-linux-x86_64.zip && \
    cd bowtie-1.2.2-linux-x86_64 && \
    cp bowtie bowtie-align-l bowtie-align-s bowtie-build bowtie-build-l bowtie-build-s bowtie-inspect bowtie-inspect-l bowtie-inspect-s /usr/local/bin/ && \
    cd /usr/local/bin/ && \
    chmod ugo+x bowtie bowtie-align-l bowtie-align-s bowtie-build bowtie-build-l bowtie-build-s bowtie-inspect bowtie-inspect-l bowtie-inspect-s && \
    cd /tmp && \
    rm -rf bowtie-1.2.2-linux-x86_64 bowtie-1.2.2-linux-x86_64.zip

# Deploy bowtie2
RUN cd /tmp && \
    wget -q https://github.com/BenLangmead/bowtie2/releases/download/v2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip && \
    unzip bowtie2-2.3.4.1-linux-x86_64.zip && \
    cd bowtie2-2.3.4.1-linux-x86_64 && \
    cp bowtie2 bowtie2-align-l bowtie2-align-s bowtie2-build bowtie2-build-l bowtie2-build-s bowtie2-inspect bowtie2-inspect-l bowtie2-inspect-s /usr/local/bin/ && \
    cd /usr/local/bin/ && \
    chmod ugo+x bowtie2 bowtie2-align-l bowtie2-align-s bowtie2-build bowtie2-build-l bowtie2-build-s bowtie2-inspect bowtie2-inspect-l bowtie2-inspect-s && \
    cd /tmp && \
    rm -rf bowtie2-2.3.4.1-linux-x86_64 bowtie2-2.3.4.1-linux-x86_64.zip

# Compile samtools
RUN mkdir -p /tmp/samtools && \
    cd /tmp/samtools && \
    export SAMTOOLS_VERSION=1.7 && \
    curl -fsSL https://github.com/samtools/htslib/archive/${SAMTOOLS_VERSION}.tar.gz -o htslib.tar.gz && \
    tar xf htslib.tar.gz && \
    cd htslib-${SAMTOOLS_VERSION} && \
    autoheader && \
    autoconf && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    curl -fsSL https://github.com/samtools/samtools/archive/${SAMTOOLS_VERSION}.tar.gz -o samtools.tar.gz && \
    tar xf samtools.tar.gz && \
    cd samtools-${SAMTOOLS_VERSION}/ && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure && \
    make && \
    make install && \
    cd /tmp && \
    rm -rf samtools

# Compile bedtools
RUN mkdir -p /tmp/bedtools && \
    cd /tmp/bedtools && \
    export BEDTOOLS_VERSION=2.27.1 && \
    curl -fsSL https://github.com/arq5x/bedtools2/archive/v${BEDTOOLS_VERSION}.tar.gz -o bedtools.tar.gz && \
    tar xf bedtools.tar.gz && \
    cd bedtools2-${BEDTOOLS_VERSION} && \
    make && \
    make install && \
    cd /tmp && \
    rm -rf bedtools

# Install Python packages
RUN for p in numpy scipy h5py biom-format pandas redis matplotlib statsmodels; do pip install ${p}; pip3 install ${p}; done

# Create user docker with password docker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse docker && \
    echo `echo "docker\ndocker\n" | passwd docker` && \
    mkdir /data /config ${REF_DIR} && \
    chown docker:docker /data /config ${REF_DIR} && \
    chmod -R 755 /data /config ${REF_DIR}

# Change user (CLI: su - docker)
USER docker

# Update Path variables
ENV PATH=$PATH:/home/docker/bin
ENV HOME=/home/docker

# Get the pipeline scripts
RUN mkdir ~/scripts && \
    cd ~/scripts && \
    git clone --recurse-submodules https://github.com/ivasilyev/bowtie-tools.git && \
    git clone https://github.com/ivasilyev/biopipelines-docker.git && \
    mv bowtie-tools/nBee.py bowtie-tools/cook_the_reference.py bowtie-tools/sam2coverage.py bowtie-tools/verify_coverages.py biopipelines-docker/bwt_filtering_pipeline/worker/pipeline_wrapper.py biopipelines-docker/bwt_filtering_pipeline/worker/queue_handler.py . && \
    rm -rf bowtie-tools biopipelines-docker && \
    chmod -R +x ~/scripts && \
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
VOLUME ["/data", "/config", "/reference"]

# Overwrite this with 'CMD []' in a dependent Dockerfile
# CMD ["/bin/bash"]
CMD ["python3", "/home/docker/scripts/pipeline_wrapper.py"]

# Setup the default directory
WORKDIR /home/docker

# MANUAL BUILD COMMAND:
# docker build -t bwt_filtering_pipeline_worker . && docker tag bwt_filtering_pipeline_worker ivasilyev/bwt_filtering_pipeline_worker:latest && docker push ivasilyev/bwt_filtering_pipeline_worker:latest
