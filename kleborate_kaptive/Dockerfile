FROM ivasilyev/conda-base:latest
# Base image CLI: docker run --rm --net=host -it ivasilyev/conda-base:latest bash

RUN echo "Software install" && \
    apt-get clean all && \
    apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y  \
            git

RUN echo "Cleanup" && \
    apt-get autoremove --purge -y && \
    apt-get clean all && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN echo "Software install" && \
    conda install --channel bioconda --quiet -y biopython kleborate kaptive && \
    cd /opt/ && \
    git clone --recursive https://github.com/katholt/Kleborate.git && \
    cd Kleborate && \
    python setup.py install && \
    mv /opt/conda/bin/kleborate /opt/conda/bin/_kleborate && \
    pip3 install biopython && \
    sed -i 's|^#!/usr/bin/env python$|#!/usr/bin/env python3|' kleborate-runner.py && \
    ln -s $(pwd)/kleborate-runner.py /usr/local/bin/Kleborate && \
    Kleborate --version

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["Kleborate"]

# Setup the default directory
WORKDIR /root
