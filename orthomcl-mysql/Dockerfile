FROM ivasilyev/conda-base:latest
# Base image CLI:
# IMG=ivasilyev/conda-base:latest && docker pull ${IMG} && docker run --rm --net=host -it ${IMG} bash

RUN echo "Software install" && \
    apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y  \
            cpanminus \
            mysql-server \
            libmysqlclient-dev

RUN echo "Cleanup" && \
    apt-get autoremove --purge -y && \
    apt-get clean all && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda config --add channels bioconda && \
    conda upgrade -y conda

RUN echo "Install software via Conda" && \
    conda install -c bioconda -y diamond blast mcl orthomcl && \
    conda update --all -y

RUN echo "Install modules via CPANM" && \
    cpanm DBI DBD::mysql

# update-alternatives: using /etc/mysql/mysql.cnf to provide /etc/mysql/my.cnf (my.cnf) in auto mode
COPY orthomcl-mysql/my.cnf /etc/mysql/

# Add scripts to the image
RUN mkdir -p /opt/my_tools
COPY orthomcl-mysql/run_orthomcl.py orthomcl-mysql/pipeline_handler.py orthomcl-mysql/orthomcl.config /opt/my_tools/
RUN chmod a+x /opt/my_tools/*
ENV PATH=${PATH}:/opt/my_tools

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["pipeline_handler.py"]