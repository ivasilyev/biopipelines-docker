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

RUN echo "Install Perl modules via CPAN" && \
    cpanm DBI DBD::mysql Bio::Perl

RUN echo "Install Python 3 packages" && \
    pip3 install PyMySQL numpy pandas biopython

# update-alternatives: using /etc/mysql/mysql.cnf to provide /etc/mysql/my.cnf (my.cnf) in auto mode
COPY orthomcl-mysql/my.cnf /etc/mysql/
RUN mkdir -p /nonexistent && \
    chown root:root /etc/mysql/my.cnf && \
    chmod 600 /etc/mysql/my.cnf

# Add scripts to the image
RUN mkdir -p /opt/my_tools
COPY orthomcl-mysql/*.py \
     orthomcl-mysql/orthomcl.config \
     /opt/my_tools/
RUN chmod a+x /opt/my_tools/*
ENV PATH=${PATH}:/opt/my_tools

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD service mysql restart; \
    python3 /opt/my_tools/pipeline_handler.py

# To use MySQL, start it explicitly, e.g:
#export IMG=ivasilyev/orthomcl-mysql:latest && \
#docker pull ${IMG} && \
#docker run --rm --net=host -it ${IMG} \
#bash -c \
#    'service mysql restart
#     echo foo
#     echo bar
#     echo baz'
