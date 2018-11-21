# docker run -it --rm --net=host -p 8888:8888 qiime2/core:latest bash

ENV DEBIAN_FRONTEND noninteractive

RUN echo "Software update" && \
    apt-get -y update && \
    apt-get -y install software-properties-common apt-utils dpkg gnupg2 curl

RUN curl -sL https://deb.nodesource.com/setup_11.x | bash - && \
    apt-get install -y nodejs npm

RUN apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN Q2_VER=2018.11.0 && \
    curl -fsSL https://github.com/qiime2/q2studio/archive/${Q2_VER}.tar.gz | tar -xvzf - && \
    cd q2studio-${Q2_VER}

RUN pip install -e . && \
    npm install

# npm run dev --host 0.0.0.0 --port 8888
