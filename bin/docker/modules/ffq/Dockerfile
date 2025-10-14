## BUILD FFQ
FROM ubuntu:22.04 AS FFQ
ENV DEBIAN_FRONTEND noninteractive
ENV FFQ_VER 0.3.1
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        jq \
        time \
        pigz \
        build-essential \
        make \
        wget \
        python3 \
        python3-pip \
    && wget --no-check-certificate \
        https://github.com/pachterlab/ffq/archive/refs/tags/v${FFQ_VER}.tar.gz \
    && tar -xzf v${FFQ_VER}.tar.gz \
    && python3 -m pip install -e ffq-${FFQ_VER} \
    && apt-get -y clean && apt-get -y autoclean \
    && rm -rf /var/lib/apt/lists/*
ENV HOME="/home"
RUN chmod -R 777 /home
ENV TIME="\n\nres memory (KB)     %M\ntime (HH:mm:ss)     %E"
ENTRYPOINT ["/usr/bin/time"]