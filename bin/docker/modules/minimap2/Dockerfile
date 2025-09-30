## REQUIRED NOTICE: Copyright (c) 2020-2023, Regents of the University of California
## All rights reserved. https://polyformproject.org/licenses/noncommercial/1.0.0
## 
## This software was developed by the Daniel Kim lab at the University of California, Santa Cruz.
## Authors: Roman E. Reggiardo, Vikas Peddu, Alex D. Hill
## 
## The licensor grants you a copyright license for the software to do everything you might do with
## the software that would otherwise infringe the licensor’s copyright in it for any permitted
## purpose.
## 
## As far as the law allows, the software comes as is, without any warranty or condition, and the
## licensor will not be liable to you for any damages arising out of these terms or the use or
## nature of the software, under any kind of legal claim.

## BUILD SAMTOOLS
FROM ubuntu:22.04 as SAMTOOLS
ENV DEBIAN_FRONTEND noninteractive
ENV SAMTOOLS_VER 1.19
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        time \
        pigz \
        build-essential \
        make \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        autoconf \
        wget \
    && wget --no-check-certificate \
        https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 \
    && tar -xvf samtools-1.19.tar.bz2 \
    && cd samtools-1.19 \ 
    && ./configure --prefix=/usr/local/bin/ \
    && make \
    && make install \
    && cp samtools /usr/local/bin/

## BUILD MINIMAP2
FROM ubuntu:22.04 as MINIMAP2
ENV DEBIAN_FRONTEND noninteractive
ENV MINIMAP2_VER 2.26
COPY --from=SAMTOOLS /usr/local/bin/samtools /usr/local/bin/
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        time \
        pigz \
        build-essential \
        make \
        zlib1g-dev \
        wget \
    && wget --no-check-certificate \
        https://github.com/lh3/minimap2/archive/refs/tags/v${MINIMAP2_VER}.tar.gz \
    && tar -xvf v2.26.tar.gz \
    && cd minimap2* \
    && make \
    && cp minimap2 /usr/local/bin/ \
    && apt-get -y clean && apt-get -y autoclean \
    && rm -rf /var/lib/apt/lists/*
ENV HOME="/home"
RUN chmod -R 777 /home
ENV TIME="\n\nres memory (KB)     %M\ntime (HH:mm:ss)     %E"
ENTRYPOINT ["/usr/bin/time"]