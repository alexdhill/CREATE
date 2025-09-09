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

## BUILD BEDTOOLS
FROM ubuntu:22.04 as BEDTOOLS
ENV DEBIAN_FRONTEND noninteractive
ENV BEDTOOLS_VER 2.31.0
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        build-essential \
        zlib1g-dev \
        libghc-bzlib-dev \
        liblzma-dev \
        wget \
        python3 \
        time \
        pigz \
    && wget --no-check-certificate -q \
        https://github.com/arq5x/bedtools2/archive/refs/tags/v${BEDTOOLS_VER}.tar.gz \
    && tar -xzf v${BEDTOOLS_VER}.tar.gz \
    && cd bedtools2-${BEDTOOLS_VER} \
    && make \
    && make install \
    && apt-get -y clean && apt-get -y autoclean \
    && rm -rf /var/lib/apt/lists/*
ENV HOME="/home"
RUN chmod -R 777 /home
ENV TIME="\n\nres memory (KB)     %M\ntime (HH:mm:ss)     %E"
ENTRYPOINT ["/usr/bin/time"]