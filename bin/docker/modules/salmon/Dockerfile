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

## BUILD SALMON
FROM ubuntu:22.04 as SALMON
ENV DEBIAN_FRONTEND noninteractive
ENV SALMON_VER 1.10.1
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        make \
        wget \
        curl \
        unzip \
        cmake \
        libgff-dev \
        libtbb-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        build-essential \
        ca-certificates \
        libboost-all-dev \
        libcurl4-openssl-dev \
    && wget --no-check-certificate \
        https://github.com/COMBINE-lab/salmon/archive/refs/tags/v${SALMON_VER}.tar.gz \
    && tar -xvf v${SALMON_VER}.tar.gz \
    && cd salmon-${SALMON_VER} \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install
RUN ls -la /salmon*/