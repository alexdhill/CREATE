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

## BUILD FASTP
FROM ubuntu:22.04 as FASTP
ENV FASTP_VER 0.23.4
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential automake autogen autoconf libtool \
        wget git ca-certificates \
        nasm yasm cmake \
    && git clone https://github.com/intel/isa-l.git \
    && cd isa-l \
    && ./autogen.sh \
    && ./configure \
    && make -j8 \
    && make -j8 install \
    && cd .. \
    && git clone https://github.com/ebiggers/libdeflate.git \
    && cd libdeflate \
    && cmake -B build \
    && cmake --build build \
    && cmake --install build \
    && cd .. \
    && wget https://github.com/OpenGene/fastp/archive/refs/tags/v${FASTP_VER}.tar.gz \
    && tar -xvf v${FASTP_VER}.tar.gz \
    && cd fastp* \
    && make -j8 \
    && make -j8 install \
    && cd ..