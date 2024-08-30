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

## BUILD STAR
FROM ubuntu:22.04 as STAR
ENV DEBIAN_FRONTEND noninteractive
ENV STAR_VER 2.7.11a
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        build-essential \
        make \
        wget \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        autoconf \
        xxd \
    && wget --no-check-certificate \
        https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VER}.tar.gz \
    && tar -xvf ${STAR_VER}.tar.gz \
    && cd STAR-${STAR_VER}/source \
    && make STAR
RUN ls -lahR /STAR*/bin/