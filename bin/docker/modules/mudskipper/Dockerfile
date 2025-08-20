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

## BUILD MUDSKIPPER
FROM ubuntu:22.04 as MUDSKIPPER
ENV MUDSKIPPER_VER effd3fa
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        build-essential \
        make \
        curl \
        unzip \
        cmake \
        pkg-config \
        zlib1g-dev \
        ca-certificates \
    && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y \
    && curl -kL \
        https://github.com/OceanGenomics/mudskipper/archive/refs/tags/v${MUDSKIPPER_VER}.zip \
        -o v${MUDSKIPPER_VER}.zip \
    && unzip v${MUDSKIPPER_VER}.zip \
    && cd mudskipper-${MUDSKIPPER_VER} \
    && $HOME/.cargo/bin/cargo update \
    && $HOME/.cargo/bin/cargo build --release