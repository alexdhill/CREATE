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
    && wget --no-check-certificate -q \
        https://github.com/arq5x/bedtools2/archive/refs/tags/v${BEDTOOLS_VER}.tar.gz \
    && tar -xzf v${BEDTOOLS_VER}.tar.gz \
    && cd bedtools2-${BEDTOOLS_VER} \
    && make \
    && make install

## BUILD SAMTOOLS
FROM ubuntu:22.04 as SAMTOOLS
ENV DEBIAN_FRONTEND noninteractive
ENV SAMTOOLS_VER 1.19
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
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
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        build-essential \
        make \
        zlib1g-dev \
        wget \
    && wget --no-check-certificate \
        https://github.com/lh3/minimap2/archive/refs/tags/v${MINIMAP2_VER}.tar.gz \
    && tar -xvf v2.26.tar.gz \
    && cd minimap2* \
    && make

## BUILD FLAIR
FROM ubuntu:22.04 as FLAIR
ENV DEBIAN_FRONTEND noninteractive
ENV FLAIR_VER d777b30
COPY --from=BEDTOOLS /usr/local/bin/bedtools /usr/local/bin/bedtools
COPY --from=SAMTOOLS /usr/local/bin/samtools /usr/local/bin/samtools
COPY --from=MINIMAP2 /minimap2*/minimap2 /usr/local/bin/minimap2
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        git \
        wget \
        unzip \
        gpg-agent \
        r-base-dev \
        build-essential \
        python3-dev \
        python3-pip \
        python3-setuptools \
        software-properties-common \
        libmkl-dev \
        libblis-dev \
        liblapack-dev \
        libopenblas-dev \
        libcurl4-openssl-dev \
    && add-apt-repository universe \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        libatlas-base-dev \
    && python3 -m pip install --upgrade pip \
    && python3 -m pip install \
        tqdm \
        ncls \
        rpy2 \
        mappy \
        numpy \
        pysam \
        pandas \
        seaborn \
        pipettor \
        matplotlib \
        pybedtools \
    && Rscript -e " \
        install.packages(c( \
            'ggplot2', \
            'dplyr', \
            'tidyr', \
            'qqman', \
            'BiocManager' \
        ))\
    "\
    && Rscript -e "\
        BiocManager::install(c( \
            'DESeq2', \
            'DRIMSeq', \
            'stageR' \
        )) \
    "
RUN git clone https://www.github.com/BrooksLabUCSC/FLAIR \
    && cd FLAIR \
    && git checkout ${FLAIR_VER}
ENV PATH=$PATH:/FLAIR/bin