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

## BUILD PYTHON
FROM ubuntu:22.04 as PYTHON
ENV PYTHON_VER 3.10
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        wget \
        time \
        pigz \
        ca-certificates \
        python${PYTHON_VER}-dev \
        python3-pip \
    && python3 -m pip install --no-cache-dir --upgrade pip \
    && python3 -m pip install --no-cache-dir patchelf \
    && python3 -m pip install --no-binary :all --no-cache-dir \
        "multiqc==$MULTIQC_VER" \
    && python3 -m pip install \
        tqdm \
        ncls \
        rpy2 \
        mappy \
        numpy \
        pysam \
        pandas \
        seaborn \
        cutadapt \
        biopython \
        matplotlib \
        pybedtools \
        fasta_digest \
    && python3 -m pip install -e /ffq \
    && apt-get -y clean && apt-get -y autoclean \
    && rm -rf /var/lib/apt/lists/*
ENV HOME="/home"
RUN chmod -R 777 /home
ENV TIME="\n\nres memory (KB)     %M\ntime (HH:mm:ss)     %E"
ENTRYPOINT ["/usr/bin/time"]