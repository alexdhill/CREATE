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

## BUID MULTIQC
FROM ubuntu:22.04 as MULTIQC
ENV MULTIQC_VER 1.8
RUN apt-get update \
    && apt-get install -y --no-install-recommends wget \
    && wget --no-check-certificate \
        https://github.com/jgm/pandoc/releases/download/2.7.3/pandoc-2.7.3-linux.tar.gz \
    && tar -xvf pandoc-2.7.3-linux.tar.gz \
    && apt-get -y clean && apt-get -y autoclean \
    && rm -rf /var/lib/apt/lists/*