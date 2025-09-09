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

## BUILD R
FROM ubuntu:22.04 as R
ENV R_VER 4.5
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        wget \
        time \
        pigz \
        ca-certificates \
        r-base-dev=${R_VER}-1ubuntu2 \
    Rscript -e \
       "install.packages(c( \
            'devtools', \
            'tidyr', \
            'dplyr', \
            'readr', \
            'stringr', \
            'remotes', \
            'qqman', \
            'ggplot2', \
            'ggforce', \
            'ggupset', \
            'ggrepel', \
            'ggdist', \
            'ggraph', \
            'ggh4x', \
            'ggbump', \
            'ggfortify', \
            'ggdag', \
            'ggridges', \
            'gggenes', \
            'ggpointdensity', \
            'ggdensity', \
            'ggcorrplot', \
            'ggheatmap', \
            'ggmuller', \
            'ggpubr', \
            'ggsignif', \
            'ggsci', \
            'gganimate', \
            'ggthemes', \
            'patchwork', \
            'RColorBrewer', \
            'argparse', \
            'import', \
            'Matrix', \
            'reticulate', \
            'languageserver' \
       )); \
       remotes::install_version('dbplyr', version='2.3.4'); \
       remotes::install_version('matrixStats', version='1.1.0'); \
       install.packages( \
           'BiocManager', \
           repos='http://cran.us.r-project.org' \
       ); \
       BiocManager::install(c( \
           'tidyverse', \
           'biocLite', \
           'BiocIO', \
           'Biostrings', \
           'GenomicFeatures', \
           'SummarizedExperiment', \
           'SingleCellExperiment', \
           'HDF5Array', \
           'DESeq2', \
           'tximeta', \
           'tximport', \
           'eisaR', \
           'pcaMethods', \
           'multtest', \
           'metap', \
           'scp', \
           'scater', \
           'scuttle', \
           'scran', \
           'scRNAseq', \
           'fishpond', \
           'velociraptor', \
           'BSgenome', \
           'DECIPHER', \
           'DRIMSeq', \
            'stageR' \
       )); \
       remotes::install_github('satijalab/seurat'); \
       remotes::install_github('satijalab/seurat-wrappers'); \
       remotes::install_github('mojaveazure/seurat-disk'); \
       remotes::install_github('stuart-lab/signac'); \
       remotes::install_github('velocyto-team/velocyto.R');" \
    && apt-get -y clean && apt-get -y autoclean \
    && rm -rf /var/lib/apt/lists/*
ENV HOME="/home"
RUN chmod -R 777 /home
ENV TIME="\n\nres memory (KB)     %M\ntime (HH:mm:ss)     %E"
ENTRYPOINT ["/usr/bin/time"]