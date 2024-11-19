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


## BUILD FFQ
FROM ubuntu:22.04 AS FFQ
ENV DEBIAN_FRONTEND noninteractive
ENV FFQ_VER 0.3.1
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        build-essential \
        make \
        wget \
        python3 \
        python3-pip \
    && wget --no-check-certificate \
        https://github.com/pachterlab/ffq/archive/refs/tags/v${FFQ_VER}.tar.gz \
    && tar -xzf v${FFQ_VER}.tar.gz \
    && python3 -m pip install -e ffq-${FFQ_VER}


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


## BUID GFFREAD
FROM ubuntu:22.04 as GFFREAD
ENV DEBIAN_FRONTEND noninteractive
ENV GFFREAD_VER 0.12.1
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        build-essential \
        make \
        wget \
        git \
    && git clone https://github.com/gpertea/gffread --recursive \
    && cd ./gffread \
    && make release \
    && ls -l


## BUILD FASTQC
FROM ubuntu:22.04 as FASTQC
ENV FASTQC_VER 0.12.1
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        unzip \
        wget \
        perl \
        default-jre \
        procps \
    && wget --no-check-certificate \
        https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VER}.zip \
    && unzip fastqc_v${FASTQC_VER}.zip \
    && rm fastqc_v${FASTQC_VER}.zip \
    && chmod +x FastQC/fastqc


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


## BUILD TRIMGALORE
FROM ubuntu:22.04 as TRIMGALORE
ENV TRIMGALORE_VER 0.6.10
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
    && curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz \
        -o trim_galore.tar.gz \
    && tar -xvzf trim_galore.tar.gz \
    && mv /TrimGalore-0.6.10/ /TrimGalore


## BUILD PORECHOP
FROM ubuntu:22.04 as PORECHOP
ENV PORECHOP_VER 0.2.4
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        wget \
        python3-pip \
        python3-distutils \
        python3 \
        make \
        g++ \
    && apt-get autoclean && rm -rf /var/lib/apt/lists/* \
    && pip3 install setuptools \
    && wget -q --no-check-certificate\
        https://github.com/rrwick/Porechop/archive/refs/tags/v${PORECHOP_VER}.tar.gz \
    && tar -xf v${PORECHOP_VER}.tar.gz \
    && rm v${PORECHOP_VER}.tar.gz \
    && cd /Porechop-${PORECHOP_VER}


## BUILD CHOPPER
FROM ubuntu:22.04 as CHOPPER
ENV CHOPPER_VER 0.9.0
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
        https://github.com/wdecoster/chopper/archive/refs/tags/v${CHOPPER_VER}.zip \
        -o v${CHOPPER_VER}.zip \
    && unzip v${CHOPPER_VER}.zip \
    && cd chopper-${CHOPPER_VER} \
    && $HOME/.cargo/bin/cargo build --release


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


## BUILD ALEVIN
FROM ubuntu:22.04 as ALEVIN
ENV DEBIAN_FRONTEND noninteractive
ENV ALEVIN_VER 0.10.0
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        build-essential \
        make \
        cmake \
        curl \
        ca-certificates \
        unzip \
    && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y \
    && curl -kL \
        https://github.com/COMBINE-lab/alevin-fry/archive/refs/tags/v${ALEVIN_VER}.zip \
        -o v${ALEVIN_VER}.zip \
    && unzip v${ALEVIN_VER}.zip \
    && cd alevin-fry-* \
    && $HOME/.cargo/bin/cargo build --release


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


## BUILD SEQTK
FROM ubuntu:22.04 as SEQTK
ENV DEBIAN_FRONTEND noninteractive
ENV SEQTK_VER 1.4
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        build-essential \
        make \
        wget \
        unzip \
        libz-dev \
    && wget --no-check-certificate \
        https://github.com/lh3/seqtk/archive/refs/tags/v${SEQTK_VER}.zip \
    && unzip v${SEQTK_VER}.zip \
    && cd seqtk-${SEQTK_VER} \
    && make \
    && cp seqtk /usr/local/bin


## BUILD FLAIR
FROM ubuntu:22.04 as FLAIR
ENV DEBIAN_FRONTEND noninteractive
ENV FLAIR_VER 2.0.0
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        wget \
        unzip \
    && wget --no-check-certificate \
        https://github.com/BrooksLabUCSC/flair/archive/refs/tags/${FLAIR_VER}.zip \
    && unzip ${FLAIR_VER}.zip \
    && chmod +x flair-${FLAIR_VER}/flair.py


## Build Final
FROM ubuntu:22.04 as app
COPY --from=FFQ /ffq-0.3.1 /ffq
COPY --from=SAMTOOLS /usr/local/bin/samtools /usr/local/bin/samtools
COPY --from=BEDTOOLS /usr/local/bin/bedtools /usr/local/bin/bedtools
COPY --from=GFFREAD /gffread/gffread /usr/local/bin/gffread
COPY --from=FASTQC /FastQC /usr/local/bin/FastQC
COPY --from=MULTIQC /pandoc-*-linux/bin/* /usr/local/bin
COPY --from=FASTP /usr/local/bin/fastp /usr/local/bin/fastp
COPY --from=FASTP /usr/local/lib/libdeflate.* /usr/lib/
COPY --from=FASTP /usr/lib/libisal.* /usr/lib/
COPY --from=TRIMGALORE /TrimGalore /usr/local/bin/TrimGalore
COPY --from=PORECHOP /Porechop-* /Porechop
COPY --from=CHOPPER /chopper-*/target/release/chopper /usr/local/bin/chopper
COPY --from=MINIMAP2 /minimap2*/minimap2 /usr/local/bin/minimap2
COPY --from=SALMON /salmon*/bin/salmon /usr/local/bin/salmon
COPY --from=ALEVIN /alevin-fry-*/target/release/alevin-fry /usr/local/bin/alevin-fry
COPY --from=STAR /STAR*/bin/*static/* /usr/local/bin/
COPY --from=SEQTK /usr/local/bin/seqtk /usr/local/bin/seqtk
COPY --from=FLAIR /flair-* /usr/local/bin/flair
ENV DEBIAN_FRONTEND noninteractive
ENV FFQ_VER 0.3.1
ENV SAMTOOLS_VER 1.19
ENV BEDTOOLS_VER 2.31.0
ENV FASTQC_VER 0.12.1
ENV MULTIQC_VER 1.8
ENV FASTP_VER 0.20.1
ENV TRIMGALORE_VER 0.6.10
ENV PORECHOP_VER 0.2.4
ENV CHOPPER_VER 0.9.0
ENV MINIMAP2_VER 2.26
ENV SALMON_VER 1.10.1
ENV ALEVIN_VER 0.8.2
ENV STAR_VER 2.7.11a
ENV SEQTK_VER 1.4
ENV LC_ALL="C.UTF-8" LANG="C.UTF-8"
ENV PATH="$PATH:/usr/local/bin/FastQC/:/usr/local/bin/TrimGalore/:/usr/local/bin/flair/bin/:/usr/local/bin/flair/src/flair/"
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        wget \
        gpg-agent \
        software-properties-common \
    && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
        | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository ppa:openjdk-r/ppa \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
    && add-apt-repository universe \
    && apt-get update -y \
    && apt-get install -y --no-install-recommends \
        time \
        perl \
        pigz \
        locales \
        fftw3-dev \
        r-base-dev \
        gfortran-10 \
        python3-dev \
        python3-pip \
        mysql-server \
        texlive-xetex \
        openjdk-11-jre \
        build-essential \
        texlive-science \
        python3-setuptools \
        libz-dev \
        libmkl-dev \
        libpng-dev \
        libbz2-dev \
        libtbb-dev \
        libssl-dev \
        liblzma-dev \
        libxml2-dev \
        libblas-dev \
        libjpeg-dev \
        libblis-dev \
        libtiff5-dev \
        liblapack-dev \
        libfribidi-dev \
        libharfbuzz-dev \
        libfreetype6-dev \
        libboost-all-dev \
        libatlas-base-dev \
        libpoppler-cpp-dev \
        libgfortran-10-dev \
        libfontconfig1-dev \
        libcurl4-gnutls-dev \
    && cp /usr/lib/x86_64-linux-gnu/libgfortran* /usr/lib/ \
    && cp /usr/lib/gcc/x86_64-linux-gnu/10/libgfortran* /usr/lib/ \
    && cp /usr/bin/gfortran-10 /usr/bin/gfortran \
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
    && Rscript -e \
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
    && cd /Porechop \
    && python3 setup.py install \
    && cd / \
    && rm -rf /Porechop \
    && apt-get -y clean && apt-get -y autoclean \
    && rm -rf /var/lib/apt/lists/*
ENV HOME="/home"
RUN chmod -R 777 /home
ENV TIME="\n\nres memory (KB)     %M\ntime (HH:mm:ss)     %E"
ENTRYPOINT ["/usr/bin/time"]