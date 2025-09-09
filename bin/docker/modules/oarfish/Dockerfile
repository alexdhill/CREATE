## BUILD OARFISH
FROM ubuntu:22.04 as OARFISH
ENV DEBIAN_FRONTEND noninteractive
ENV OARFISH_VER 0.9.0
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        time \
        pigz \
        build-essential \
        make \
        cmake \
        curl \
        ca-certificates \
        unzip \
    && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y \
    && curl -kL \
        https://github.com/COMBINE-lab/oarfish/archive/refs/tags/v${OARFISH_VER}.zip \
        -o v${OARFISH_VER}.zip \
    && unzip v${OARFISH_VER}.zip \
    && cd oarfish-${OARFISH_VER} \
    && $HOME/.cargo/bin/cargo build --release \
    && apt-get -y clean && apt-get -y autoclean \
    && rm -rf /var/lib/apt/lists/*
ENV HOME="/home"
RUN chmod -R 777 /home
ENV TIME="\n\nres memory (KB)     %M\ntime (HH:mm:ss)     %E"
ENTRYPOINT ["/usr/bin/time"]