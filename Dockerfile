FROM ubuntu-22.04

LABEL stage=builder

LABEL maintainer="tjakobi@arizona.edu"

ARG MAKEFLAGS="-j32"

# Set non-interactive frontend to avoid prompts during installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Phoenix

# Update packages and install necessary dependencies
RUN apt update && \
    apt-get install -y \
    git \
    wget \
    unzip \
    python3-pip \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libncurses-dev \
    libcurl4-openssl-dev

RUN mkdir /build/

# Download and install BEDTools
RUN cd /build/ && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static && \
    chmod +x bedtools.static && \
    mv bedtools.static /usr/local/bin/bedtools

# Download and install pblat
RUN cd /build/ && \
    git clone https://github.com/icebert/pblat.git && \
    cd pblat && \
    make && \
    cp pblat /usr/local/bin/

# Download and install samtools
RUN cd /build/ && \
    wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar xvf samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && \
    make && \
    make install

RUN rm /build/ -rf

# Download and install Python libraries
RUN pip install pyyaml argparse nanofilt pybedtools tqdm requests --break-system-packages -v

ADD . /app/

RUN ln -s /app/long_read_circRNA /usr/local/bin/

# Check installation of tools
RUN NanoFilt --version
RUN bedtools --version
RUN samtools --version
RUN pblat
RUN long_read_circRNA check

RUN mkdir /data/

ENTRYPOINT ["long_read_circRNA"]
