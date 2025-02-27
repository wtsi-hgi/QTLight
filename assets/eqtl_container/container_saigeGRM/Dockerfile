FROM condaforge/miniforge3:latest

# Set environment variables
ENV PATH="/opt/conda/bin:${PATH}" \
    DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    NB_USER=container_user

# Update and install base system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    unzip \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libgsl-dev \
    libboost-all-dev \
    dirmngr \
    apt-transport-https \
    ca-certificates \
    software-properties-common \
    gnupg2 \
    curl \
    git \
    g++ \
    gcc \
    gfortran \
    make \
    autoconf \
    automake \
    libtool \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    lbzip2 \
    libblas-dev \
    libx11-dev \
    libreadline-dev \
    libxt-dev \
    libpcre2-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN apt-get update && \
    apt-get install --yes -y --no-install-recommends dirmngr apt-transport-https ca-certificates software-properties-common gnupg2 \
    build-essential curl git \ 
    g++ gcc gfortran make autoconf automake libtool \
    zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
    libblas-dev libx11-dev \
    libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*

RUN conda create -y -n saige -c conda-forge -c bioconda "r-base>=4.0" r-saige

RUN echo "conda activate saige" >> ~/.bashrc

# Activate the Conda environment
ENV PATH="/opt/conda/envs/saige/bin:$PATH"
SHELL ["conda", "run", "-n", "saige", "/bin/bash", "-c"]

RUN createSparseGRM.R  --help
# RUN which conda