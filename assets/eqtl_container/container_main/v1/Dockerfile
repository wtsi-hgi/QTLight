# Use the latest available Miniforge image
FROM condaforge/miniforge3:24.7.1-2


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
# RUN apt install libcurl4-gnutls-dev --yes
# install.packages(c("devtools", "roxygen2", "pkgdown", "covr")); devtools::install_dev_deps()
# RUN apt install r-cran-ragg r-cran-xml2
# RUN apt install libtiff-dev  --yes
# RUN apt-get install libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev build-essential libcurl4-openssl-dev libxml2-dev libssl-dev libfontconfig1-dev --yes
# Copy environment.yml and create the Conda environment

COPY environment.yml /tmp/environment.yml
RUN conda config --set always_yes yes --set changeps1 no && \
    conda env create -f /tmp/environment.yml && \
    conda clean -afy && \
    echo "conda activate eqtl-env" >> ~/.bashrc

# Activate the Conda environment
ENV PATH="/opt/conda/envs/eqtl-env/bin:$PATH"
ENV PKG_CONFIG_PATH=/opt/conda/envs/eqtl-env/lib/pkgconfig
SHELL ["conda", "run", "-n", "eqtl-env", "/bin/bash", "-c"]

# Install Python libraries
RUN pip install --no-cache-dir \
    Pgenlib \
    git+https://github.com/saketkc/pysctransform.git

# Install R packages
RUN R --no-echo -e "install.packages(c(  'remotes'), repos='http://cran.us.r-project.org')"
RUN R --no-echo -e "remotes::install_version('matrixStats', repos='http://cran.us.r-project.org')"
RUN apt-get update && apt-get install -y \
    gfortran \
    build-essential \
    libreadline-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    wget
RUN pip install matplotlib==3.5.3
RUN pip install --upgrade scanpy seaborn pandas
RUN R --no-echo -e "devtools::install_version('spam', version = '2.9-0', repos = 'http://cran.us.r-project.org')"
# # Install BCFtools
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 \
  && tar jxvf bcftools-1.19.tar.bz2 \
  && cd bcftools-1.19 \
  && ./configure \
  && make \
  && make install

# Install HTSlib
RUN mkdir -p /opt/htslib && \
    cd /opt/htslib && \
    wget -q https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2 && \
    tar -xjf htslib-1.3.1.tar.bz2 && \
    cd htslib-1.3.1 && make && \
    rm -f /opt/htslib/htslib-1.3.1.tar.bz2

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libgsl-dev \
    libboost-all-dev \
    libx11-dev \
    libxt-dev \
    libreadline-dev \
    libpcre2-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Rmath library
RUN mkdir -p /opt/rmath && \
    cd /opt/rmath && \
    wget https://cran.r-project.org/src/base/R-3/R-3.2.2.tar.gz && \
    tar xzvf R-3.2.2.tar.gz && \
    cd R-3.2.2 && \
    ./configure --without-x --enable-R-shlib && \
    cd src/nmath/standalone && \
    make

RUN mkdir -p /opt/rmath && \
    cd /opt/rmath && \
    wget https://cran.r-project.org/src/base/R-3/R-3.2.2.tar.gz && \
    tar xzvf R-3.2.2.tar.gz && \
    cd R-3.2.2 && \
    ./configure && \
    cd src/nmath/standalone && \
    make

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update && \
    apt-get install -y gcc-11 g++-11


RUN R --no-echo -e "install.packages(c(  'SeuratObject'), repos='http://cran.us.r-project.org')"

RUN apt-get update && apt-get install -y --no-install-recommends \
    libpng-dev \            
    libglpk-dev \   
    libgmp-dev \  
    libhdf5-dev libzstd-dev libcairo2-dev \ 
    python3-dev python3-pip \
    cmake \  
    && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN pip install leidenalg

RUN R --version

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libpng-dev libjpeg-dev libtiff-dev libfreetype6-dev \
    libharfbuzz-dev libfribidi-dev libcairo2-dev \
    libglpk-dev libgmp-dev libzstd-dev libhdf5-dev \
    libgit2-dev libx11-dev libxt-dev libfontconfig1-dev \
    python3-dev libpython3-dev zlib1g-dev build-essential gfortran cmake

RUN R --no-echo -e "install.packages(c(  'curl'), repos='http://cran.us.r-project.org')"
RUN R --no-echo -e 'library(remotes); remotes::install_github("satijalab/seurat", "seurat5", quiet = FALSE)'
RUN R --no-echo -e "remotes::install_github(c('mojaveazure/seurat-disk'), dependencies = TRUE)"

# Install QTLtools
RUN mkdir -p /opt/qtltools && \
    cd /opt/qtltools && \
    wget -q https://qtltools.github.io/qtltools/binaries/QTLtools_1.3.1_source.tar.gz && \
    tar -xzf QTLtools_1.3.1_source.tar.gz && \
    cd qtltools && \
    sed -i 's|BOOST_INC=.*|BOOST_INC=/usr/include|' Makefile && \
    sed -i 's|BOOST_LIB=.*|BOOST_LIB=/usr/lib/x86_64-linux-gnu|' Makefile && \
    sed -i 's|RMATH_INC=.*|RMATH_INC=/opt/rmath/R-3.2.2/src/include|' Makefile && \
    sed -i 's|RMATH_LIB=.*|RMATH_LIB=/opt/rmath/R-3.2.2/src/nmath/standalone|' Makefile && \
    sed -i 's|HTSLD_INC=.*|HTSLD_INC=/opt/htslib/htslib-1.3.1|' Makefile && \
    sed -i 's|HTSLD_LIB=.*|HTSLD_LIB=/opt/htslib/htslib-1.3.1|' Makefile && \
    make && \
    rm -rf /opt/qtltools/QTLtools_1.3.1_source.tar.gz

ENV PATH="/opt/qtltools/qtltools/bin:${PATH}"


# Set library paths
ENV RMATH_INC=/opt/rmath/R-3.2.2/src/include \
    RMATH_LIB=/opt/rmath/R-3.2.2/src/nmath/standalone \
    HTSLD_INC=/opt/htslib/htslib-1.3.1 \
    HTSLD_LIB=/opt/htslib/htslib-1.3.1
RUN QTLtools --help 
# Clean up temporary files
RUN rm -rf /opt/htslib/htslib-1.3.1.tar.bz2 /opt/rmath/R-3.2.2.tar.gz

# Set default command
# RUN apt-get update
# RUN apt-get update && \
#     apt-get install -y software-properties-common && \
#     add-apt-repository ppa:ubuntu-toolchain-r/test && \
#     apt-get update && \
#     apt-get install -y gcc-11 g++-11
RUN strings /usr/lib/x86_64-linux-gnu/libstdc++.so.6 | grep GLIBCXX

RUN apt-get autoremove -y && apt-get clean
ENV LD_LIBRARY_PATH="/opt/libstdc++-old:$LD_LIBRARY_PATH"


RUN strings /usr/lib/x86_64-linux-gnu/libstdc++.so.6 | grep GLIBCXX
RUN strings /usr/lib/x86_64-linux-gnu/libstdc++.so.6 | grep CXXABI

RUN R --no-echo -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='http://cran.us.r-project.org')"
RUN R --no-echo -e "BiocManager::install(c('org.Hs.eg.db', 'AnnotationDbi', 'fgsea', 'biomaRt'))"
RUN R --no-echo -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz', repos = NULL, type = 'source')"

CMD ["bash"]
