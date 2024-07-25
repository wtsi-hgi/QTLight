FROM ubuntu:20.04
RUN cat /etc/os-release

# RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install --yes -y --no-install-recommends dirmngr apt-transport-https ca-certificates software-properties-common gnupg2 \
    build-essential curl git \ 
    g++ gcc gfortran make autoconf automake libtool \
    zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
    libblas-dev libx11-dev \
    libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*
ENV NB_USER container_user
# RUN apt-get install wget
ENV HOME /home/${NB_USER}
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
# RUN which R && R --version
RUN apt-get install wget
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh -O /tmp/miniconda.sh && chmod +x  /tmp/miniconda.sh && bash /tmp/miniconda.sh -b -p ${HOME}/conda
# RUN bash Miniconda3-latest-Linux-x86_64.sh
ENV PATH ${HOME}/conda/bin:$PATH



# set up the user directory structure
RUN mkdir ${HOME}/tmp2 && \
    mkdir ${HOME}/analysis
    
WORKDIR ${HOME}

RUN git clone https://github.com/single-cell-genetics/limix_qtl.git
RUN ls -lh ${HOME}
RUN chmod +x ${HOME}/limix_qtl/Limix_QTL/run_QTL_analysis.py
RUN chmod +x ${HOME}/limix_qtl/Limix_QTL/post_processing/minimal_postprocess.py
RUN chmod +x ${HOME}/limix_qtl/Limix_QTL/run_QTL_analysis.py
RUN ln -s ${HOME}/limix_qtl/Limix_QTL/run_QTL_analysis.py /usr/bin/limix_run
RUN ln -s ${HOME}/limix_qtl/Limix_QTL/post_processing/minimal_postprocess.py /usr/bin/limix_aggregate

RUN which limix_run

RUN conda config --set always_yes yes --set changeps1 no && \
    conda list --name root && \
    conda clean -atipy 
RUN which conda    
RUN pip install glimix-core
RUN git clone https://github.com/broadinstitute/tensorqtl.git
RUN cd tensorqtl && pip install -r install/requirements.txt .
RUN pip install Pgenlib
RUN apt install zip --yes
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && mv bedtools.static.binary /usr/bin/bedtools && chmod a+x /usr/bin/bedtools
RUN cd /usr/bin/ && wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip && unzip plink_linux_x86_64_20231211.zip
RUN rm -r tensorqtl
RUN ls
RUN rm -r tmp2 analysis
RUN rm -r /tmp/miniconda.sh
RUN pip install gtfparse==1.2.1
RUN pip install pandas==1.1.3
RUN which python
RUN conda install -c conda-forge r-base=4.4.0   
RUN which R
RUN pip install h5py numpy pip install pandas-plink bgen-reader bgen matplotlib-venn
RUN pip install anndata seaborn
RUN pip install -U scikit-learn
RUN pip install statsmodels==0.12.2
RUN pip install numba==0.53.1 scanpy debugpy tensorqtl numpy_sugar
RUN pip install tables igraph leidenalg glimix-core==3.1.12
RUN pip install rpy2
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20240704.zip && unzip plink2_linux_x86_64_20240704.zip
RUN mv ./plink2 /usr/bin/
RUN R --no-echo -e "install.packages('ggfortify',repos = 'http://cran.us.r-project.org')"
RUN R --no-echo -e "install.packages('dplyr',repos = 'http://cran.us.r-project.org')"
RUN R --no-echo -e "install.packages('readr',repos = 'http://cran.us.r-project.org')"
RUN R --no-echo -e "install.packages('RColorBrewer',repos = 'http://cran.us.r-project.org')"
RUN R --no-echo -e "install.packages('devtools',repos = 'http://cran.us.r-project.org')"
RUN R --no-echo -e "install.packages('rlang',repos = 'http://cran.us.r-project.org')"
RUN pip install pandas==1.5 
RUN R --no-echo -e "install.packages('BiocManager',repos = 'http://cran.us.r-project.org')"
RUN R --no-echo -e "BiocManager::install('edgeR')"
RUN R --no-echo -e "BiocManager::install('pcatools')"
RUN R --no-echo -e "BiocManager::install('multtest')"
RUN R --no-echo -e "BiocManager::install('rhdf5')"
RUN R --no-echo -e "BiocManager::install('DESeq2')"
RUN R --no-echo -e "BiocManager::install('PCAtools')"
RUN R --no-echo -e 'install.packages("Matrix", repos = "https://cran.r-project.org")'
RUN R --no-echo -e "remotes::install_version('matrixStats', repos = 'http://cran.us.r-project.org')"
RUN R --no-echo -e "BiocManager::install('qvalue')"
RUN pip install --upgrade scanpy
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 \
  && tar jxvf bcftools-1.19.tar.bz2 \
  && cd bcftools-1.19 \
  && ./configure \
  && make \
  && make install