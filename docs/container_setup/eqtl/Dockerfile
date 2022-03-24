FROM ubuntu
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

WORKDIR ${HOME}/tmp2

COPY limix_fork ${HOME}
RUN ls -lh ${HOME}
RUN chmod +x ${HOME}/limix_qtl/Limix_QTL/run_QTL_analysis.py
RUN chmod +x ${HOME}/limix_qtl/Limix_QTL/post_processing/minimal_postprocess.py
RUN chmod +x ${HOME}/limix_qtl/Limix_QTL/run_QTL_analysis.py
RUN ln -s ${HOME}/limix_qtl/Limix_QTL/run_QTL_analysis.py /usr/bin/limix_run
RUN ln -s ${HOME}/limix_qtl/Limix_QTL/post_processing/minimal_postprocess.py /usr/bin/limix_aggregate

RUN which limix_run
COPY environment.yml ${HOME}/tmp2/environment.yml
RUN conda config --set always_yes yes --set changeps1 no && \
    conda env update -v -n root --file ${HOME}/tmp2/environment.yml && \
    conda list --name root && \
    conda clean -atipy && \
    rm ${HOME}/tmp2/environment.yml
RUN which conda    
RUN which conda 
RUN which limix
RUN which R && R --version
RUN python -V

USER ${NB_USER}
