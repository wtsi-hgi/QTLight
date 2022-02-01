FROM conda/miniconda3
RUN cat /etc/os-release
RUN apt-get update && \
    apt-get install --yes --no-install-recommends dirmngr apt-transport-https ca-certificates software-properties-common gnupg2 \
    build-essential curl git python-pip procps \ 
    g++ gcc gfortran make autoconf automake libtool \
    zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
    libblas-dev libx11-dev \
    libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
# RUN which R && R --version

ENV HOME /home

# set up the user directory structure
RUN mkdir ${HOME}/tmp && \
    mkdir ${HOME}/analysis



WORKDIR ${HOME}/tmp
COPY environment.yml ${HOME}/tmp/environment.yml
COPY limix_fork ${HOME}
RUN ls -lh ${HOME}
RUN chmod +x /home/limix_qtl/Limix_QTL/run_QTL_analysis.py
RUN chmod +x /home/limix_qtl/Limix_QTL/post_processing/minimal_postprocess.py
RUN chmod +x /home/limix_qtl/Limix_QTL/run_QTL_analysis.py
RUN ln -s /home/limix_qtl/Limix_QTL/run_QTL_analysis.py /usr/bin/limix_run
RUN ln -s /home/limix_qtl/Limix_QTL/post_processing/minimal_postprocess.py /usr/bin/limix_aggregate

RUN which limix_run
RUN conda config --set always_yes yes --set changeps1 no && \
    conda update conda --yes && \
    conda env update -v -n root --file ${HOME}/tmp/environment.yml && \
    conda list --name root && \
    conda clean -atipy && \
    rm ${HOME}/tmp/environment.yml
RUN conda update readline --yes 
RUN which limix
RUN which R && R --version


