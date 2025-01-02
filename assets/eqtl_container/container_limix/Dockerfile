# Use the latest available Miniforge image
FROM condaforge/miniforge3

# Set environment variables
ENV PATH="/opt/conda/bin:${PATH}" \
    DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libgit2-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy the environment.yml file into the container
COPY environment.yml /tmp/environment.yml

# Create the Conda environment from the YAML file
RUN conda config --set always_yes yes --set changeps1 no && conda env create -f /tmp/environment.yml && \
    conda clean -afy

# Install R and necessary R packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends r-base && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    R --slave -e 'install.packages(c("BiocManager"))' && \
    R --slave -e 'BiocManager::install(c("qvalue","multtest","rhdf5","dplyr","readr"))'

# Activate the environment in every shell session
ENV PATH="/opt/conda/envs/limix-env/bin:$PATH"

# Set default Conda environment for Python
RUN echo "conda activate limix-env" >> ~/.bashrc
SHELL ["conda", "run", "-n", "limix-env", "/bin/bash", "-c"]


RUN pip install --no-cache-dir pgenlib
RUN pip install --no-cache-dir h5py
RUN pip install bgen-reader


RUN conda install scikit-learn dask -c conda-forge
RUN conda install -c conda-forge glimix-core

# RUN pip install glimix-core
# Clone the Limix QTL repository
RUN git clone https://github.com/single-cell-genetics/limix_qtl.git
ENV PYTHONPATH="/limix_qtl/Limix_QTL:${PYTHONPATH}"
RUN which python
RUN which pip
RUN conda install debugpy
# Set the default command to /bin/bash
CMD ["/bin/bash"]
