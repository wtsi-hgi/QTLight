# Use the latest available Miniforge image
FROM mercury/eqtl:14_07_2025
RUN conda create -n bgenix_env -c bioconda bgenix
ENV PATH="$PATH:/opt/conda/envs/bgenix_env/bin"
RUN which bgenix

RUN git clone https://github.com/jeffreypullin/quasar.git && cd quasar && mkdir -p build && cd build && cmake .. && make
ENV PATH="$PATH:$PWD/quasar/build"
RUN which quasar
CMD ["bash"]
