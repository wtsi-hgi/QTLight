FROM wzhou88/saigeqtl:0.1.0

RUN R --no-echo -e "install.packages('BiocManager',repos = 'http://cran.us.r-project.org')"
RUN R --no-echo -e "BiocManager::install('qvalue')"