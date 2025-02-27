FROM wzhou88/saige:1.3.3

RUN R --no-echo -e "install.packages('BiocManager',repos = 'http://cran.us.r-project.org')"
RUN R --no-echo -e "BiocManager::install('qvalue')"