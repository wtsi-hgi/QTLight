FROM nvcr.io/nvidia/jax:23.08-py3

# Avoid rich UI and reduce thread usage in pip/git
ENV PIP_PROGRESS_BAR=off
ENV GIT_HTTP_LOW_SPEED_LIMIT=0
ENV GIT_HTTP_LOW_SPEED_TIME=999999
# ENV GIT_CONFIG_PARAMETERS='protocol.version=2'
COPY jaxqtl /opt/jaxqtl
# Upgrade pip and install with low thread pressure
RUN pip install --upgrade pip && \
        pip install \
        numpy==1.25.2 \
        scipy==1.11.2 \
        ml-dtypes==0.2.0 \
        jax==0.4.14 \
        jaxlib==0.4.14 \
        jaxtyping==0.2.25 \
        lineax \
        qtl \
        equinox
RUN pip install lineax qtl
RUN pip install /opt/jaxqtl

CMD ["bash"]
