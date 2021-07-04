FROM rocker/rstudio:4.0.2

COPY ["environment.yml", "./"]

RUN apt update && \
    apt install --assume-yes wget nano

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH /opt/miniconda3/condabin/:/opt/miniconda3/bin:$PATH

RUN conda init bash && \
    conda install mamba -n base -c conda-forge && \
    mamba env update --name base -f environment.yml && \
    mamba clean -afy
