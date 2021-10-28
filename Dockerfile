FROM mambaorg/micromamba:0.15.3

USER root 

COPY --chown=micromamba:micromamba environment.yml environment.yml

RUN apt update && \
    DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends tzdata wget nano && \ 
    apt-get clean && \
    micromamba install -y -n base -f environment.yml && \   
    micromamba clean --all --yes