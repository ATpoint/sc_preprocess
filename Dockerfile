FROM condaforge/mambaforge:4.10.3-1

COPY ["environment.yml", "./"]

RUN apt update && \
    DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends tzdata wget nano && \
    apt-get clean

RUN mamba env update --name base --file environment.yml

CMD echo "This is the atac_chip_preprocess image with versions:" $(mamba --version | tr "\n" "\ ")