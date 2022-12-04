FROM condaforge/mambaforge:4.14.0-0

COPY ["environment.yml", "./"]

RUN apt update && \
    DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends tzdata procps nano && \
    apt-get clean && \
    wget -q https://github.com/COMBINE-lab/salmon/releases/download/v1.9.0/salmon-1.9.0_linux_x86_64.tar.gz && \
    tar zxf salmon-1.9.0_linux_x86_64.tar.gz

ENV PATH=/salmon-1.9.0_linux_x86_64/bin:$PATH

RUN mamba env update --name base --file environment.yml