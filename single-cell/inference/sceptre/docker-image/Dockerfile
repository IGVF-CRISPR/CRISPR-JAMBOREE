# Use Ubuntu 20.04 LTS as base image
FROM ubuntu:20.04

# Set DEBIAN_FRONTEND to noninteractive to skip any interactive post-install configuration steps.
ENV DEBIAN_FRONTEND noninteractive

# Install R
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    gnupg2 \
    software-properties-common \
    dirmngr \
    make \
    curl \
    gcc \
    build-essential \
    libcurl4-openssl-dev \
    zlib1g-dev && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 0E98404D386FA1D9 && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt-get update && \
    apt-get install -y --no-install-recommends r-base

# Cleanup
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Bioconductor dependencies
RUN Rscript -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(version = "3.18")'
RUN Rscript -e 'BiocManager::install("MultiAssayExperiment")'
RUN Rscript -e 'BiocManager::install("SingleCellExperiment")'
RUN Rscript -e 'BiocManager::install("rhdf5")'

# GitHub dependencies
RUN Rscript -e 'install.packages("remotes", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'remotes::install_github("ilia-kats/MuData")'
RUN Rscript -e 'remotes::install_github("katsevich-lab/sceptre")'
RUN Rscript -e 'remotes::install_github("IGVF-CRISPR/sceptreIGVF")'