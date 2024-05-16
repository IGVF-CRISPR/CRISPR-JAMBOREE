# Use Ubuntu 20.04 LTS as base image
FROM rocker/r-ver:4.4.0

# Set DEBIAN_FRONTEND to noninteractive to skip any interactive post-install configuration steps.
ENV DEBIAN_FRONTEND noninteractive

# System dependencies
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
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# R dependencies
RUN Rscript -e 'install.packages("renv", repos = "http://cran.us.r-project.org")'
COPY renv.lock renv.lock
RUN Rscript -e 'renv::restore()'
