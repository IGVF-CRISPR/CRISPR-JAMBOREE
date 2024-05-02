ARG PYTHON_VERSION=3.11.8-slim-bookworm

FROM python:${PYTHON_VERSION}

RUN apt-get update && apt-get install --no-install-recommends -qq wget ca-certificates make gcc g++ \
    r-base-core git build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev zlib1g-dev libfontconfig-dev \
    libgit2-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev gfortran \
    libblas3 libblas-dev liblapack3 liblapack-dev

WORKDIR /cleanser

COPY . .

RUN Rscript -e 'install.packages("tidyverse", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("devtools", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install("MultiAssayExperiment")'; Rscript -e 'BiocManager::install("SingleCellExperiment")'; Rscript -e 'BiocManager::install("SummarizedExperiment")'; Rscript -e 'BiocManager::install("MuData")'
RUN Rscript -e 'devtools::install_github("katsevich-lab/sceptre")'; Rscript -e 'devtools::install_github("IGVF-CRISPR/sceptreIGVF")'

RUN pip install mudata; pip install "git+https://github.com/siyansusan/CLEANSER.git"; install_cmdstan
