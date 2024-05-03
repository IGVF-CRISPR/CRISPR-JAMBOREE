ARG PYTHON_VERSION=3.11.8-slim-bookworm

FROM python:${PYTHON_VERSION}

RUN apt-get update && apt-get install --no-install-recommends -qq wget ca-certificates make gcc g++ \
    r-base-core git build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev zlib1g-dev libfontconfig-dev \
    libgit2-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev gfortran \
    libblas3 libblas-dev liblapack3 liblapack-dev

WORKDIR /cleanser

COPY . .

RUN Rscript -e 'install.packages("dplyr", repos = "http://cran.us.r-project.org")'; Rscript -e 'install.packages("Matrix", repos = "http://cran.us.r-project.org")'; \
    Rscript -e 'install.packages("tibble", repos = "http://cran.us.r-project.org")'; Rscript -e 'install.packages("crayon", repos = "http://cran.us.r-project.org")'; \
    Rscript -e 'install.packages("ggplot", repos = "http://cran.us.r-project.org")'; Rscript -e 'install.packages("scales", repos = "http://cran.us.r-project.org")'; \
    Rscript -e 'install.packages("ondisc", repos = "http://cran.us.r-project.org")'; Rscript -e 'install.packages("data.table", repos = "http://cran.us.r-project.org")'; \
    Rscript -e 'install.packages("purrr", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'; Rscript -e 'BiocManager::install(version = "3.18")'
RUN Rscript -e 'install.packages("remotes", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install("MultiAssayExperiment")'; Rscript -e 'BiocManager::install("SingleCellExperiment")'; Rscript -e 'BiocManager::install("SummarizedExperiment")';
RUN Rscript -e 'remotes::install_github("grimbough/rhdf5@537cd703a3a27683f3e5224ead62d00972d426d6")'; Rscript -e 'remotes::install_github("ilia-kats/MuData@4f4f4aeffc83343ed4d15c3db8cda44684306688")'
RUN Rscript -e 'remotes::install_github("katsevich-lab/sceptre@8abf87aaa390dcf5188fbf5bbe7e700b05e53697")'; Rscript -e 'remotes::install_github("IGVF-CRISPR/sceptreIGVF@87a35d6736a5ddbf9c2975323d55e4423bd88219")'

RUN pip install mudata==0.2.3; pip install "git+https://github.com/siyansusan/CLEANSER.git@4f3c331058d2e6cc237666523f29b75c768443d4"; install_cmdstan
