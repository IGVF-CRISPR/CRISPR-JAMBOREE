ARG PYTHON_VERSION=3.11.8-slim-bookworm

FROM python:${PYTHON_VERSION}

RUN apt-get update && apt-get install --no-install-recommends -qq wget ca-certificates make gcc g++ \
    git build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev zlib1g-dev libgit2-dev

WORKDIR /cleanser

COPY . .

RUN pip install "git+https://github.com/Gersbachlab-Bioinformatics/CLEANSER.git@197ebea49ae7ec02e45cbea2b3420e9b26fc4a81"; install_cmdstan
