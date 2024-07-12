ARG PYTHON_VERSION=3.11.8-slim-bookworm

FROM python:${PYTHON_VERSION}

RUN apt-get update && apt-get install --no-install-recommends -qq wget ca-certificates make gcc g++ \
    git build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev zlib1g-dev libgit2-dev

WORKDIR /cleanser

COPY . .

RUN pip install mudata==0.2.3; pip install "git+https://github.com/siyansusan/CLEANSER.git@4f3c331058d2e6cc237666523f29b75c768443d4"; install_cmdstan
