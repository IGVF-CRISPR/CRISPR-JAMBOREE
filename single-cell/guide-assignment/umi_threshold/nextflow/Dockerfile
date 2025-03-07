ARG PYTHON_VERSION=3.11.8-slim-bookworm

FROM python:${PYTHON_VERSION}

RUN apt-get update && apt-get install --no-install-recommends -qq wget ca-certificates make gcc g++ \
    git build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev zlib1g-dev libgit2-dev

WORKDIR /umi_threshold

COPY . .

RUN pip install "mudata==0.2.3" "scipy==1.15.2"
