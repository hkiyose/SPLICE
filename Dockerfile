FROM ubuntu:20.04

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update \
  && apt-get -y --no-install-recommends install \
    python3 \
    minimap2 \
    vim \
  && apt-get -y clean \
  && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /SPLICE/data
COPY ./ /SPLICE

WORKDIR /SPLICE

