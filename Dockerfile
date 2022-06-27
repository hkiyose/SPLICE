FROM ubuntu:20.04

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update \
  && apt-get -y --no-install-recommends install \
    python3 \
    python3-pip \
    minimap2 \
    vim \
    luigi \
  && apt-get -y clean \
  && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /SPLICE/data
RUN pip3 install luigi
COPY ./ /SPLICE

WORKDIR /SPLICE

