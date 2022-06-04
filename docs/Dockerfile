# syntax=docker/dockerfile:1
#
# Base Dockerfile for running Jekyll
#
FROM ruby

# apt-utils line due to https://github.com/phusion/baseimage-docker/issues/319
RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils apt-transport-https gnupg && \
    apt-get -y upgrade && \
    apt-get -y install bzip2 build-essential wget

RUN curl -fsSL https://deb.nodesource.com/setup_16.x | bash -
RUN apt-get install -y nodejs
RUN npm install --location=global surge
RUN npm install --location=global npm

WORKDIR /PharmCAT

