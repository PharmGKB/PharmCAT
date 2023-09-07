# syntax=docker/dockerfile:1
#
# Base Dockerfile for running Jekyll
#
FROM ruby

ENV NODE_MAJOR_VERSION 20

RUN mkdir -p /etc/apt/keyrings
RUN curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg
RUN echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_${NODE_MAJOR_VERSION}.x nodistro main" | tee /etc/apt/sources.list.d/nodesource.list

# apt-utils line due to https://github.com/phusion/baseimage-docker/issues/319
RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils apt-transport-https gnupg && \
    apt-get -y upgrade && \
    apt-get -y install bzip2 build-essential wget nodejs

RUN npm install --location=global surge
RUN npm install --location=global npm

COPY Makefile Gemfile Gemfile.lock ./
RUN make setup && \
    rm -f Makefile Gemfile Gemfile.lock

WORKDIR /PharmCAT
