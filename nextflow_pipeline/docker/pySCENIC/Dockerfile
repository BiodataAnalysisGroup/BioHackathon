FROM python:3.10.6-slim-bullseye AS compile-image

ENV DEBIAN_FRONTEND=noninteractive

# Install git to be able to clone the pySCENIC repository.
RUN apt-get update && \
    apt-get -y --no-install-recommends install git

# Clone the pySCENIC repository from GitHub.
RUN git clone https://github.com/aertslab/pySCENIC.git /tmp/pySCENIC

# Install dependencies from the pySCENIC requirements.txt.
RUN pip install --no-cache-dir -r /tmp/pySCENIC/requirements.txt

# Install pySCENIC from the cloned repository.
RUN cd /tmp/pySCENIC && \
    pip install . && \
    cd .. && rm -rf pySCENIC

FROM python:3.10.6-slim-bullseye AS build-image

RUN apt-get -y update && \
    apt-get -y upgrade && \
    apt-get -y --no-install-recommends install \
        # Need to run ps
        procps \
        libxml2 \
        less && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

# Copy the installed packages from the compile-image
COPY --from=compile-image /usr/local/lib/python3.10/site-packages /usr/local/lib/python3.10/site-packages
COPY --from=compile-image /usr/local/bin /usr/local/bin

# Ensure the path is set correctly for the installed binaries
ENV PATH="/usr/local/bin:$PATH"

EXPOSE 8787
