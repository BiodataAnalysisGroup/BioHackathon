# Use an official Python image as a base
FROM python:3.8-slim

# Set environment variables to make Python run in unbuffered mode and define frontend as noninteractive
ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential libcurl4-openssl-dev zlib1g-dev libfftw3-dev libc++-dev git wget hdf5-tools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory in the container
WORKDIR /app

# Clone the GitHub repository
RUN git clone https://github.com/BhattacharyaLab/scVIDR.git .

# Install virtualenv and create a virtual environment
RUN pip install virtualenv
RUN python -m virtualenv env

# Set the PATH to use the virtualenv for all subsequent commands
ENV PATH="/app/env/bin:$PATH"

# Check if the virtual environment is being used
RUN which python
RUN python --version
RUN which pip
RUN pip --version

# First, install PyTorch and related packages
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir torch==1.8.1+cu111 torchaudio torchvision -f https://download.pytorch.org/whl/torch_stable.html

# Then install Python dependencies including geomloss
RUN pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir geomloss==0.2.5

# Verify that packages are installed in the virtual environment
#RUN pip list

