# Use an official Python runtime as a parent image
FROM python:3.7-slim

# Set the working directory in the container
WORKDIR /app

# Copy the requirements directly into the Dockerfile for better readability
RUN pip install \
    pytest==6.2.5 \
    pandas==1.2.5 \
    numpy==1.20.3 \
    scipy==1.6.3 \
    setuptools==56.2.0 \
    typer==0.4.0 \
    PyYAML==5.4.1 \
    ray==1.8.0 \
    scikit-learn==0.24.2 \
    tensorly==0.6.0 \
    requests==2.26.0 \
    seaborn==0.11.1 \
    matplotlib==3.4.3 \
    networkx==2.6.3 \
    scanpy==1.7.2 \
    protobuf==3.20.*

# Install scTenifoldpy package
RUN pip install scTenifoldpy

# Install the missing importlib_metadata for Python 3.7 compatibility
RUN pip install importlib-metadata

# Modify umap package to use "importlib_metadata" instead of "importlib.metadata"
RUN sed -i 's/importlib.metadata/importlib_metadata/' /usr/local/lib/python3.7/site-packages/umap/__init__.py

# Default command
#CMD ["python3"]

