# Use the official Python image with version 3.9.18
FROM python:3.9.18

# Set the working directory in the container
WORKDIR /app

# Install DSPIN and other required Python packages
RUN pip install --no-cache-dir \
    dspin \
    anndata==0.10.3 \
    matplotlib==3.8.2 \
    scanpy==1.9.6 \
    tqdm==4.65.0 \
    leidenalg==0.10.1 \
    igraph==0.10.8

# Start a shell when the container runs
CMD ["/bin/bash"]
