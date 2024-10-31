# Start with the official Python 3.9.5 base image
FROM python:3.9.5

# Set environment variables to prevent Python from writing .pyc files and to buffer stdout/stderr for logging
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Install necessary system dependencies
RUN apt-get update && apt-get install -y \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Install specific versions of packages using pip
RUN pip install --no-cache-dir \
    setuptools==49.6.0 \
    numpy==1.19.5 \
    matplotlib==3.3.3 \
    scanpy==1.8.1 \
    jupyterlab==3.0.16 \
    ipython==7.18.1 \
    anndata==0.7.6 \
    ml-collections==0.1.0 \
    pandas==1.2.5 \
    scipy==1.8.1 \
    scikit-learn==1.1.1 \
    seaborn==0.11.1 \
    sphinx==4.2.0 \
    sphinx-press-theme==0.8.0 \
    torch==1.11.0

# Clone and install cellot from GitHub
RUN git clone https://github.com/bunnech/cellot.git /app/cellot && \
    pip install /app/cellot

# Set the default command to open a bash shell
CMD ["/bin/bash"]
