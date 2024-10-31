# Use the official Python 3.9 image as the base image
FROM python:3.8

# Set the working directory inside the container
WORKDIR /app

# Install system dependencies for building Python packages
RUN apt-get update && \
    apt-get install -y git build-essential python3-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Clone the repository
RUN git clone https://github.com/JaneJiayiDong/scPreGAN.git

# Set the working directory to the cloned repository
WORKDIR /app/scPreGAN

# Upgrade pip and setuptools to ensure compatibility with complex dependencies
RUN pip install --upgrade pip setuptools wheel

# Install numpy and scipy first, as they often have complex build requirements
RUN pip install numpy==1.19.2 scipy==1.5.3

# Install the remaining Python packages
RUN pip install anndata==0.7.4 h5py==2.10.0 pandas==1.1.3 scanpy==1.6.0 scapy==2.4.4 \
    scikit-learn==0.23.2 tensorboard==1.15.0 tensorboardX==2.2 torch==1.8.1 \
    torchaudio==0.8.1 torchvision==0.9.1

# (Optional) Expose any ports if necessary
# EXPOSE <port_number>

# Command to run the application (adjust as necessary)
CMD ["python", "your_script.py"]
