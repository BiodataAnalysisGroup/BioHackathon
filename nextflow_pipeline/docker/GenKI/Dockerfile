# Use the official Python 3.10 image from Docker Hub
FROM python:3.10

# Set the working directory
WORKDIR /app

# Install Git
RUN apt-get update && apt-get install -y git

# Install the required Python packages
RUN pip install --no-cache-dir git+https://github.com/yjgeno/GenKI.git
RUN pip install --no-cache-dir torch torchaudio torch_geometric
RUN pip install --no-cache-dir tensorboard

# Set the default command to bash
#CMD ["bash"]
