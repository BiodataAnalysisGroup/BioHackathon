# Start from an official Python 3.6 image
FROM python:3.6

# Set a working directory in the container
WORKDIR /usr/src/app

# Install required packages
RUN pip install --no-cache-dir adjustText==0.7.3 trvae

# Copy your application files (if any) into the container
# COPY . .

# Define the default command to run when starting the container
CMD ["python"]
