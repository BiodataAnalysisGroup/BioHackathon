# Use the official Python 3.8 image from the Docker Hub
FROM python:3.8

# Set the working directory inside the container
WORKDIR /app

# Copy any necessary files (like requirements.txt) if you have them
# COPY requirements.txt .

# Install scpram package from the specified PyTorch CPU index URL
RUN pip install --no-cache-dir scpram --extra-index-url https://download.pytorch.org/whl/cpu

# Optional: Copy your application code (if you have any)
# COPY . .

# Command to run your application (modify as needed)
# CMD ["python", "your_script.py"]

# Expose any ports if necessary (uncomment and specify port)
# EXPOSE 5000

