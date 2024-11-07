# Use the official Python 3.7 image from Docker Hub
FROM python:3.7

# Set the working directory in the container
WORKDIR /app

# Upgrade pip and install necessary dependencies
RUN pip install --upgrade pip \
    && pip install numpy \
    && pip install pandas==1.3.5 \
    && pip install cython \
    && pip install h5py \
    && pip install igraph \
    && pip install louvain \
    && pip install scanpy \
    && pip install anndata \
    && pip install joblib \
    && pip install scikit-learn \
    && pip install velocyto \
    && pip install pyarrow \
    && pip install seaborn \
    && pip install genomepy \
    && pip install numba \
    && pip install celloracle

# If you have a Python script to run, copy it into the container
# COPY your_script.py .

# Command to run your application
# Replace 'your_script.py' with your actual script name if needed
# CMD ["python", "your_script.py"]

