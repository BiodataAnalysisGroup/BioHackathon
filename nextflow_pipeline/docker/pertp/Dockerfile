# Use the official Python 3.10.15 image as the base
FROM python:3.10.15

# Set the working directory in the container
WORKDIR /app

# Copy the project files into the container
COPY . .

# Upgrade pip, wheel, and setuptools to ensure compatibility
RUN pip install --upgrade pip wheel setuptools

# Install required build system and main dependencies
RUN pip install hatchling && \
    pip install pertpy

# Install additional dependencies
#RUN pip install rich scanpy[leiden] muon requests scikit-misc scipy scvi-tools adjusttext decoupler \
#                ott-jax sparsecca openpyxl pubchempy pyarrow blitzgsea lamin_utils

# Uncomment the lines below if you want to install additional optional dependencies
# RUN pip install toytree arviz ete3 pyqt5
# RUN pip install formulaic pydeseq2
# RUN pip install pre-commit
# RUN pip install docutils sphinx furo myst-nb myst-parser sphinxcontrib-bibtex sphinx-gallery \
#                sphinx-autodoc-typehints sphinx-copybutton sphinx-remove-toctrees sphinx-design \
#                sphinx-last-updated-by-git sphinx-automodapi sphinxext-opengraph pygments \
#                pyenchant nbsphinx nbsphinx-link ipykernel ipython
# RUN pip install pytest coverage

# Start a shell when the container runs
#CMD ["/bin/bash"]
