FROM continuumio/miniconda3:24.9.2-0

RUN conda install -y -n base -c conda-forge -c bioconda python=3.11 pip procps-ng gxx git python-annoy datrie && \
    conda clean --all --yes

RUN git clone https://github.com/aertslab/scenicplus 

RUN cd scenicplus && pip install .

