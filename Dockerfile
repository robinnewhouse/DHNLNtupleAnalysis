FROM continuumio/miniconda3
COPY python /ntuple_analysis/python
COPY data /ntuple_analysis/data

# Build the image as root user
USER root
# Install the necessary packages
RUN conda install -c conda-forge root -y && \
    pip install uproot3

# Add user "docker"
RUN useradd -ms /bin/bash docker

WORKDIR /ntuple_analysis/run

# Start the image as user "docker"
USER docker