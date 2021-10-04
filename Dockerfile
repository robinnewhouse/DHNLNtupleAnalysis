FROM rootproject/root:6.24.06-conda
COPY python /ntuple_analysis/python
COPY data /ntuple_analysis/data

# Build the image as root user
USER root
# Install the necessary packages
RUN pip install uproot3

# Add user "docker"
RUN useradd -ms /bin/bash docker

WORKDIR /ntuple_analysis/run

# Start the image as user "docker"
USER docker
CMD [ "/bin/bash" ]