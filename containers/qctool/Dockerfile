# Base image
FROM continuumio/anaconda

# Update the system
RUN apt-get update && apt-get upgrade -y

# Install GCC 5.5.0
RUN apt-get install -y gcc-5=5.5.0-12ubuntu1

# Set environment variables for Anaconda
ENV PATH="/opt/conda/bin:${PATH}"

# Install tar
RUN apt-get install -y tar

# Copy environment.yaml file
COPY container/environment.yaml .

# Install additional packages if needed
RUN conda env create --file environment.yaml

# Download the file using wget
RUN wget https://www.well.ox.ac.uk/~gav/resources/qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz -O qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz

# Unzip the file
RUN tar xvzf qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz

# Remove old file
RUN rm qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz

# Copy your project files to the container (if necessary)
COPY qctool*/qctool /usr/bin

