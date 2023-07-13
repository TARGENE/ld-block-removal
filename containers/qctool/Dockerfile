# Base image
FROM continuumio/miniconda3

# Update the system and install necessary packages
RUN apt-get update && apt-get install -y build-essential

# Install GCC 9
RUN apt-get install -y gcc-9 tar unzip liblapack3

# Set environment variables for Anaconda
ENV PATH="/opt/conda/bin:${PATH}"

# Copy environment.yaml file
COPY ./environment.yaml /env/

# Install additional packages if needed
RUN conda env update -n root --file /env/environment.yaml && conda clean -a

# Set working directory
WORKDIR /env/

# Download the file using wget
RUN wget https://www.well.ox.ac.uk/~gav/resources/qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz -O qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz

# Unzip the file
RUN mkdir qctool_v2.0.8 && tar xvzf qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz -C qctool_v2.0.8 --strip-components=1

# Remove old file
RUN rm qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz

# Add this to our path
ENV PATH="${PATH}:/env/qctool_v2.0.8/"