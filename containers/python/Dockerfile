FROM nfcore/base:1.12.1

RUN apt-get update && apt-get -y install unzip

COPY ./environment.yaml /env/
RUN conda env update -n root --file /env/environment.yaml && conda clean -a