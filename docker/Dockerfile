################## BASE IMAGE ######################

FROM continuumio/miniconda3:latest

################## METADATA ######################
LABEL base.image="continuumio/miniconda3:latest"
LABEL version="1"
LABEL software="rna_lib_design"
LABEL software.version="0.1.0"
LABEL about.summary="designs RNA libraries and barcodes them "
LABEL about.tags="RNA bioinformatics,RNA structure"

################## MAINTAINER ######################
MAINTAINER Joe Yesselman <jyesselm@unl.edu>

################## INSTALLATION ######################

USER root

RUN /bin/bash -c "mkdir /data && mkdir /data2 && mkdir /data3 && mkdir /data4 "


# replace with installing vienna rna from github
#RUN /bin/bash -c "conda install -c bioconda viennarna"


RUN /bin/bash -c "mkdir -p /bin/rna_lib_design"
COPY . /bin/rna_lib_design
RUN /bin/bash -c "cd /bin/rna_lib_design && pip install ."
RUN /bin/bash -c "cd /bin/rna_lib_design && pytest"

WORKDIR /data

