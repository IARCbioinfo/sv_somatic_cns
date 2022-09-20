################# BASE IMAGE #####################
FROM continuumio/miniconda3:4.7.12
##site to test docker configuration files
# https://labs.play-with-docker.com/
################## METADATA #######################
LABEL base_image="continuumio/miniconda3"
LABEL version="4.7.12"
LABEL software="sv_somatic_cns-nf"
LABEL software.version="1.1"
LABEL about.summary="Container image containing all requirements for sv_somatic_cns"
LABEL about.home="https://github.com/IARCbioinfo/sv_somatic_cns-nf"
LABEL about.documentation="https://github.com/IARCbioinfo/sv_somatic_cns-nf/README.md"
LABEL about.license_file="https://github.com/IARCbioinfo/sv_somatic_cns-nf/LICENSE.txt"
LABEL about.license="MIT"

################## MAINTAINER ######################
MAINTAINER Nicolas Alcala <alcalan@iarc.who.int>
################## INSTALLATION ######################
COPY environment.yml /
#the next command install the ps command needed by nexflow to collect run metrics
RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps
RUN conda env create -n sv_somatic_cns -f /environment.yml && conda clean -a
# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/sv_somatic_cns/bin:$PATH
# Dump the details of the installed packages to a file for posterity
RUN conda env export --name sv_somatic_cns > sv_somatic_cns-v1.1.yml
