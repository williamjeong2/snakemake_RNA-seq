FROM ubuntu:18.04
MAINTAINER Jinwoo jeong <jinwoo5480@snu.ac.kr>

# ENV R_BASE_VERSION 4.0.0

RUN sed -i.bak -re "s/([a-z]{2}.)?archive.ubuntu.com|security.ubuntu.com/mirror.kakao.com/g" /etc/apt/sources.list
RUN apt-get update \
    && apt-get install -y python3 python3-pip wget git less vim \
    && pip3 install pandas 
RUN apt-get install -y libxml2-dev curl libcurl4-gnutls-dev libssl-dev
RUN echo 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y r-base r-base-dev r-recommended

####################
# install R
####################

RUN apt-get update && apt-get install -y --no-install-recommends locales && \
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.UTF-8 && \
    LC_ALL=en_US.UTF-8 && \
    LANG=en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8 && \
    TERM=xterm && \
    apt-get install -y --no-install-recommends \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libgsl0-dev \
    libcairo2-dev \
    libhdf5-dev \
    && apt-get update

RUN Rscript -e 'install.packages("optparse", repos="http://cran.us.r-project.org", dependencies=TRUE)'
RUN Rscript -e 'install.packages("data.table", repos="http://cran.us.r-project.org", dependencies=TRUE)'
RUN Rscript -e 'install.packages("tidyr", repos="http://cran.us.r-project.org", dependencies=TRUE)'
RUN Rscript -e 'install.packages("dplyr", repos="http://cran.us.r-project.org", dependencies=TRUE)'
RUN Rscript -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org", dependencies=TRUE); \
                BiocManager::install(version = "3.11");  \
                BiocManager::install("DESeq2"); library(DESeq2)\
                BiocManager::install("ggplot2"); library(ggplot2)\
                BiocManager::install("clusterProfiler"); library(clusterProfiler)\
                BiocManager::install("ReactomePA"); library(ReactomePA)\
                BiocManager::install("DOSE"); library(DOSE)\
                BiocManager::install("KEGG.db"); library(KEGG.db)\
                BiocManager::install("org.Mm.eg.db"); library(org.Mm.eg.db)\
                BiocManager::install("org.Hs.eg.db"); library(org.Hs.eg.db)\
                BiocManager::install("pheatmap"); library(pheatmap)\
                BiocManager::install("genefilter"); library(genefilter)\
                BiocManager::install("RColorBrewer"); library(RColorBrewer)\
                BiocManager::install("GO.db"); library(GO.db)\
                BiocManager::install("topGO"); library(topGO)\
                BiocManager::install("gage"); library(gage)\
                BiocManager::install("ggsci"); library(ggsci)\
                BiocManager::install("biomaRt"); library(biomaRt)'

####################
# pull git repo and install conda env for rna sequencing
####################

ENV PATH /opt/conda/bin:${PATH}

RUN mkdir /opt/tmp && cd /opt/tmp/ \
    && wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm Miniconda3-latest-Linux-x86_64.sh \
    && git clone https://github.com/WilliamJeong2/snakemake_RNA-seq.git
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc
RUN cd opt/tmp/snakemake_RNA-seq \
    && conda env create -n rnaseq -f envs/docker.yml \
    && /bin/bash -c "source activate rnaseq"
RUN echo "conda activate rnaseq" > ~/.bashrc
ENV PATH /opt/conda/envs/rnaseq/bin:${PATH}

## Clean up
RUN cd / && \
    rm -rf /tmp/* && \
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean

WORKDIR /opt/tmp/snakemake_RNA-seq
