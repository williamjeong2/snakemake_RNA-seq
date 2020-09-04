FROM ubuntu:18.04
MAINTAINER Jinwoo jeong <jinwoo5480@snu.ac.kr>

SHELL ["/bin/bash", "-c"]
ARG DEBIAN_FRONTEND=noninteractive

# RUN sed -i.bak -re "s/([a-z]{2}.)?archive.ubuntu.com|security.ubuntu.com/mirror.kakao.com/g" /etc/apt/sources.list
RUN apt-get update \
    && apt-get install -y python3 python3-pip wget git less vim \
    && pip3 install pandas
RUN apt-get install -y libxml2-dev curl libssl-dev libgsl0-dev libcairo2-dev libhdf5-dev libcurl4-openssl-dev
RUN echo 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
RUN apt-get update && apt-get install -y r-base r-base-dev r-recommended

####################
# install R
####################

RUN apt-get update && apt-get install -y --no-install-recommends locales apt-utils && \
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.UTF-8 && \
    LC_ALL=en_US.UTF-8 && \
    LANG=en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8 && \
    TERM=xterm && \
    apt-get update

RUN Rscript -e 'install.packages("optparse", repos="http://cran.rstudio.com/", dependencies=TRUE)'
RUN Rscript -e 'install.packages("data.table", repos="http://cran.rstudio.com/", dependencies=TRUE)'
RUN Rscript -e 'install.packages("tidyr", repos="http://cran.rstudio.com/", dependencies=TRUE)'
RUN Rscript -e 'install.packages("dplyr", repos="http://cran.rstudio.com/", dependencies=TRUE)'
RUN Rscript -e 'install.packages("BiocManager", repos="http://cran.rstudio.com/", dependencies=TRUE);'
RUN Rscript -e 'BiocManager::install(version = "3.11"); '
RUN Rscript -e 'BiocManager::install("DESeq2"); library(DESeq2);'
RUN Rscript -e 'BiocManager::install("ggplot2"); library(ggplot2);'
RUN Rscript -e 'BiocManager::install("clusterProfiler"); library(clusterProfiler);'
RUN Rscript -e 'BiocManager::install("DOSE"); library(DOSE);'
RUN Rscript -e 'BiocManager::install("KEGG.db"); library(KEGG.db);'
RUN Rscript -e 'BiocManager::install("org.Mm.eg.db"); library(org.Mm.eg.db);'
RUN Rscript -e 'BiocManager::install("org.Hs.eg.db"); library(org.Hs.eg.db);'
RUN Rscript -e 'BiocManager::install("pheatmap"); library(pheatmap);'
RUN Rscript -e 'BiocManager::install("genefilter"); library(genefilter);'
RUN Rscript -e 'BiocManager::install("RColorBrewer"); library(RColorBrewer);'
RUN Rscript -e 'BiocManager::install("GO.db"); library(GO.db);'
RUN Rscript -e 'BiocManager::install("topGO"); library(topGO);'
RUN Rscript -e 'BiocManager::install("gage"); library(gage);'
RUN Rscript -e 'BiocManager::install("ggsci"); library(ggsci);'
RUN Rscript -e 'BiocManager::install("curl"); library(curl);'
RUN Rscript -e 'BiocManager::install("biomaRt"); library(biomaRt)'

# error occured
# RUN Rscript -e 'BiocManager::install("reactome.db");'
# RUN Rscript -e 'BiocManager::install("ReactomePA"); library(ReactomePA);'

# install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

# pull git repo and install conda env for rna sequencing
WORKDIR /home/
RUN git clone https://github.com/WilliamJeong2/snakemake_RNA-seq.git
RUN conda env create -n rnaseq -f snakemake_RNA-seq/envs/docker.yml

# Clean up
RUN cd / && \
  rm -rf /tmp/* && \
  apt-get autoremove -y && \
  apt-get autoclean -y && \
  rm -rf /var/lib/apt/lists/* && \
  apt-get clean

RUN cat /home/snakemake_RNA-seq/docker/conda_init.sh >> /root/.bashrc
RUN echo "conda activate rnaseq" >> /root/.bashrc

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "base", "/bin/bash", "-c"]
CMD /bin/bash -c "source /root/.bashrc"

WORKDIR /home/snakemake_RNA-seq