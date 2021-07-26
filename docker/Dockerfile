FROM ubuntu:18.04
MAINTAINER Jinwoo jeong <wjd5480@alopaxalgo.com>
LABEL org.opencontainers.image.source https://github.com/williamjeong2/snakemake_RNA-seq

SHELL ["/bin/bash", "-c"]
ARG DEBIAN_FRONTEND=noninteractive

# RUN sed -i.bak -re "s/([a-z]{2}.)?archive.ubuntu.com|security.ubuntu.com/mirror.kakao.com/g" /etc/apt/sources.list
RUN apt-get update \
    && apt-get install -y python3 python3-pip wget git less vim \
    && pip3 install pandas awscli
RUN apt-get install -y libssl-dev libcurl4-openssl-dev libxml2-dev curl libcairo2-dev libxt-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libproj-dev libproj-dev libgdal-dev
RUN echo 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
RUN apt-get update && apt-get install -y r-base r-base-dev r-recommended

####################
# install R
####################

RUN apt-get update -y && \
    apt-get install -y apt-utils 2>&1 | grep -v "debconf: delaying package configuration, since apt-utils is not installed" && \
    apt-get install -y --no-install-recommends locales && \
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
RUN Rscript -e 'install.packages("BiocManager", repos="http://cran.rstudio.com/", dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install("DESeq2", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("ggplot2", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("clusterProfiler", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("DOSE", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("KEGG.db", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("org.Mm.eg.db", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("org.Hs.eg.db", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("pheatmap", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("genefilter", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("RColorBrewer", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("GO.db", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("topGO", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("gage", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("ggsci", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("curl", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("biomaRt", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("EnhancedVolcano", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("fgsea", suppressUpdates = TRUE)'
RUN Rscript -e 'BiocManager::install("DEGreport", suppressUpdates = TRUE)'

# install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

# pull git repo and install conda env for rna sequencing
WORKDIR /home/
RUN git clone https://github.com/WilliamJeong2/snakemake_RNA-seq.git
RUN conda install -c conda-forge mamba
RUN mamba env create -f snakemake_RNA-seq/envs/global_env.yaml

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
RUN chmod -R 777 /home/snakemake_RNA-seq/
RUN /bin/bash -c "source /root/.bashrc"

WORKDIR /home/snakemake_RNA-seq