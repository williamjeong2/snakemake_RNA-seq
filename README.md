# snakemake_RNA-seq
This repo is forked from [KoesGroup/Snakemake_hisat-DESeq](https://github.com/KoesGroup/Snakemake_hisat-DESeq)

A snakemake pipeline for the analysis of RNA-seq data that makes use of [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) and [Stringtie](https://ccb.jhu.edu/software/stringtie/).

# Aim
To align. count, normalize counts and compute DEG between conditions using single-end or paired-end Illumina RNA-seq data.

# Content
- `Snakefile`:
- `config.yaml`:
- `data/`:
- `envs/`:
- `samples.tsv`:

# Usage

## Download or clone the Github repository
You will need a local copy of the `Snakemake_RNA-seq` on your machine.
You can either:
1. use git in the shell: `git clone git@github.com:WilliamJeong2/snakemake_RNA-seq.git
2. click on "Clone or download" and select `download`

## Installing and activating a virtual environment
First, you need to create an environment where `Snakemake` and the python `pandas` package and something else will be installed. To do that, we will use the conda package manager.
1. Create a virtual environment named `rna-seq` using the `global_env.yaml` file with the folling command: `conda env create --name rna-seq --file envs/global_env.yaml`
2. Activate this virtual environment with source activate rna-seq

The Snakefile will then take care of installing and loading the packages and softwares required by each step of the pipeline.

## Configuration file
Make sure you have changed the parameters in the `config.yaml` file that specifies where to find the sample data file, the genomic and transcriptomic referece fasta files to use and the parameters for certains rules etc.
This file is used so the `Snakefile` does not need to be changed when locations or parameters need to be changed.

## Snakemake execution
The Snakemake pipeline/workflow management system reads a master file (often called `Snakefile`) to list the steps to be executed and defining their order. It has many rich features. Read more [here](https://snakemake.readthedocs.io/en/stable/)

## Dry run (recommend)
From the folder containing the `Snakefile`, use the command `snakemake --use-conda -np` to perform a dry run that prints out the rules and commands.

## Real run
Simply type `snakemake --use-conda` and provide the number of cores with `--cores 60` for the cores for instance. 

# output files
- the RNA-seq read alignment files : ***.bam** (in temp dir)
- the fastqc report files : ***.html** (in results dir)
- the unscaled RNA-seq read counts : **counts.txt** (in results dir)
- gene/transcript level RPKM or FPKM : **gene_FPKM.csv** (in results dir)
