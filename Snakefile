#########################################
# Snakemake pipeline for RNA-Seq analysis
#########################################

###########
# Libraries
###########
import pandas as pd

###############
# Configuration
###############

configfile: "/data/config.yaml" # where to find parameters
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]
THREADS = config["threads"]

########################
# Samples and conditions
########################

# read the tabulated separated table containing the sample, condition and fastq file information∂DE
units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

# create lists containing the sample names and conditions
SAMPLES = units.index.get_level_values('sample').unique().tolist()
samples = pd.read_csv(config["units"], dtype=str,index_col=0,sep="\t")
samplefile = config["units"]

###########################
# Input functions for rules
###########################

def sample_is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    if "fq2" not in samples.columns:
        return True
    else:
        return pd.isnull(samples.loc[(sample), "fq2"])

def get_fastq(wildcards):
    """ This function checks if the sample has paired end or single end reads
    and returns 1 or 2 names of the fastq files """
    if sample_is_single_end(wildcards.sample):
        return samples.loc[(wildcards.sample), ["fq1"]].dropna()
    else:
        return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_trimmed(wildcards):
    """ This function checks if sample is paired end or single end
    and returns 1 or 2 names of the trimmed fastq files """
    if sample_is_single_end(wildcards.sample):
        return WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz"
    else:
        return [WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz", WORKING_DIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq.gz"]

#################
# Desired outputs
#################
rule all:
    input:
        #expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES),
        RESULT_DIR + "counts.txt",
        # expand(WORKING_DIR + "mapped/{sample}.sorted.bam", sample=SAMPLES),
        expand(WORKING_DIR + 'stringtie/{sample}/transcript.gtf', sample=SAMPLES),
        RESULT_DIR + 'gene_FPKM.csv'
#        RESULT_DIR + "result.csv",
#        RESULT_DIR + "plotSelection.txt",
        #clusts = WORKING_DIR + "results/clusters.txt",
        #plots  = RESULT_DIR + "plots.pdf",
        #final = RESULT_DIR + "final.txt"
    message:
        "Job done! Removing temporary directory"

#######
# Rules
#######

##################################
# Fastp
##################################

rule fastp:
    input:
        get_fastq
    output:
        fq1  = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}.html"
    message:"trimming {wildcards.sample} reads"
    threads: THREADS
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sampleName = "{sample}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    run:
        if sample_is_single_end(params.sampleName):
            shell("fastp --thread {threads}  --html {output.html} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --in1 {input} --out1 {output.fq1} \
            2> {log}; \
            touch {output.fq2}")
        else:
            shell("fastp --thread {threads}  --html {output.html} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --detect_adapter_for_pe \
            --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2}; \
            2> {log}")

#########################
# RNA-Seq read alignement
#########################
if config["indexed"].upper().find("TRUE") == -1:
    if config["organism"].upper().find("HOMO") or config["organism"].upper().find("HUMAN") >= 0:
        rule ref_download:
            input:
                version = config["ref"]["hg_release_ver"] # release version
            output:
                directory(WORKING_DIR + 'genome')
            shell:"""
            mkdir {output} && \
            wget ftp://ftp.ensembl.org/pub/release-{input.version}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O genome.fa.gz | tar -zxvf -C {output} &&\
            wget ftp://ftp.ensembl.org/pub/release-{input.version}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{input.version}.gtf.gz -O genome.gtf.gz | tar -zxvf -C {output}"""
    elif config["organism"].upper().find("MUS") or config["organism"].upper().find("MOUSE") >= 0:
        rule ref_download:
            input:
                version = config["ref"]["mm_release_ver"] # release version
            output:
                directory(WORKING_DIR + 'genome')
            shell:"""
            mkdir {output} && \
            wget ftp://ftp.ensembl.org/pub/release-{input.version}/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz -O genome.fa.gz | tar -zxvf -C {output} &&\
            wget ftp://ftp.ensembl.org/pub/release-{input.version}/gtf/mus_musculus/Mus_musculus.GRCm38.{input.version}.gtf.gz -O genome.gtf.gz | tar -zxvf -C {output}"""


if config["aligner"].upper().find("HISAT2") >= 0:
    if config["indexed"].upper().find("TRUE") == -1:
        rule hisat_index:
            input:
                fasta = WORKING_DIR + "genome/genome.fa",
                gtf = WORKING_DIR + "genome/genome.gtf"
            output:
                [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
            message:
                "indexing genome"
            threads: THREADS
            shell:
                "cp scripts/make_grch38_tran.sh genome/ && sh genome/make_grch38_tran.sh"

    rule hisat_mapping:
        input:
            get_trimmed,
            indexFiles = [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
        output:
            bams = WORKING_DIR + "mapped/{sample}.sorted.bam",
            sum  = RESULT_DIR + "logs/{sample}_sum.txt",
            met  = RESULT_DIR + "logs/{sample}_met.txt"
        params:
            indexName = WORKING_DIR + "genome/genome",
            sampleName = "{sample}"
        message:
            "mapping reads to genome to bam files."
        threads: THREADS
        run:
            if sample_is_single_end(params.sampleName):
                shell("hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -q -x {params.indexName} \
                -U {input[0]} | samtools view -@ {threads} -Sb -F 4 | samtools sort -@ {threads} -o {output.bams}; \
                samtools index {output.bams}")
            else:
                shell("hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -q -x {params.indexName} \
                -1 {input[0]} -2 {input[1]} | samtools view -@ {threads} -Sb -F 4 | samtools sort -@ {threads} -o {output.bams}; \
                samtools index {output.bams}")
elif config["aligner"].upper().find("STAR") >= 0:
    if config["indexed"].upper().find("TRUE") == -1:
        rule star_index:
            input:
                fasta = WORKING_DIR + "genome/genome.fa", 
                gtf  = WORKING_DIR + "genome/genome.gtf"
            output:
                directory(WORKING_DIR + 'genome')
            message:
                "indexing genome"
            threads: THREADS
            shell:"""
            STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 100
            """
    rule star_mapping:
        input:
            get_trimmed
        output:
            bams = WORKING_DIR + "mapped/{sample}.sorted.bam",
        params:
            gtf = WORKING_DIR + 'genome/genome.gtf',
            index = WORKING_DIR + 'genome/'
        threads: THREADS
        shell:"""
        STAR --runThreadN {threads} \
        --sjdbOverhang 100 \
        --outSAMunmapped Within \
        --outputSAMtype BAM Unsorted \
        --outStd BAM_Unsorted \
        --sjdbGTFfile {params.gtf} \
        --genomeDir {params.index}
        --readFilesIn {input} \
        --readFilesCommand zcat \
        | samtools sort -@ {threads} -O bam -o {output.bam}
        """


#########################################
# Get table containing the RPKM or FPKM
#########################################

rule stringtie:
    input:
        bams = WORKING_DIR + "mapped/{sample}.sorted.bam"
    output:
        r1 = WORKING_DIR + "stringtie/{sample}/transcript.gtf",
        r2 = WORKING_DIR + "stringtie/{sample}/gene_abundances.tsv",
        r3 = WORKING_DIR + "stringtie/{sample}/cov_ref.gtf"
    message:
        "assemble RNA-Seq alignments into potential transcripts."
    threads: THREADS
    params:
        gtf = WORKING_DIR + "genome/genome.gtf"
    shell:
        "stringtie -p {threads} -G {params.gtf} --rf -e -B -o {output.r1} -A {output.r2} -C {output.r3} --rf {input.bams}"

rule create_PKM_table:
    input:
        WORKING_DIR
    output:
        r1 = RESULT_DIR + "gene_FPKM.csv",
        outdir = RESULT_DIR
    params:
        dataset = config["merge_PKM"]["organism"]
    message:
        "create gene and transcript FPKM(if single-end reads, RPKM)."
    conda:
        "envs/merge_fpkm.yaml"
    shell:
        "Rscript scripts/merge_RFPKM.r -i {input} -o {output.outdir} -d {params.dataset}"

#########################################
# Get table containing the raw counts
#########################################

rule create_counts_table:
    input:
        bams = expand(WORKING_DIR + "mapped/{sample}.sorted.bam", sample = SAMPLES),
        gff  = WORKING_DIR + "genome/genome.gtf"
    output:
        RESULT_DIR + "counts.txt"
    message:
        "create read count talbe"
    threads: THREADS
    shell:
        "featureCounts -T {threads} -a {input.gff} -t exon -g gene_id -o {output} {input.bams}"
