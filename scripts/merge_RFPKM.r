#!/usr/bin/env Rscript
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table",repos = "http://cran.us.r-project.org")
library(data.table)
if (!requireNamespace("tidyr", quietly = TRUE))
  install.packages("tidyr",repos = "http://cran.us.r-project.org")
library(tidyr)
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr",repos = "http://cran.us.r-project.org")
library(dplyr)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repos = "http://cran.us.r-project.org")
# update.packages(ask = FALSE,repos = "http://cran.us.r-project.org")
#BiocManager::install(version = "3.11")
if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")
# BiocManager::install("digest")
library(biomaRt)
library(optparse)

option_list = list(
    make_option(c("-d", "--dataset"), type="character", default="hsapiens_gene_ensembl", metavar="character"),
    make_option(c("-i", "--indir"), type="character", default="temp/", metavar="character"), # wpath : ~/temp/
    make_option(c("-o", "--outdir"), type="character", default="results/", metavar="character"), # result_path : ~/results/
    make_option("--trans", type="character", default="gmIDs_hg.tsv", metavar="character"),
    make_option("--gene", type="character", default="gmIDs_g_hg.tsv", metavar="character")
)

# parse the command-line arguments and pass them to a list called 'opt'
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#mart = useEnsembl("ENSEMBL_MART_ENSEMBL", mirror = "www")
#mart = useMart(biomart = "ensembl", dataset = opt$dataset, host = "www.ensembl.org")

bmIDs = fread(opt$trans)
bmIDs_g = fread(opt$gene)

#bmIDs = getBM(attributes=c('ensembl_transcript_id','description'),mart = mart)
#bmIDs_g = getBM(attributes = c('ensembl_gene_id', 'description'), mart = mart)

result_path <- paste0(getwd(), "/", opt$outdir)
wpath <- opt$indir
src_dir <- paste0(wpath, "stringtie/")
setwd(src_dir)
wdir <- getwd()

src_files <- list.files(getwd(), pattern = "transcript.gtf", recursive = T, include.dirs = T)
src_files_cnt <- length(src_files)
#src_files <- src_files[1:src_files_cnt]
#src_files_cnt <- length(src_files)

for(i in 1:src_files_cnt){
  print(paste0("working ", i, " in ", src_files_cnt))
  # wdir <- paste0(src_dir, src_files[i])
  # setwd(paste0(wdir, src_files[i]))
  print(src_files[i])
  transcript <- fread(src_files[i], sep = "\t", header = F, stringsAsFactors = T)

  transcript <- transcript[, -c(1, 2, 6, 7, 8)]

  names(transcript)[1] <- c("Name")
  names(transcript)[2] <- c("Start")
  names(transcript)[3] <- c("End")
  transcript <- subset(transcript, Name=="transcript")

  trans_sep <- separate(transcript, V9, c("ensembl_gene_id", "ensembl_transcript_id", "ref_gene_name", "cov", "FPKM", "TPM"), sep = "\\;")
  trans_sep <- trans_sep[, -c(7, 9)]
  trans_sep$ensembl_gene_id <- gsub("gene_id ", "", trans_sep$ensembl_gene_id)
  trans_sep$ensembl_gene_id <- gsub("\"", "", trans_sep$ensembl_gene_id)
  trans_sep$ensembl_transcript_id <- gsub(" transcript_id ", "", trans_sep$ensembl_transcript_id)
  trans_sep$ensembl_transcript_id <- gsub("\"", "", trans_sep$ensembl_transcript_id)
  trans_sep$ref_gene_name <- gsub(" ref_gene_name ", "", trans_sep$ref_gene_name)
  trans_sep$ref_gene_name <- gsub(" gene_name ", "", trans_sep$ref_gene_name)
  trans_sep$ref_gene_name <- gsub("\"", "", trans_sep$ref_gene_name)
  trans_sep$FPKM <- gsub("\"", "", trans_sep$FPKM)
  trans_sep$FPKM <- gsub(" FPKM ", "", trans_sep$FPKM)
  mode(trans_sep$FPKM) <- "double"
  names(trans_sep)[7] <- paste0(src_files[i], "_FPKM")

  if(i ==1) {
    sprintf("It is the %s out of %s.", i, src_files_cnt)
    trans_sep_final <- trans_sep
  }else if(1 < i && i < src_files_cnt) {
    sprintf("It is the %s out of %s.", i, src_files_cnt)
    trans_sep <- trans_sep[, c(5, 7)]
    trans_sep_final <- merge(x = trans_sep_final, y = trans_sep, by = "ensembl_transcript_id", all = FALSE)
  }else {
    sprintf("It is the %s out of %s.", i, src_files_cnt)
    trans_sep <- trans_sep[, c(5, 7)]
    trans_sep_final <- merge(x = trans_sep_final, y = trans_sep, by = "ensembl_transcript_id", all = FALSE)
    trans_sep_final <- merge(x = trans_sep_final, y = bmIDs, by = "ensembl_transcript_id", all = FALSE)
  }
  print(paste0(i, "done"))
}
setDF(trans_sep_final)
ttemp <- trans_sep_final[,7:(ncol(trans_sep_final)-1)]
trans_sep_final["avg"] <- apply(ttemp, 1, mean)
trans_sep_final <- subset(trans_sep_final, avg!=0)
trans_sep_final <- trans_sep_final[,-c(ncol(trans_sep_final))]
setwd(result_path)
write.csv(trans_sep_final, "transcript_FPKM.csv", row.names = FALSE)

trans <- trans_sep_final %>% dplyr::select(5:(ncol(trans_sep_final)-1))
trans <- trans %>% group_by(ensembl_gene_id, ref_gene_name) %>% summarise_all(funs(sum))

g_df_final <- merge(x = trans, y = bmIDs_g, by = "ensembl_gene_id", all = FALSE)
setwd(result_path)
write.csv(g_df_final, "gene_FPKM.csv", row.names = FALSE)

rm()
quit(save="no")
