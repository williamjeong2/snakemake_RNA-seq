#!/usr/bin/env Rscript
list.of.packages <- c("data.table", "tidyr", "dplyr", "optparse")
list.of.bio.packages <- c("")
not.installed.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(not.installed.packages)) install.packages(not.installed.packages, repos="http://cran.rstudio.com/", dependencies=TRUE)

if(length(not.installed.bio.packages)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(not.installed.bio.packages, suppressUpdates = TRUE)
}
lapply(list.of.bio.packages, require, character.only = T
lapply(list.of.packages, require, character.only = TRUE)

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

for(i in 1:src_files_cnt){
  print(paste0("working ", i, " in ", src_files_cnt))
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
colnames(trans_sep_final) <- gsub("/transcript.gtf", "", colnames(trans_sep_final), fixed = T)
write.csv(trans_sep_final, "transcript_FPKM.csv", row.names = FALSE)

trans <- trans_sep_final %>% dplyr::select(5:(ncol(trans_sep_final)-1))
trans <- trans %>% group_by(ensembl_gene_id, ref_gene_name) %>% summarise_all(funs(sum))

g_df_final <- merge(x = trans, y = bmIDs_g, by = "ensembl_gene_id", all = FALSE)
colnames(g_df_final) <- gsub("/transcript.gtf", "", colnames(g_df_final), fixed = T)
setwd(result_path)
write.csv(g_df_final, "gene_FPKM.csv", row.names = FALSE)

rm()
quit(save="no")
