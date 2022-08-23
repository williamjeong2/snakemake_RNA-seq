list.of.packages <- c("optparse", "readxl", "data.table", "ggplot2", "tidyr", "dplyr", "tibble")
list.of.bio.packages <- c("DESeq2", "DEGreport", "clusterProfiler", "DOSE", "org.Hs.eg.db", "pheatmap", "org.Mm.eg.db", "umap",
                          "genefilter", "RColorBrewer", "GO.db", "topGO", "gage", "ggsci", "curl", "biomaRt", "fgsea", "EnhancedVolcano")

# Install Require packages
not.installed.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(not.installed.packages)) install.packages(not.installed.packages, repos="http://cran.rstudio.com/", dependencies=TRUE, Ncpus = 6)

not.installed.bio.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[, "Package"])]
if(length(not.installed.bio.packages)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(not.installed.bio.packages, update = TRUE, Ncpus = 6, suppressUpdates = TRUE)
}
lapply(list.of.bio.packages, require, character.only = TRUE)
lapply(list.of.packages, require, character.only = TRUE)

# arguments to provide
option_list = list(
  # default params
  make_option("--count", type="character", default="results/counts.txt", metavar="character"), # read count data 
  make_option("--metadata", type="character", default="data/sample-list.xlsx", metavar="character"), # "samples.tsv"
  make_option("--outdir", type="character", default="results/visualization/", metavar="character"), # dir to save
  # PCA params
  make_option("--pcax", type="character", default="1", metavar="character"), # string
  make_option("--pcay", type="character", default="2", metavar="character"), # string
  # heatmap parmas
  make_option("--fdrval", type="double", default=0.05, metavar="number"),
  make_option("--ntopgene", type="integer", default= 30, metavar="number"),
  make_option("--hmapcolor", type="character", default="YlOrBr", metavar="character"), # RdBu
  # GSEA
  make_option("--gseafdr", type="double", default=0.25, metavar="number"),
  make_option("--gseapval", type="double", default=0.01, metavar="number")
  # make_option("--gmtpath", type="character", default="", metavar="character"),
)

# parse the command-line arguments and pass them to a list called 'opt'
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

detectGroups <- function (x){  # x are col names
  tem <- gsub("[0-9]*$","",x) # Remove all numbers from end
  tem <- gsub("_$","",tem); # remove "_" from end
  return( tem )
}

savePlot <- function(path, plot, width = 7, height = 7){
  for(i in c("svg")) {
    ggsave(filename = paste0(path, ".", i),
           plot = plot,
           device = i, scale = 2,
           width = width, height = height, units = "in",
           dpi = 320)
  }
}

####################
## Prepare data
####################

# Import counts table
countdata <- read.table(opt$count, header = TRUE, row.names = 1, sep = "\t", check.names = F)

# Remove .bam from column identifiers
colnames(countdata) <- gsub(".sorted.bam", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("temp.mapped.", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("_REP$", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("_rep$", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("_Rep$", "", colnames(countdata), fixed = T)

if (grepl('^ENSG', row.names(countdata)[1])) {
  organism <- org.Hs.eg.db
}
if (grepl('^ENSMUSG', row.names(countdata)[1])) {
  organism <- org.Mm.eg.db
}

# Remove chr, start, end, strand, length columns
countdata <- countdata[, c(-1)]

comparison <- read_excel(opt$metadata, skip = 2, sheet = "comparisons")
comp = as.data.frame(paste(paste(comparison$ctrl_cellType,comparison$ctrl_group, sep = "__"), paste(comparison$treat_cellType,comparison$treat_group, sep = "__"), sep = "-"))
colnames(comp) <- "comparison_list"

samples <- read_excel(opt$metadata, sheet = "samples")
colnames(samples)[1] <- "samples"
samples$samples = gsub(".gz", "", samples$samples)
samples$samples = gsub(".fq", "", samples$samples)
samples$samples = gsub(".fastq", "", samples$samples)
rownames(samples) = samples$samples

for(i in 1:nrow(comp)){
  ifelse(!dir.exists(paste0(opt$outdir, comp[i, ])), dir.create(paste0(opt$outdir, comp[i, ]), showWarnings = FALSE, recursive = T), print("Directory already created"))
}

metad<- samples[colnames(countdata), ]
metad<- metad[, 'cell_type']

dds <- DESeqDataSetFromMatrix(
  countData = countdata[,1:ncol(countdata)],
  colData = samples[colnames(countdata), ],
  design = ~ treatment + cell_type
)
dds <- DESeq(dds, parallel = TRUE)

dds_norm <- vst(dds)

## PCA
pcaData <- plotPCA(dds_norm, intgroup=c("treatment", "cell_type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
shapes = c(16,17,15,3,7,8,   1,2,4:6,9:15,18:25)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), aspect.ratio=1)
p <- ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=cell_type, label=row.names(pcaData))) + scale_shape_manual(values= shapes) + geom_text(aes(colour = treatment), vjust = -0.5) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("Principal component analysis (PCA)") + 
  coord_fixed()  + theme
savePlot(paste0(opt$outdir, "/PCA_all"), p)

## UMAP
normalized_counts <- assay(dds_norm) %>% t()  # We need to transpose this data so each row is a sample

umap_results <- umap::umap(normalized_counts) # Now perform UMAP on the normalized data
umap_plot_df <- data.frame(umap_results$layout) %>% 
  tibble::rownames_to_column("samples") %>%
  dplyr::inner_join(samples[colnames(countdata),] %>% select(c("samples", "treatment", "cell_type")), by = "samples")
rownames(umap_plot_df) = umap_plot_df$samples
p <- ggplot(umap_plot_df, aes(x = X1, y = X2, color = treatment, shape = cell_type, label=row.names(umap_plot_df))) + scale_shape_manual(values= shapes) + geom_text(aes(colour = treatment), vjust = -0.5) +
  geom_point(size=3) +
  xlab(paste0("UMAP Dimension 1")) +
  ylab(paste0("UMAP Dimension 2")) +
  ggtitle("UMAP of All Samples") + 
  coord_fixed()  + theme
savePlot(paste0(opt$outdir, "/UMAP_all"), p)



for(i in 1:nrow(comp)){
  # get comparison value
  selectComparisons = strsplit(comp[i, ], "-")
  comp1 = selectComparisons[[1]][1]
  comp2 = selectComparisons[[1]][2]
  comp1_tb = countdata[, subset(samples, paste(cell_type, treatment, sep="__") == comp1)$samples]
  comp2_tb = countdata[, subset(samples, paste(cell_type, treatment, sep="__") == comp2)$samples]
  comp_tb = data.frame(c(comp1_tb, comp2_tb), row.names = rownames(countdata), check.names = F)
  rm(comp1_tb, comp2_tb)
  
  metadata = samples[colnames(comp_tb), ]
  metadata = unite(metadata, condition, c(cell_type, treatment), remove = T)
  metadata$condition <- factor(metadata$condition)
  metadata <- metadata[,"condition"]
  
  # make DESeq2 object from counts and metadata
  # - countData : count dataframe
  # - colData : sample metadata in the dataframe with row names as sampleID's
  # - desigh : The design of the comparisons to use.
  #            Use (~) befor the name of the column variable to compare
  ddsMat <- DESeqDataSetFromMatrix(countData = comp_tb,
                                   colData = metadata,
                                   design = ~condition)
  
  
  # Find differential expressed genes
  library("BiocParallel")
  # register(MulticoreParam(((detectCores(all.tests = FALSE, logical = TRUE))*2)%/%3))
  ddsMat <- DESeq(ddsMat, parallel = TRUE)
  
  ####################
  ## PCA Plot
  ####################
  rld <- vst(ddsMat)
  rld_mat <- assay(rld)
  pca <- prcomp(t(rld_mat))
  PCAxy <- c(as.integer(opt$pcax),as.integer(opt$pcay)) # selected principal components
  percentVar=round(100*summary(pca)$importance[2, PCAxy],0)
  
  # meta <- data.frame(detectGroups(colnames(comp_tb)), row.names = colnames(comp_tb))
  meta = data.frame(metadata, row.names = colnames(comp_tb))
  colnames(meta)[1] <- "group"
  df <- cbind(meta, pca$x[, PCAxy])
  
  if(ncol(comp_tb)<20){ # change size depending of # samples
    p = ggplot(df, aes(x=PC1, y=PC2+.5, color=group, label=row.names(df)), size=5) + geom_point(aes(x=PC1, y= PC2, color = group), size=5) + geom_text(aes(colour = group), vjust = -0.5)
  }else if(ncol(comp_tb)<50){
    p = ggplot(df, aes(x=PC1, y=PC2+.5, color=group, label=row.names(df)), size=3) + geom_point(aes(x=PC1, y= PC2, color = group), size=3) + geom_text(aes(colour = group), vjust = -0.5)
  }else{
    p = ggplot(df, aes(x=PC1, y=PC2+.5, color=group, label=row.names(df)), size=2) + geom_point(aes(x=PC1, y= PC2, color = group), size=2) + geom_text(aes(colour = group), vjust = -0.5)
  }
  
  # for showing shapes on ggplot2. The first 6 are default. Default mapping can only show 6 types.
  shapes = c(16,17,15,3,7,8,   1,2,4:6,9:15,18:25)
  theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), aspect.ratio=1)
  pca.group <- p + scale_shape_manual(values= shapes) + 
    xlab(paste0("PC", opt$pcax ,": ", percentVar[1],"% variance")) +
    ylab(paste0("PC", opt$pcay ,": ", percentVar[2],"% variance")) + 
    ggtitle("Principal component analysis (PCA)") + 
    coord_fixed(ratio=1.0) + theme
  
  savePlot(paste0(opt$outdir, comp[i, ], "/PCA_by_condition"), pca.group)
  
  ####################
  ## Heatmap Plot
  ####################
  
  # Get results from testing with FDR adjust pvalues
  # results <- results(ddsMat, pAdjustMethod = "fdr", alpha = opt$fdrval)
  results <- results(ddsMat)

  # Add ENSEMBL
  results$ensembl <- mapIds(x = organism,
                            keys = row.names(results),
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")
  
  # Subset for only significant genes( q < 0.05)
  results_sig <- subset(results, pvalue < opt$fdrval &
                          (log2FoldChange < -0.7 | log2FoldChange > 0.7))
  comp_tb$row <- rownames(comp_tb)
  results_sig$row <- rownames(results_sig)
  sig_tb = merge(x = as.data.frame(results_sig), y = comp_tb, by = "row", all=F)
  colnames(sig_tb)[1] <- "geneid"
  fwrite(x = sig_tb,
         file = paste0(opt$outdir, comp[i, ], "/results_gene_annotated_significant.tsv"),
         sep = '\t', row.names = F)

  # 
  # mat_t <- biomaRt::select(organism, keys=row.names(mat), columns = c("SYMBOL", "ENSEMBL"), keytype = "ENSEMBL")
  # mat_t$SYMBOL[is.na(mat_t$SYMBOL)] <- mat_t$ENSEMBL[is.na(mat_t$SYMBOL)]
  # mat_t = subset(mat_t, !duplicated(subset(mat_t, select=c("ENSEMBL"))))
  # rownames(mat_t) <- mat_t$ENSEMBL
  # rownames(mat) <- mat_t$SYMBOL
  
  # Make Heatmap with pheatmap function
  pdf(NULL)
  temp <- vst(ddsMat)
  resSig <- subset(results(ddsMat), padj < 0.05 & 
                     (log2FoldChange < -1 | log2FoldChange > 1))
  if (nrow(as.data.frame(resSig)) != 0) {
    if (length(rownames(resSig[order(resSig$padj),])) > opt$ntopgene){
      select <- rownames(resSig[order(resSig$padj),])[1:opt$ntopgene]
    }else {
      select <- rownames(resSig[order(resSig$padj),])
    }
    if (opt$ntopgene > 200){
      show_rownames <- FALSE
    }else {
      show_rownames <- TRUE
    }
    if (dim(assay(temp))[2] > 200){
      show_colnames <- FALSE
    }else {
      show_colnames <- TRUE
    }
    heatmap.plot <- pheatmap(assay(temp)[select, ], 
            cluster_rows=T,
            cluster_cols=T,
            color = colorRampPalette(rev(brewer.pal(n = 9, name = opt$hmapcolor)))(255),
            breaks = seq(-2, 2, length.out = 255),
            show_rownames = show_rownames,
            show_colnames = show_colnames,
            scale="row",
            border_color = NA,
            angle_col = 315,
            fontsize_col = 5, 
            fontsize_row = 6)
    savePlot(paste0(opt$outdir, comp[i, ], "/heatmap"), heatmap.plot)

  }
  
  
  ####################
  ## Volcano Plot
  ####################
  
  # Gather Log-Fild change and FDR-corrected pvalues from DESeq2 results
  data <- data.frame(gene = row.names(results),
                     padj = results$padj,
                     lfc = results$log2FoldChange,
                     symbol = results$ensembl)
  
  # Replace any emsembl to geneId that have NA 
  data$ensembl[is.na(data$ensembl)] <- data$gene[is.na(data$ensembl)]
  data <- na.omit(data)
  
  # Color the points which are up or down
  ## If fold-change > 2 and pvalue > 0.05 (Increased significant)
  ## If fold-change < -2 and pvalue > 0.05 (Decreased significant)
  data <- data %>%
    mutate(color = case_when(data$lfc > 1 & data$padj < 0.05 ~ "Increased",
                             data$lfc < -1 & data$padj < 0.05 ~ "Decreased",
                             data$lfc > -1 | data$lfc < 1 | data$padj > 0.05 ~ "No Significant"))
  
  data <- data %>% mutate(genelabels = "")
  data <- data %>% arrange(padj)
  
  data$genelabels[1:30] <- as.character(data$symbol[1:30])
  
  # # Make a basic ggplot2 object with x-y values
  # vol <- ggplot(data, aes(x = lfc, y = pval, color = color))
  
  # # Add ggplot2 layers
  # vol + 
  #   ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  #   geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  #   scale_color_manual(name = "Directionality",
  #                     values = c(Increased = "#008B00", Decreased = "#CD4F39", nosignificant = "darkgrey")) +
  #   theme_bw(base_size = 14) +
  #   theme(legend.position = "right") +
  #   xlab(expression(log[2]("treat"/"wt"))) +
  #   ylab(expression(-log[10]("adjusted p-value"))) +
  #   geom_hline(yintercept = 1.3, colour = "darkgrey") +
  #   scale_y_continuous(trans = "log1p")
  # dev.off()
  
  
  
  volcano.plot <- EnhancedVolcano(results,
                                  lab = results$ensembl,
                                  subtitle = "Differential expression",
                                  x = 'log2FoldChange',
                                  y = 'pvalue',
                                  selectLab = data$genelabels[order(data$padj, decreasing = F)[1:10]],
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  pCutoff = 10e-14,
                                  FCcutoff = 2.0,
                                  pointSize = 2.5,
                                  labSize = 6.0,
                                  labCol = 'black',
                                  labFace = 'bold',
                                  boxedLabels = TRUE,
                                  colAlpha = 4/5,
                                  legendPosition = 'right',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 1.0,
                                  # col = c('grey', '#56abec', '#355C7D', '#FD6666'),
                                  col = c('grey', 'grey', 'grey', '#FD6666'),
                                  gridlines.major = FALSE,
                                  gridlines.minor = FALSE)
  
  savePlot(paste0(opt$outdir, comp[i, ], "/volcano"), volcano.plot)
  
  ####################
  ## GSEA
  ####################
  
  results$ens <- row.names(results)
  ens2symbol <- AnnotationDbi::select(organism,
                                      keys = results$ens,
                                      columns = "SYMBOL",
                                      keytype = "ENSEMBL")
  ens2symbol <- as_tibble(ens2symbol)
  results_sig = as.data.frame(results_sig)
  results_sig$row = row.names(results_sig)
  
  res = inner_join(results_sig, ens2symbol, by=c("row"="ENSEMBL"))
  res2 <- res %>%
    dplyr::select(SYMBOL, stat) %>%
    na.omit() %>% 
    distinct() %>%
    group_by(SYMBOL) %>%
    summarize(stat=mean(stat))
  
  ranks <- deframe(res2)
  
  pathways.gobp <- gmtPathways("scripts/MSigDB/c5.go.bp.v7.5.1.symbols.gmt")
  pathways.reactome <- gmtPathways("scripts/MSigDB/ReactomePathways.gmt")
  # pathways.reactome <- gmtPathways("scripts/MSigDB/c2.cp.reactome.v7.5.1.symbols.gmt")
  
  fgseaRes.gobp <- fgsea(pathways=pathways.gobp, stats=ranks, min = 15, max = 500)
  fgseaRes.reactome <- fgsea(pathways=pathways.reactome, stats=ranks, min = 15, max = 500)
  
  # GOBP
  fgseaResTidy.gobp <- fgseaRes.gobp %>%
    as_tibble() %>%
    arrange(desc(abs(NES))) 
  
  fgseaResTidy.gobp$leadingEdge <- vapply(fgseaResTidy.gobp$leadingEdge, paste, collapse = ", ", character(1L))
  fwrite(fgseaResTidy.gobp[fgseaResTidy.gobp$ES > 0 & fgseaResTidy.gobp$padj < opt$gseafdr & fgseaResTidy.gobp$pval < opt$gseapval, ], file = paste0(opt$outdir, comp[i, ], "/UP_", comp[i, ], "_GO_pathway.tsv"), sep = "\t")
  fwrite(fgseaResTidy.gobp[fgseaResTidy.gobp$ES < 0 & fgseaResTidy.gobp$padj < opt$gseafdr & fgseaResTidy.gobp$pval < opt$gseapval, ], file = paste0(opt$outdir, comp[i, ], "/DOWN_", comp[i, ], "_GO_pathway.tsv"), sep = "\t")
  
  gsea.subtitle = paste0("FDR < ", as.character(opt$gseafdr), " & p-value < ", as.character(opt$gseapval), " & abs(NES) > 2")
  
  gobp.fig1 <- ggplot(fgseaResTidy.gobp[fgseaResTidy.gobp$padj < opt$gseafdr & fgseaResTidy.gobp$pval < opt$gseapval & abs(fgseaResTidy.gobp$NES) > 2, ][1:10, ], aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill= padj < opt$gseafdr & pval < opt$gseapval & abs(NES) > 2)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title = "GO pathways NES from GSEA",
         subtitle = "Top 10 abs(NES)",
         fill = gsea.subtitle) +
    theme_minimal()
  
  savePlot(paste0(opt$outdir, comp[i, ], "/GO_pathways_", comp[i, ], "_NES_from_GSEA"), gobp.fig1, 7, 4)
  
  # Reactome
  fgseaResTidy.reactome <- fgseaRes.reactome %>%
    as_tibble() %>%
    arrange(desc(abs(NES))) 
  
  fgseaResTidy.reactome$leadingEdge <- vapply(fgseaResTidy.reactome$leadingEdge, paste, collapse = ", ", character(1L))
  fwrite(fgseaResTidy.reactome[fgseaResTidy.reactome$ES > 0 & fgseaResTidy.reactome$padj < opt$gseafdr & fgseaResTidy.reactome$pval < opt$gseapval, ], file = paste0(opt$outdir, comp[i, ], "/UP_", comp[i, ], "_Reactome_pathway.tsv"), sep = "\t")
  fwrite(fgseaResTidy.reactome[fgseaResTidy.reactome$ES < 0 & fgseaResTidy.reactome$padj < opt$gseafdr & fgseaResTidy.reactome$pval < opt$gseapval, ], file = paste0(opt$outdir, comp[i, ], "/DOWN_", comp[i, ], "_Reactome_pathway.tsv"), sep = "\t")
  
  reactome.fig1 <- ggplot(fgseaResTidy.reactome[fgseaResTidy.reactome$padj < opt$gseafdr & fgseaResTidy.reactome$pval < opt$gseapval & abs(fgseaResTidy.reactome$NES) > 2, ][1:10, ], aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill= padj < opt$gseafdr & pval < opt$gseapval & abs(NES) > 2)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title = "Reactome pathways NES from GSEA",
         subtitle = "Top 10 abs(NES)",
         fill = gsea.subtitle) +
    theme_minimal()
  
  savePlot(paste0(opt$outdir, comp[i, ], "/Reactome_pathways_", comp[i, ], "_NES_from_GSEA"), reactome.fig1, 7, 4)
}
