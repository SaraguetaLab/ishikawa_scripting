#!/usr/bin/Rscript
########################################################
## Packages
library(tidyverse)
library(DESeq2)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(org.Hs.eg.db)
library(vsn)
library(genefilter)
library(ggrepel)
library(ggpubr)
library(ggExtra)
library(cowplot)
library(gplots)
library(knitr)
library(magrittr)
library(BiocParallel)
library(forcats)
########################################################

## Custom functions

fileNaming <- function(file_prefix, var, var1, ext, dir_prefix) {
  if (missing(dir_prefix)) {
    name = paste0(file_prefix, "_", var,"_", var1, ".", ext)
  } else {
    name = paste0(dir_prefix, "_", var, "/", file_prefix, "_", var,"_", var1,".", ext)
  }
  return(name)
}

go_deseq <- function(input_counts, input_metadata, job_title, feature, 
                      cond1, cond2) {
  ## prepare data
  raw_counts <- read.csv(input_counts, header = T)
  rownames(raw_counts) <- raw_counts$geneName
  raw_counts$geneName <- NULL
  raw_counts[is.na(raw_counts)] <- 0
  raw_counts<-raw_counts[ , order(names(raw_counts))]

  ## prepare conditions
  conditions <- read.csv(input_metadata, header = T, na.strings=c("","NA"),
    stringsAsFactors = F) %>% 
    dplyr::rename(filename=file, condition = feature) %>%
    tidyr::drop_na(condition) %>%
    dplyr::mutate(condition=as.factor(condition),type=as.factor(type))
  rownames(conditions) <- conditions$filename
  

  ## check data vs conditions
  if (all(colnames(raw_counts) == rownames(conditions)) == FALSE) {
    id_vec <- intersect(as.vector(conditions$filename),as.vector(colnames(raw_counts)))
    id_vec = id_vec[!(id_vec %in% id_remove)]
    raw_counts<- raw_counts[ names(raw_counts)[names(raw_counts) %in% id_vec] ]
    conditions <- subset(conditions, filename %in% id_vec)
    # check again
    test <- all(colnames(raw_counts) == rownames(conditions))
    if (test == TRUE) {
      print('all good!')
    } else {
      print('something is wrong with sample labels...check manually!')
    }
  }
  
  ## perform core DESeq2 process
  ddsMatrix <- DESeqDataSetFromMatrix(countData = raw_counts,
                                      colData = conditions,
                                      design = ~ condition)
  ddsMatrix <- ddsMatrix[ rowSums(counts(ddsMatrix)) > 100, ]
  ddsMatrix

  dds <- DESeq(ddsMatrix,parallel = TRUE)
  
  ## extract DE results
  resdf <- results(dds, tidy = TRUE, contrast = c("condition", cond1, cond2)) %>%
    dplyr::filter(pvalue < 0.05) %>%
    dplyr::filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
    dplyr::arrange(desc(log2FoldChange)) %>%
    tibble::as_tibble()
  
  ## export DE list for "feature"
  name_de = fileNaming("DEresults_deseq2",feature, job_title,"csv")
  write.csv(resdf, file = name_de, row.names = F, quote = F)
  
  
  ## Transformation of normalized counts
  vsd <- vst(dds, blind = FALSE)
  df_vsd <- as.data.frame(assay(vsd))
  df_vsd$gene_name <- rownames(df_vsd)
  name_dfvsd = fileNaming("count_matrix_vstransformation", feature, job_title,"csv")
  write.csv(df_vsd, name_dfvsd, row.names = F)

  ## pax genes profile
  df_vsd_m <- melt(df_vsd)
  paxGenes <- df_vsd_m %>% dplyr::filter(grepl("^PAX[0-9]*$", gene_name)) %>% mutate(name=as.character((variable)))
  pax <- merge(paxGenes, clin_tcga_filt[, c("name", "stage")], by="name")

  ggplot(pax, aes(x=stage,y=value,fill=gene_name)) + geom_boxplot()+
  #scale_fill_manual(values = c("lightgoldenrod4","coral","dodgerblue","gray40"))+
  facet_wrap(~ gene_name, ncol=3, scales = "free_y")+ 
  labs(x="", y="Normalized gene expression")+
  theme(legend.title = element_blank(), legend.position = "none",
        axis.text.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(angle = 60, hjust = 1,vjust = 1))



  # plot PCA (DESeq2 function)
  name_pca = fileNaming("pca_vsd_normcounts", feature, job_title,"pdf")
  pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE,ntop=500)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) + geom_vline(xintercept = 0, linetype = "longdash") +
    geom_hline(yintercept = 0, linetype = "longdash") +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()
  
  ggsave(name_pca, width = 6, height = 6)
  
  ## plot components in desity graphs
  theme_set(theme_bw(base_size=20) + theme(strip.background = element_blank(),panel.grid =element_blank()))    
  name_dens = fileNaming("dens_pc1-pc2_vsd_normcounts", feature, job_title,"pdf")
  pc1.data <- as.data.frame(pcaData$PC1)
  rownames(pc1.data) <- rownames(pcaData$name)
  pc1.data$stage <- clin_tcga_filt$stage
  pc1.data$type <- clin_tcga_filt$type
  pc1.data <- pc1.data %>% dplyr::select(value="pcaData$PC1", stage, type)
  pc1.data$PC <- "PC1"


  pc2.data <- as.data.frame(pcaData$PC2)
  rownames(pc2.data) <- rownames(pcaData$name)
  pc2.data$stage <- clin_tcga_filt$stage
  pc2.data$type <- clin_tcga_filt$type
  pc2.data <- pc2.data %>% dplyr::select(value="pcaData$PC2", stage,type)
  pc2.data$PC <- "PC2"

  pcs <- rbind(pc1.data,pc2.data)

  ggplot(pcs, aes(x=value, color=stage, group=stage)) + geom_density() + facet_wrap(~PC, scales = "free_y") +
  theme(legend.position = "bottom", legend.title = element_blank()) + geom_vline(xintercept = 0, linetype = "longdash")
  
  ggsave(name_dens, width = 6, height = 6)

  
  # annotation data for heatmaps
  annotation_col <- data.frame(name=conditions$filename,
                               condition = conditions$condition,
                               type=conditions$type)
  annotation_col$name <- as.character(annotation_col$name)
  
  annotation_col <- annotation_col %>% arrange(condition)
   
  annot <- as.vector(annotation_col$name)
  rownames(annotation_col) <- annotation_col$name
  annotation_col$name <- NULL
  

  mat_complete <- assay(vsd)
  
  ## plot distance matrix
  # Distance matrix (with rld norm data)
  sampleDists <- dist(t(mat_complete))
  sampleDistMatrix <- as.matrix( sampleDists )
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  name_dist = fileNaming("distanceMatrix",feature,job_title,"pdf")
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors, annotation_row = annotation_col,
           filename = name_dist,
           cellwidth = 20,
           cellheight = 20)
  
  
  ## plot DE genes for feature of interest
  # filter by most DE genes (if different to filter in go_deseq function)
  
  
  # extract genes DE from deseq results table
  de_genes <- as.vector(res_top$row)
  
  
  
  
  hmcols <- colorRampPalette(c("blue","gray90", "darkred"))(50)
  mat <- mat_complete[row.names(mat_complete)%in%de_genes,]
  mat  <- mat - rowMeans(mat)
  mat <- mat[ de_genes,annot]
  
  name_topde = fileNaming("DEexpressingGenes",feature,job_title,"pdf")
  pheatmap(mat, col = hmcols,
           main = "",
           cluster_cols = F,
           border_color = NA,
           cluster_rows = F,
           show_rownames = F, show_colnames = F,
           gaps_col = 0,
           cellwidth = 1.5,
           cellheight = 1, annotation_col = annotation_col,
           scale = "row",
           fontsize = 4,
           filename = name_topde
  )
  
    list_of_results <- list("dds" = dds, "resdf" = resdf,
                            "vsd" = vsd,"norms"=df,"conditions" = conditions)

  return(list_of_results)
}

########################################################
## Set theme for plots
theme_set(theme_bw(base_size = 16) +
            theme(panel.background = element_blank(),
                  panel.grid = element_blank(), panel.border = element_blank()))
########################################################


########################################################
## run deseq2 with interesting feature
########################################################
## Input arguments
args = commandArgs(trailingOnly = TRUE)

# Path to raw counts file
path_to_counts <- args[1]

# Path to where metadata is stored.
path_to_metadata <- args[2]

# Job description ("BBK", "CCT", "BBK-CCT"....)
job_title <- args[3]

# Set feature of interest
feature <- args[4]
condition1 <- args[5]
condition2 <- args[6]

# number of cores for parallel processing
cores <- 6
register(MulticoreParam(cores))


# Run go_deseq to perform deseq and optional normalization methods (vsd or rld)
resDeseq <- go_deseq(path_to_counts, path_to_metadata, job_title, feature,
                      condition1, condition2)
