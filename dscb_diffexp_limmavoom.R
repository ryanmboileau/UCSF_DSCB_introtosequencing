# install.packages("dplyr")
# install.packages("ggplot")
# if (!requireNamespace("BiocManager", quietly = TRUE)) +
#   install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("limma","edgeR"))
# 
# install.packages("data.table")
# install.packages("ggplot")
# install.packages("dplyr")

library(limma)
library(edgeR)
library(ggplot2)
library(dplyr)
library(data.table)



# -----------------------------------------------------------
# -----------------Pre set code for Limma-Voom -------------------

#setwd("~/Desktop/2021.DSCB257.discussion/20201014_DSCB_discussion_download")

# Read gene id and gene symbol files
# Replace with download location on your computer
#geneid_symbol <- read.csv('~/Github/Blellochlab/DSCB257/download/geneid_symbol.csv')
geneid_symbol <- read.csv('geneid_symbol.csv')
geneLength <- read.csv(file = "geneLength.csv", 
                       row.names = 1)
# Load counts data
# Replace with download location on your computer
counts <- read.delim(file = "counts_discussion.txt", 
                     sep = "\t", header = T)
dim(counts)

# Set up the sample table
sampleTable <- data.frame(row.names = colnames(counts),
                          stage = as.factor(c('E3.5_ICM', 'E3.5_ICM', 'E3.5_ICM',
                                              'E4.5_Epi', 'E4.5_Epi', 'E4.5_Epi',
                                              'E5.5_Epi', 'E5.5_Epi', 'E5.5_Epi')))
#Merge counts and sampleTable using the DGEList function from edgeR
x <- DGEList(counts=counts,samples=sampleTable)

#Now we normalize the data using TMM (again using edgeR)
x_norm <- calcNormFactors(x, method = "TMM")
#We can look at the normalization factors here
x_norm$samples$norm.factors

#Log transform the normalized data 
#Convert to rpkm function from edgeR. Also, log transform to generate lcpm.
rpkm <- rpkm(x, gene.length = geneLength$Length, log = T)
rpkm_norm <- rpkm(x_norm, gene.length = geneLength$Length, log = T)

#Vizualize the effect of normalization 
par(mfrow = c(1,2))
boxplot(rpkm, las=2,  main="")
title(main="Unnormalised data",ylab="rpkm")

boxplot(rpkm_norm, las=2,  main="")
title(main="TMM Normalised data",ylab="rpkm")


## Plot PCA (Vizualize the dimensions of this data that contain the most information)
##Center the normalized rpkm data
rpkm_norm_centered <- as.data.frame(rpkm_norm)
rpkm_norm_centered$average <- rowMeans(rpkm_norm_centered, dims = 1)
rpkm_norm_centered <- rpkm_norm_centered %>% mutate(across(1:9) - average)
rpkm_norm_centered <- rpkm_norm_centered[1:9]
#run PCA analysis using prcomp and plot PC1 variance and screeplot
rpkm_pca <- prcomp(rpkm_norm_centered, scale. = TRUE)
rpkm_pca_df <- as.data.frame(rpkm_pca$x)

#calculate percentage of variance explained by each PC
percentage <- round(rpkm_pca$sdev / sum(rpkm_pca$sdev) * 100, 2)
percentage <- paste(colnames(rpkm_pca_df), "(", paste(as.character(percentage), "%", ")", sep=""))
percentage
#use output of prcomp to generate dataframe and add construct object/names of samples to the dataframe as samplename
rotationdf <- as.data.frame(rpkm_pca$rotation)
rotationdf$Stage <- sampleTable$stage

#Plot the data using the coordinates for the first two Principal Components
ggplot(as.data.frame(rotationdf), aes(x=PC1, y=PC2, color = Stage)) +
  geom_point() +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  scale_color_manual(breaks = c('E3.5_ICM', 'E4.5_Epi', 'E5.5_Epi'),
                     values = c("royalblue", "green", "orange"))

#Plot the data using the coordinates for the third and fourth Principal Components
ggplot(as.data.frame(rotationdf), aes(x=PC3, y=PC4, color = Stage)) +
  geom_point() +
  xlab(percentage[3]) +
  ylab(percentage[4]) +
  scale_color_manual(breaks = c('E3.5_ICM', 'E4.5_Epi', 'E5.5_Epi'),
                     values = c("royalblue", "green", "orange"))

##Principal Components have certain limitations. Plot the data using UMAP, which collapses many dimensions of data
# into 2D
BiocManager::install("umap")
library(umap)
###Transposing rpkm_norm_centered data
rpkm_norm_centered_t <- (t(rpkm_norm_centered))

#run UMAP package on df_lcpm, batchcorrected if necessary
rpkm_umap <- umap(rpkm_norm_centered_t, n_neighbors = 4)
head(rpkm_umap$layout)
rpkm_umap_df <- as.data.frame(rpkm_umap$layout)
colnames(rpkm_umap_df) <- c("one", "two")
rpkm_umap_df$Stage <- sampleTable$stage
ggplot(rpkm_umap_df, mapping = aes(x=one , y=two,color = Stage)) +
  geom_point() +
  xlab("UMAP 1") +
  ylab("UMAP 2")

## Run Limma-Voom to calculate the differentially expressed transcripts between E3.5, E4.5, and E5.5
stage <- sampleTable$stage
design <- model.matrix(~0 + stage)
colnames(design) <- gsub("stage", "", colnames(design))

# Compare E4.5 vs E3.5 gene expression changes
contrast <- makeContrasts(
  e4.5vs3.5 = E4.5_Epi - E3.5_ICM, 
  levels = colnames(design))
contrast


# Perform Differential Expression using Limma-Voom
#filter out the lowest expressed genes using edgeR and update the normalization (both using edgeR)
keep.exprs <- filterByExpr(x, group=stage)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
x_norm <- calcNormFactors(x, method = "TMM")

#Voom converts counts to log-CPM
dge_voom <- voom(x_norm, design)

#Then fits a linear model and compares based on the given contrasts 
#In this case, E4.5_Epi vs. E3.5_ICM
dge_fit <- lmFit(dge_voom, design)
stage_fit <- contrasts.fit(dge_fit, contrasts=contrast)
efit <- eBayes(stage_fit)
dim(efit)

# Summary of differentially expressed genes
summary(decideTests(efit))

# Convert output to a data frame for further manipulation
?topTable
E4.5_vs_E3.5_voom_dataframe <- as.data.frame(topTable(efit, number=Inf))
E4.5_vs_E3.5_voom_dataframe <- setDT(E4.5_vs_E3.5_voom_dataframe, keep.rownames = TRUE)
colnames(E4.5_vs_E3.5_voom_dataframe) <- c("geneid", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")

# Add gene symbol to results data frame
E4.5_vs_E3.5_voom_dataframe <- dplyr::inner_join(E4.5_vs_E3.5_voom_dataframe, geneid_symbol, by = 'geneid')

# --------------------------------------------------------------------

# --------------------Interactive code ----------------------------------

# Exploring E4.5 vs E3.5 gene expression changes

# head command allows you to view the top several lines of a data frame
head(E4.5_vs_E3.5_voom_dataframe)

# This is the results of the RNA-seq data
# We are most interested in log2FoldChange, padj (the adjusted p-value), and the gene symbol
# This data frame contains significant and non significant results
# Lets filter for significant changes
significant_voom <- filter(E4.5_vs_E3.5_voom_dataframe, adj.P.Val < 0.05)

# nrow gives the number of rows in a data frame
age_dataframe
nrow(age_dataframe)

# 1) How many significant genes change? (hint: each row in significant dataframe is a gene)

# 2) How many significant genes are upregulated and how many are downregulated?
# (hint: filter on logFC, or go back to the summary)

# We can order data using the arrange command, which takes a data frame and what we want to order by
arrange(age_dataframe, age)
# If we want to reverse the order 
arrange(age_dataframe, -age)

# 3) What are the 10 most upregulated genes by fold change? What about by adjusted p value?

# 4) Looking at the code for the results section above, generate a data frame of
# gene expression changes for E5.5 vs E4.5

# 5) Generate a data frame of the genes of interest and plot the expression in a heatmap

# Heatmaps can be a useful way to visualize expression of genes of interest
#Convert rpkm_norm Ensembl IDs to Gene Symbols
geneid_symbol$geneid <- as.character(geneid_symbol$geneid)
rpkm_norm_symbols <- as.data.frame(cbind(rpkm_norm))
rpkm_norm_symbols$geneid <- cbind(rownames(rpkm_norm))
rpkm_norm_symbolsgeneid <- dplyr::left_join(rpkm_norm_symbols, geneid_symbol, by = "geneid", .keep = ALL)
#Choose what genes you would like to display on the heatmap, case sensitive
E3.5genesofinterest <- c("Pou5f1", "Nanog", "Sox2")
E4.5genesofinterest <- c("Klf4", "Zfp42", "Esrrb")
E5.5genesofinterest <- c("Otx2", "Fgf5", "Pou3f1")
#othergenesofinterest <- c()
genesofinterest <- as.data.frame(c(E3.5genesofinterest, E4.5genesofinterest, E5.5genesofinterest))
colnames(genesofinterest) <- c("geneid")
#Filter your RPKM data frame for select genes in genesofinterest
rpkm_norm_filtered <- dplyr::filter(rpkm_norm_symbolsgeneid, symbol %in% genesofinterest$geneid)
#convert dataframe to matrix 
rownames(rpkm_norm_filtered) <- rpkm_norm_filtered$symbol
rpkm_norm_matrix <- as.matrix(rpkm_norm_filtered[1:9])
#load package for plotting heatmaps and visualize data across all replicates
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
Heatmap(rpkm_norm_matrix)

