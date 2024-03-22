###############
#   Differential gene expression analysis with DESeq2
#   Lea Huber
#   Adapted from Rachel Stuart: Analysis of RNAseq data, 2024 Workshop on Genomics, Cesky Krumlov
#   February 2024

setwd("/faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/analysis/testsamples/DGE")
# Packages
library("DESeq2")
library("topGO")
library("tidyverse")
library("RColorBrewer")
library("ggpubr")
library("ashr")
library("org.Gg.eg.db")


#### load, convert and normalise data ---------------------------------------------------------------------------

# make sample data for small sample set, I think a data.frame?
#sampledata <- data.frame(ID = c("B2002","B2009","B2019","B2038"), condition = c("thermoneutral","heat","heat","thermoneutral"))

IDcond = read.table("SampleNames.txt", col.names = c("ID","extraction"))
IDcond$extraction = rep(c("later", "tempus"), 5)
countmatrixnames = read.table("../analysis/testsamples/mapping/matrixFileNamesCopy.txt", col.names = c("filename"))


## construct count matrix using DESeq2
# sample table for DESeqDataSet (dds) function
dds_table = data.frame(ID = IDcond$ID,counts = countmatrixnames$filename,
                       treatment = as.factor(IDcond$extraction))
## make dds
dds <- DESeqDataSetFromHTSeqCount(sampleTable = dds_table, 
                                          directory = "/faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/analysis/testsamples/mapping/",
                                          design = ~ extraction)
print(dds)

## normalisation
# Estimate the size factors
dds_norm <- estimateSizeFactors(dds)
print("Size factors:")
print(sizeFactors(dds_norm)) # Check the size factors: everything between 0.74 and 1.37 so good
print("max size factors:")
print(max(sizeFactors(dds_norm)))

#### QC and exploratory analysis --------------------------------------------

### unsupervised clustering
# transformation needed for unsupervised clustering, here variance stabilizing (vst)
vst <- vst(dds_norm, blind=TRUE) # blind means without considering design formula, so without considering condition
# set blind to False if conditions are expected to be VERY different, like almost not comparable

## PCA with 1000 most variable genes
PCAplot <- plotPCA(vst, intgroup=c("extraction"), ntop = 1000) +
  theme_pubr() + # plotPCA depends on ggplot2
  scale_color_brewer(palette = "Dark2")

## hierarchical clustering
# Create the distance matrix
vst_dists <- dist(t(assay(vst)))
vst_mat <- as.matrix(vst_dists)
# run clustering analysis
vst_hc <- hclust(vst_dists)
# visualise distance btw samples using heatmap
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) # create a colour palette 
heatMap <- heatmap(vst_mat, Rowv=as.dendrogram(vst_hc), symm=TRUE, col = rev(hmcol),
                   margin=c(13, 13), cexRow = 0.75, cexCol = 0.75)
# colour names by treatment, or have treatment names instead of sample names

#### Differential gene expression analysis ---------------------------------------------

# DESeq() does all in one including normalisation, so use raw dds object (or perhaps a filtered one)
dds <- DESeq(dds)
print("dds results:")
print(results(dds))
### hot vs. benign
## identifying DE genes using Wald test
# check what contrasts are even possible for identifying DE genes
print("possible contrasts:")
print(resultsNames(dds))
# extract result for specified contrasts
res_table <- results(dds, name = "extraction_later_vs_tempus")
# sort results based on the adjusted p-value
res_table <- res_table[order(res_table$padj),]
# You will now have the smallest padj values at the top of your table
print("head of res table")
print(head(res_table))
# save your results 
write_tsv(as.data.frame(res_table), "res_table_L_T.tsv")
