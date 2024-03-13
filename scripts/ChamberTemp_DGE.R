###############
#   Differential gene expression analysis with DESeq2
#   Lea Huber
#   Adapted from Rachel Stuart: Analysis of RNAseq data, 2024 Workshop on Genomics, Cesky Krumlov
#   February 2024

setwd("/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/data/working/")
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

IDcond = read.table("SampleNames.txt", col.names = c("ID","temperature","solution"))
IDcond$rep = sub("B...._(rep.)","\\1" ,IDcond$ID, perl = TRUE)
countmatrixnames = read.table("matrixFileNames.txt", col.names = c("filename","nix"))

## construct count matrix using DESeq2
# sample table for DESeqDataSet (dds) function
dds_table = data.frame(ID = IDcond$ID,counts = countmatrixnames$filename,
                       treatment = as.factor(IDcond$treatment),  age = as.factor(IDcond$age), cloacat = IDcond$cloacat,
                       weight = IDcond$weight, NeckLength = IDcond$NeckLength, growthrate = IDcond$growthrate,
                       chickno = as.factor(IDcond$chickno), sex = as.factor(IDcond$sex), group = as.factor(IDcond$group),
                       Date = IDcond$Date, time = IDcond$time, run = as.factor(IDcond$run))
## make dds
dds <- DESeqDataSetFromHTSeqCount(sampleTable = dds_table, 
                                          directory = ".",
                                          design = ~ cloacat + age + weight + NeckLength +
                                    treatment:cloacat + treatment:age +  treatment:weight + treatment:NeckLength + treatment)
dds

## optional filtering: only keep features with at least 5 reads in 3 samples
# Filtering to remove features with uniformly low read abundance is straightforward
# and has been shown to improve the detection of true differential expression (Bourgon, R., Gentleman, R. & Huber, W.
# Independent filtering increases detection power for high-throughput experiments. Proc. Natl Acad. Sci. USA 107, 9456–9551 (2010).)

# determine whether genes meet our threshold
keep <- rowSums(counts(dds) >= 5 ) >= 3 # with broiler data that's everything
# filter the dds object to keep genes that meet the threshold
dds_filt <- dds[keep,]
# if you plan to use the filtered data set, make sure you use the filtered dds object in the following analysis

## normalisation
# Estimate the size factors
dds_norm <- estimateSizeFactors(dds) 
sizeFactors(dds_norm) # Check the size factors: everything between 0.74 and 1.37 so good
max(sizeFactors(dds_norm))

#### QC and exploratory analysis --------------------------------------------

### unsupervised clustering
# transformation needed for unsupervised clustering, here variance stabilizing (vst)
vst <- vst(dds_norm, blind=TRUE) # blind means without considering design formula, so without considering condition
# set blind to False if conditions are expected to be VERY different, like almost not comparable

## PCA with 1000 most variable genes
plotPCA(vst, intgroup=c("treatment"), ntop = 1000) +
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
heatmap(vst_mat, Rowv=as.dendrogram(vst_hc), symm=TRUE, col = rev(hmcol),
                   margin=c(13, 13), cexRow = 0.75, cexCol = 0.75)
# colour names by treatment, or have treatment names instead of sample names

#### Differential gene expression analysis ---------------------------------------------

# DESeq() does all in one including normalisation, so use raw dds object (or perhaps a filtered one)
dds <- DESeq(dds)
results(dds)
### hot vs. benign
## identifying DE genes using Wald test
# check what contrasts are even possible for identifying DE genes
resultsNames(dds)
# extract result for specified contrasts
res_table_H_B <- results(dds, name = "treatment_hot_vs_benign")
# sort results based on the adjusted p-value
res_table_H_B <- res_table_H_B[order(res_table_H_B$padj),]
# You will now have the smallest padj values at the top of your table
head(res_table_H_B)
# save your results 
write_tsv(as.data.frame(res_table_H_B), "DESeq_output/res_table_H_B.tsv")

#### visualising DE genes --------------------------------------------------------------------

# visualise difference in top gene
plotCounts(dds, gene="ENSGALG00000000003", intgroup="temperature") # see that genes with a negative LFC are upregulated in heat
# MA plot of all genes
plotMA(object = dds, 
       ylim=c(-2,2), 
       alpha = 0.05)
# weird patterns at genes with low counts
# try Log fold change shrinkage for visualization and ranking to resolve it
resLFCashr <- lfcShrink(dds, coef="temperature_thermoneutral_vs_heat", type="ashr")
resLFCashr
plotMA(object = resLFCashr, 
       ylim=c(-0.3,0.3), 
       alpha = 0.05) # pattern is gone and lfc smaller

# volcano plot of all genes
ggplot() +
  geom_point(
    data = as.data.frame(res_table_H_B) %>% 
      filter(padj > 0.01 | abs(log2FoldChange) < 1), 
    aes(x = log2FoldChange, y = -log10(padj)),
    color = "grey") +
  geom_point(
    data = as.data.frame(res_table_H_B) %>%
      filter(padj < 0.01 & abs(log2FoldChange) > 1), 
    aes(x = log2FoldChange, y = -log10(padj)), 
    color = "blue") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_bw()
# save plot

# What are the names of those few blue genes?
bluegenes = rownames(filter(as.data.frame(res_table_H_B), padj < 0.01 & abs(log2FoldChange) > 1))

### cold vs. benign
dds$treatment = relevel(dds$treatment, "cold")
resultsNames(dds)
# extract result for specified contrasts
res_table_C_B <- results(dds, contrast = "treatment_cold_vs_benign") # or so
# sort results based on the adjusted p-value
res_table_C_B <- res_table_C_B[order(res_table_C_B$padj),]
# You will now have the smallest padj values at the top of your table
head(res_table_C_B)
# save your results 
write_tsv(as.data.frame(res_table_C_B), "DESeq_output/res_table_C_B.tsv")

# MA plot
pdf(file = "output/MAplot_C_B.pdf", height = 6, width = 6)
plotMA(object = dds, alpha = 0.01)
dev.off()

# volcano plot of all genes
ggplot() +
  geom_point(
    data = as.data.frame(res_table_C_B) %>% 
      filter(padj > 0.01 | abs(log2FoldChange) < 1), 
    aes(x = log2FoldChange, y = -log10(padj)),
    color = "grey") +
  geom_point(
    data = as.data.frame(res_table_C_B) %>%
      filter(padj < 0.01 & abs(log2FoldChange) > 1), 
    aes(x = log2FoldChange, y = -log10(padj)), 
    color = "blue") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_bw()
# save plot

# Identify the top 20 differentially expressed genes from the sorted results table. 
DE_genes <- rownames(head(res_table_C_B,50))
# Create a heatmap for the most expressed genes
dds_DE_B_C <- assay(vst)[rownames(Pcal_vst) %in% DE_genes,]
heatmap(dds_DE_B_C, margin = c(10,6))
# save plot
pdf(file = "output/Pcal_geneHeatmap_U_S.pdf", height = 6, width = 6)
heatmap(dds_DE_B_C,margin = c(10,6))
dev.off()

#### Compare DE genes between hot and cold
# visually: venn diagrams (e.g., gplots), euler diagrams (eulerr), upsetR plots (e.g., upsetR, ComplexHeatmap),
# alluvial plots (ggalluvial), scatterplots (e.g., ggplot2).
# same with GO terms (?)

#### Gene set enrichment analysis ------------------------------
# go annotation path (easier to use in topGO aparently)
GO_path <- "goa_chicken.gaf"
GO_annotation <- read_tsv(GO_path, comment = "!", col_names = c("DB","DB_Object_ID","DB_Object_Symbol","Qualifier","GO_ID",
                                                                "DB:Reference","Evidence Code","With (or) From","Aspect",
                                                                "DB_Object_Name","DB_Object_Synonym","DB_Object_Type",
                                                                "Taxon and Interacting taxon","Date","Assigned_By",
                                                                "Annotation_Extension","Gene_Product_Form_ID"))
# does not contain feature names as they are in dds
# how many GO terms are there
length(unique(GO_annotation$DB_Object_Name))

# here I will be only analysing GO terms with at least 5 members, as this yield more stable results.
node_size= 5

# We will focus on Biological process terms (P or BP)
GO_category="P"

# use topGO to read the functional annotation
GO2geneID <- annFUN.org("BP", mapping = "org.Gg.eg.db", ID = "genbank")

names(GO2geneID)%in%GO_annotation$GO_ID # a lot are TRUE! :D but now only get the ones who are for further analysis
GO2geneID = GO2geneID[which(names(GO2geneID)%in%GO_annotation$GO_ID)] # from 3729 down to 3670

# now somehow get geneID2GO
geneID2GO = inverseList(GO2geneID)

# define the gene universe
geneUniverse <- names(geneID2GO)
bluegenes%in%geneUniverse # shit
# the gene IDs from annFUN are not the same as from the reference genome I used although both are ensembl! :(
# so I can't connect it to the GO_annotation
table(rownames(res_table_T_H)%in%names(geneID2GO)) # see

