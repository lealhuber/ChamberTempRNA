### DGE models with glmmSeq
### Lea Huber

setwd("/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/data/working/RNAseq_countdata/")
# Packages
library(edgeR)
library(tidyverse)
library(ggpubr)
library(glmmSeq)
library(gplots)

tempcols = c(benign = "#6800a4", cold = "#377EB8", hot = "#E41A1C") # or should I do purple for benign?
agecols = c("1week" = "#1B9E77" , "8week" = "#D95F02")
sexcols = c(M = "#7570B3", F = "#E7298A")

#### exploration ---------------------------------------------------------------------------
# import ingredients
load("glmmIngredients_withchicks_rob.RData")

## do group or Date have effect on expression?
counts_sf <- counts*sizeFactors
colnames(counts_sf) == sampleInfo$SampleNr
sampleInfo$meanCount <- colMeans(counts_sf)
ggplot(sampleInfo, aes(x = treatment, y = meanCount))+
  geom_boxplot(outliers = FALSE)+
  geom_point()+
  theme_pubr()
# well it's a bit different but not that different visually

### PCA
pca_mat <- prcomp(t(counts), scale. = TRUE)
pc_scores <- as_tibble(pca_mat$x, rownames = "SampleNr")
pc_scores <- left_join(pc_scores, select(sampleInfo, c(treatment, group, SampleNr, age, sex)))
ggplot(pc_scores, aes(x = PC3, y = PC2, colour = treatment, shape = age)) +
  geom_point(size = 3)+
  scale_color_manual(values = tempcols)+
  theme_pubr(base_size = 22)


#### fit models -----------------------------------------------------------------------------------------------------------

# take out the chicks that were not in all treatments (just for categorical)
oncechicks <- c("7742","8235","8238")
SI_sub <- droplevels(filter(sampleInfo, !(chickno %in% oncechicks)))
# check
table(SI_sub$group, SI_sub$chickno)
table(SI_sub$group, SI_sub$treatment)
table(SI_sub$chickno, SI_sub$treatment)
table(SI_sub$group, SI_sub$sex)
table(SI_sub$treatment, SI_sub$sex)
# all good
counts_sub <- counts[,which(colnames(counts) %in% SI_sub$SampleNr)]
sizeFactors_sub <- sizeFactors[which(names(sizeFactors) %in% SI_sub$SampleNr)]


### full models
## using normalised raw counts
fitfull_gr <- glmmSeq(~ chamberTscaled + age + sex + age:chamberTscaled + (chamberTscaled | chickno) + (1|group), # takes about 8 minutes
                         countdata = counts, metadata = sampleInfo,
                         dispersion = disp, sizeFactors = sizeFactors, progress = TRUE)

fit_lin_noint <- glmmSeq(~ chamberTscaled + age + sex + (chamberTscaled | chickno) + (1|group), # takes about 9 minutes
                      countdata = counts, metadata = sampleInfo,
                      dispersion = disp, sizeFactors = sizeFactors, progress = TRUE)

# using same samples as categorical
fitfull_gr_sub <- glmmSeq(~ chamberTscaled + age + sex + age:chamberTscaled + (chamberTscaled | chickno) + (1|group), # takes about 8 minutes
                      countdata = counts_sub, metadata = SI_sub,
                      dispersion = disp, sizeFactors = sizeFactors_sub, progress = TRUE)

fit_lin_noint_sub <- glmmSeq(~ chamberTscaled + age + sex + (chamberTscaled | chickno) + (1|group), # takes about 9 minutes
                         countdata = counts_sub, metadata = SI_sub,
                         dispersion = disp, sizeFactors = sizeFactors_sub, progress = TRUE)



## categorical
fit_cat <- glmmSeq(~ treatment + age + sex + age:treatment + (treatment | chickno)+ (1|group), 
                   countdata = counts_sub, metadata = SI_sub,
                   dispersion = disp, sizeFactors = sizeFactors_sub, progress = TRUE)
fit_cat_par <- glmmSeq(~ treatment -1 + age + sex + age:treatment + (treatment | chickno)+ (1|group), 
                   countdata = counts_sub, metadata = SI_sub,
                   dispersion = disp, sizeFactors = sizeFactors_sub, progress = TRUE)

# test what happens without interaction to verify the way it works
fit_cat_noint <- glmmSeq(~ treatment + age + sex + (treatment | chickno)+ (1|group), 
                   countdata = counts_sub, metadata = SI_sub,
                   dispersion = disp, sizeFactors = sizeFactors_sub, progress = TRUE)
fit_cat_noint_par <- glmmSeq(~ treatment - 1 + age + sex + (treatment | chickno) + (1|group), 
                         countdata = counts_sub, metadata = SI_sub,
                         dispersion = disp, sizeFactors = sizeFactors_sub, progress = TRUE)
# check how glmmSeq/lmer interprets (treatment | chickno) (crossed vs nested random effects)


save(fitTchange, fitfull_gr, fitquad_gr, fit_lin_noint, fit_lin_noint_sub, fitfull_gr_sub, file = "glmmSeq_TchangeMods_andr.RData")
save(fitNormfull_cvb,fitNormfull_hvb,fit_cat, fit_cat_noint, fit_cat_noint_par, fit1w_cat, fit1w_cat_par, fit8w_cat, fit8w_cat_par,fit_cat_par, file = "glmmSeq_catmods_group.RData")
save(fit1w_gr, fit8w_gr, file = "glmmSeq_AgeSep_Oct.RData")

#### diagnostic plots ----------------------------------------------------------------------------------------------------------------

### LogFC of raw data
Sbgn <- as.character(droplevels(filter(sampleInfo, treatment == "benign")$SampleNr))
Shot <- as.character(droplevels(filter(sampleInfo, treatment == "hot")$SampleNr))
Scold <- as.character(droplevels(filter(sampleInfo, treatment == "cold")$SampleNr))
avlogexp <- rowMeans(log2(counts))
logFC_hvb <- log2(rowMeans(counts[,which(colnames(counts) %in% Shot)])) - log2(rowMeans(counts[,which(colnames(counts) %in% Sbgn)]))
logFC_cvb <- log2(rowMeans(counts[,which(colnames(counts) %in% Scold)])) - log2(rowMeans(counts[,which(colnames(counts) %in% Sbgn)]))
logFC_hvc <- log2(rowMeans(counts[,which(colnames(counts) %in% Shot)])) - log2(rowMeans(counts[,which(colnames(counts) %in% Scold)]))
plot(avlogexp, logFC_hvc)
plot(logFC_hvb,logFC_cvb)
abline(v=0)
abline(h=0)
