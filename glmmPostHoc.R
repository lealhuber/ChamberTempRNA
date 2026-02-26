### post hoc analysis of gene expression
## Oct 2025
## Lea Huber

library(tidyverse)
library(glmmSeq)
library(emmeans)
library(ggpubr)

tempcols = c(benign = "#6800a4", cold = "#377EB8", hot = "#E41A1C")
agecols = c("week1" = "#1B9E77" , "week8" = "#D95F02")


setwd("/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/data/working/RNAseq_countdata/")


#### gather ingredients -----------------------------------------------------------------------

## load models
load("glmmSeq_TchangeMods_andr.RData")
load("glmmSeq_catmods_group.RData")

### get genes that did not converge in all models
errgs_rgr <- names(fitfull_gr@errors)
errgs_lin_noint <- names(fit_lin_noint@errors)
errgs_cat <- names(fit_cat@errors)
errgs_cat_noint <- names(fit_cat_noint_par@errors)
errgs_2cat2lin <- unique(c(errgs_rgr, errgs_cat, errgs_cat_noint,errgs_lin_noint))

### for adding gene information
GO_path <- "../Struthio_camelus_HiC.emapper.annotations" # go annotation path (easier to use in topGO apparently)
GO_annotation <- read_tsv(GO_path, comment = "##", na = "-")
GO_annotation$COG_category <- as.factor(GO_annotation$COG_category) # oh there are 141
# remove transript information
GO_annotation$query <- str_remove(GO_annotation$query, ".t\\d" )
GO_annotation <- GO_annotation[!duplicated(GO_annotation$query), ] # each query only once

### get significant genes
# linear
fitfull_gr <- glmmQvals(fitfull_gr, cutoff = 0.01)
statsctgr <- summary(fitfull_gr)
statsct_agesig <- statsctgr[order(statsctgr[, 'P_chamberTscaled:age']), ]
statsct_agesig <- as.data.frame(statsct_agesig[1:229,])
colnames(statsct_agesig) <- c(colnames(statsct_agesig)[1:26], "q_chamberTscaled", "q_age", "q_sex", "q_chamberTscaled:age")
statsct_agesig$gene <- rownames(statsct_agesig)
statsct_agesig <- filter(statsct_agesig, !gene%in%errgs_2cat2lin2lin)
statsct_agesig <- left_join(statsct_agesig, subset(GO_annotation, select = c("query","Preferred_name","Description")), by = join_by(gene == query))
rownames(statsct_agesig) <- statsct_agesig$gene
# linear without interaction term
fit_lin_noint <- glmmQvals(fit_lin_noint, cutoff = 0.01)
statslin_noint <- summary(fit_lin_noint)
statslin_noint_sig <- statslin_noint[order(statslin_noint[, 'P_chamberTscaled']), ]
statslin_noint_sig <- as.data.frame(statslin_noint_sig[1:324,])
colnames(statslin_noint_sig) <- c(colnames(statslin_noint_sig)[1:21], "q_chamberTscaled", "q_age", "q_sex")
statslin_noint_sig$gene <- rownames(statslin_noint_sig)
statslin_noint_sig <- filter(statslin_noint_sig, !gene%in%errgs_2cat2lin2lin)
statslin_noint_sig <- left_join(statslin_noint_sig, subset(GO_annotation, select = c("query","Preferred_name","Description")), by = join_by(gene == query))
rownames(statslin_noint_sig) <- statslin_noint_sig$gene
# categorical in one model
fit_cat <- glmmQvals(fit_cat, cutoff = 0.01)
fit_cat_par <- glmmQvals(fit_cat_par)
fit_cat_noint <- glmmQvals(fit_cat_noint, cutoff = 0.01) # is it the same? No, but similar, and better according to LRT
fit_cat_noint_par <- glmmQvals(fit_cat_noint_par)
statscat_withint <- summary(fit_cat)
statscat_withint_sig <- statscat_withint[order(statscat_withint[, 'P_treatment:age']), ] # take interaction term!!
statscat_withint_sig <- as.data.frame(statscat_withint_sig[1:395,])
colnames(statscat_withint_sig) <- c(colnames(statscat_withint_sig)[1:30], "q_treatment", "q_age", "q_sex", "q_treatment:age")
statscat_withint_sig$gene <- rownames(statscat_withint_sig)
statscat_withint_sig <- filter(statscat_withint_sig, !(gene %in% errgs_2cat2lin)) # remove the ones that don't converge for other models I will use
statscat_withint_sig <- left_join(statscat_withint_sig, subset(GO_annotation, select = c("query","Preferred_name","Description")), by = join_by(gene == query))
rownames(statscat_withint_sig) <- statscat_withint_sig$gene
# categorical with no interaction
statsFull_cat <- summary(fit_cat_noint)
statsFull_cat_sig <- statsFull_cat[order(statsFull_cat[, 'P_treatment']), ]
statsFull_cat_sig <- as.data.frame(statsFull_cat_sig[1:472,])
colnames(statsFull_cat_sig) <- c(colnames(statsFull_cat_sig)[1:23], "q_treatment", "q_age", "q_sex")
statsFull_cat_sig$gene <- rownames(statsFull_cat_sig)
statsFull_cat_sig <- filter(statsFull_cat_sig, !(gene %in% errgs_2cat2lin)) # remove the ones that don't converge for other models I will use
statsFull_cat_sig <- left_join(statsFull_cat_sig, subset(GO_annotation, select = c("query","Preferred_name","Description")), by = join_by(gene == query))
rownames(statsFull_cat_sig) <- statsFull_cat_sig$gene

#### refit genes and do post hoc analyses -------------------------------------------------------------------------

### models with interaction
## linear
head(statsct_agesig)
statsct_agesig <- filter(statsct_agesig, !(gene %in% errgs_2cat2lin2lin)) # 224
# for each gene table with effect size in 1 week and 8 week, or non sig if not significant
# have to do post hoc to find out if both are significant or just one!
# but first try to understand what the coefficients mean
head(fitfull_gr@predict)
fitfull_gr@predict[1,1:36]
# I am pretty sure that the sign in chamberTscaled says whether it's overall going up or down (or just 1w)
# the sign in age8week shows whether the overall expression in 8w is higher or lower than 1w
# and adding chamberTscaled and chamberTscaled:age8week gives the sign of the slope of 8w
# want to make a table like this:
## gene # age # slope #
## g1   # 1w  #  0.2  #
## g1   # 8w  #  0.3  #
## g2   # 1w  # -0.1  #
## g2   # 8w  #  ns   #
# and I've almost got that, for slope 8w just need to add interaction term and the other is cTs!
# and I only took the ones that are sig for interaction term so that means that ages are sig different!
# but do I have to do postHoc to find out if cTs (i.e. 1w) is significant?
# this gives me the slopes and if ages are sigificantly different
# but not if ages are signficant by themselves, and if overall trend is significant
# but that I get in the summary actually! so I only need the gene fits here!
lin_int_genes <- as.list(rownames(statsct_agesig))
names(lin_int_genes) <- rownames(statsct_agesig)
GeneFitsLinInt <- lapply(lin_int_genes, function(x) suppressMessages(glmmRefit(fitfull_gr, gene = x))) # takes <1 minutes
summary(GeneFitsLinInt$g10974)
# gather ingredients for data frame
gene <- names(GeneFitsLinInt)
intercept1w <- sapply(GeneFitsLinInt, function(M) summary(M)$coefficients[1,1])
intchange8w <- sapply(GeneFitsLinInt, function(M) summary(M)$coefficients[3,1])
slope1w <- sapply(GeneFitsLinInt, function(M) summary(M)$coefficients[2,1])
slope8w <- sapply(GeneFitsLinInt, function(M) summary(M)$coefficients[5,1])
p1w <- sapply(GeneFitsLinInt, function(M) summary(M)$coefficients[2,4])
p_agediff <- sapply(GeneFitsLinInt, function(M) summary(M)$coefficients[3,4]) # also interesting if 8w is overall sig different from 1w (intercept not slope)
# but I also want to know HO the slope in 8w is 0 (which is what the chamberTscaled pval is for 1w so I got that)
# how do I get that?
library(car)
linhop <- linearHypothesis(fit_g10974, "chamberTscaled + chamberTscaled:age8week = 0")
linhop$`Pr(>Chisq)`[2] # voilà!
# for all
p8w <- sapply(GeneFitsLinInt, function(M) linearHypothesis(M, "chamberTscaled + chamberTscaled:age8week = 0")$`Pr(>Chisq)`[2]) 
linear_effects <- data.frame(gene=gene,intercept1w=intercept1w,intercept8w = intercept1w+intchange8w, slope1w=slope1w,slope8w=slope1w+slope8w,p1w=p1w,p8w=p8w,p_agediff = p_agediff)
summary(linear_effects) # looks promising
# now it's just a matter of reshaping
linear_effects_ln <- pivot_longer(linear_effects, cols = intercept1w:p8w, names_to = c("coef","age"),
                                  names_sep = -2, values_to = "value")
linear_effects <- pivot_wider(linear_effects_ln, names_from = coef, values_from = value)
linear_effects <- rename(linear_effects, pval = p)
linear_effects[which(linear_effects$pval > 0.01), "slope"] <- NA # where non sig change slope to NA
linear_effects <- filter(linear_effects, !is.na(slope)) # can also just remove those for less confusion
summary(linear_effects)
table(linear_effects$age, is.na(linear_effects$slope)) # 196 sig genes for 1w and 137 for 8w
table(linear_effects$p_agediff<0.01) # for 388/2 the intercept is sig different between ages, so most
# how many genes have significant slope in both ages?
table(duplicated(filter(linear_effects, !is.na(slope))$gene)) # 114 are in both
# still interesting to know if the slopes are sig different btw ages (should be because sig for interaction all of them)
EMres_lin <- lapply(GeneFitsLinInt, function(fit) emtrends(fit, specs = "age", var = "chamberTscaled"))
p_int_lin <- sapply(EMres_lin, function(emm) as.data.frame(pairs(emm))$p.value)
summary(p_int_lin) # yes they are indeed all very significant

save(GeneFitsLinInt, linear_effects, file = "PostHoc_linInt.RData")



## categorical
head(statscat_withint_sig)
statscat_withint_sig <- filter(statscat_withint_sig, !(gene %in% errgs_2cat2lin))
# for each gene table with effect size in 1 week and 8 week X hvc and hvb
## gene # age # hvb  # cvb #
## g1   # 1w  #  0.2 # 0.3 #
## g1   # 8w  #  0.3 # ns  #
## g2   # 1w  # -0.1 # 0.1 #
## g2   # 8w  #  ns  # ns  #

# here I need post-hoc
genelist <- as.list(rownames(statscat_withint_sig))
names(genelist) <- rownames(statscat_withint_sig)
GeneFits_catint <- lapply(genelist, function(x) suppressMessages(glmmRefit(fit_cat_par, gene = x)))  
genefit <- GeneFits_catint[[3]] # inspect
summary(genefit)
# here I want the estimates COMPARED TO BENIGN for hot in age 1, hot in age 8 etc
EMM <- emmeans(genefit, ~ treatment * age)
summary(EMM)
empairs <- pairs(EMM, simple = "treatment")    # compare treats for each age, this is exactly what I need
summary(empairs)
pairs(EMM, simple = "age")     # compare ages for each treat
pairs(EMM)                   # compares all combinations (also ones I am not interested in)
# do emmeans for all genes
EMres_cat <- lapply(GeneFits_catint, function(fit) emmeans(fit, ~ treatment * age))
EMpairs_cat <- lapply(EMres_cat, function(emm) pairs(emm, simple = "treatment"))
w1_est <- sapply(EMpairs_cat, function(emm) summary(emm)$estimate[c(1,2)])
w8_est <- sapply(EMpairs_cat, function(emm) summary(emm)$estimate[c(4,5)])
w1_est <- t(w1_est)
w8_est <- t(w8_est)
w1_pval <- sapply(EMpairs_cat, function(emm) summary(emm)$p.value[c(1,2)])
w8_pval <- sapply(EMpairs_cat, function(emm) summary(emm)$p.value[c(4,5)])
w1_pval <- t(w1_pval)
w8_pval <- t(w8_pval)
cat_effects <- data.frame(gene = rownames(w1_est), age = rep("week1",nrow(w1_est)),cvb = w1_est[,1], hvb = w1_est[,2],
                          cvb_p = w1_pval[,1],hvb_p = w1_pval[,2], row.names = NULL)
cat_effects8 <- data.frame(gene = rownames(w8_est),age = rep("week8",nrow(w8_est)), cvb = w8_est[,1], hvb = w8_est[,2],
                           cvb_p = w8_pval[,1],hvb_p = w8_pval[,2], row.names = NULL)
cat_effects <- rbind(cat_effects,cat_effects8)
summary(cat_effects)
# make one with all values wide
cat_eff_wide_all <- pivot_wider(select(cat_effects, gene:hvb), names_from = age, values_from = c(cvb,hvb))
# make non significant effects NA
cat_effects[which(cat_effects$cvb_p > 0.01),"cvb"] <- NA
cat_effects[which(cat_effects$hvb_p > 0.01),"hvb"] <- NA
nocateff <-  filter(cat_effects, is.na(cvb) & is.na(hvb))$gene
cat_effects <- filter(cat_effects, !(is.na(cvb) & is.na(hvb))) # remove rows with NA in both effects, removes a lot actually I guess thats hvc
cat_eff_wide_all <- filter(cat_eff_wide_all, !(gene %in% nocateff))
summary(cat_effects)
table(cat_effects$age, is.na(cat_effects$cvb)) # 136 for 1w and 100 for 8w
table(cat_effects$age, is.na(cat_effects$hvb)) # 204 for 1w and 128 for 8w
# how many hvb are significant in both ages?
table(duplicated(filter(cat_effects, !is.na(hvb))$gene)) # 110
# how many cvb are significant in both ages?
table(duplicated(filter(cat_effects, !is.na(cvb))$gene)) # 74

# check if ages are sig different in both heat and cold
test <- pairs(EMM, simple = "age") 
summary(test)[3,6]
EMpairs_cat_age <- lapply(EMres_cat, function(emm) pairs(emm, simple = "age"))
int_pval_cold <- sapply(EMpairs_cat_age, function(emm) summary(emm)[2,7])
int_pval_hot <- sapply(EMpairs_cat_age, function(emm) summary(emm)[3,7])
summary(int_pval_cold) # not all sig
summary(int_pval_hot) # not all sig
int_pval_cold <- int_pval_cold[which(int_pval_cold <= 0.01)]
int_pval_hot <- int_pval_hot[which(int_pval_hot <= 0.01)]

save(GeneFits_catint, cat_effects, file = "PostHoc_catInt.RData")


## how many are there?
length(unique(c(rownames(statsct_agesig), rownames(statscat_withint_sig)))) # 557



### models without interaction:
## linear
head(statslin_noint_sig)
statslin_noint_sig <- filter(statslin_noint_sig, !(gene %in% errgs_2cat2lin2lin))
lin_noint_genes <- as.list(rownames(statslin_noint_sig))
names(lin_noint_genes) <- rownames(statslin_noint_sig)
GeneFitsLinNoint <- lapply(lin_noint_genes, function(x) suppressMessages(glmmRefit(fit_lin_noint, gene = x))) # takes <1 minutes
summary(GeneFitsLinNoint[[1]])
# gather ingredients for data frame
gene <- names(GeneFitsLinNoint)
slope <- sapply(GeneFitsLinNoint, function(M) summary(M)$coefficients[2,1])
diff8w <- sapply(GeneFitsLinNoint, function(M) summary(M)$coefficients[3,1])
p_agediff <- sapply(GeneFitsLinNoint, function(M) summary(M)$coefficients[3,4]) # are 1w and 8w significantly different? That's the age p value
linear_effects_noint <- data.frame(gene = gene, slope = slope, diff8w = diff8w, p_agediff = p_agediff)
summary(linear_effects_noint)
table(linear_effects_noint$p_agediff<0.01) # 194 have significant age difference in intercept (same slope in this model)


## categorical
head(statsFull_cat_sig)
statsFull_cat_sig <- filter(statsFull_cat_sig, !(gene %in% errgs_2cat2lin))
genelist <- as.list(rownames(statsFull_cat_sig))
names(genelist) <- rownames(statsFull_cat_sig)
GeneFits_catNoint <- lapply(genelist, function(x) suppressMessages(glmmRefit(fit_cat_noint_par, gene = x))) # takes 
genefit <- GeneFits_catNoint[[1]] # inspect
summary(genefit)
# here I want the estimates COMPARED TO BENIGN for hot and cold (no matter age)
EMres_trt <- lapply(GeneFits_catNoint, function(fit) emmeans(fit, ~ treatment))
EMpairs_trt <- lapply(EMres_trt, function(emm) pairs(emm))
summary(EMpairs_trt[[1]])
cvb_est_noint <- sapply(EMpairs_trt, function(emm) summary(emm)$estimate[1])
hvb_est_noint <- sapply(EMpairs_trt, function(emm) summary(emm)$estimate[2])
cvb_p_noint <- sapply(EMpairs_trt, function(emm) summary(emm)$p.value[1])
hvb_p_noint <- sapply(EMpairs_trt, function(emm) summary(emm)$p.value[2])
gene <- names(EMpairs_trt)
cat_effects_noint <- data.frame(gene=gene,cvb=cvb_est_noint,hvb=hvb_est_noint, p_cvb=cvb_p_noint,p_hvb=hvb_p_noint)
# make df with all values
cat_effects_both_all <- full_join(cat_eff_wide_all,select(cat_effects_noint, gene:hvb))
# replace non-sig with NA
cat_effects_noint[which(cat_effects_noint$p_cvb>0.01),"cvb"] <- NA
cat_effects_noint[which(cat_effects_noint$p_hvb>0.01),"hvb"] <- NA
nocateff <- filter(cat_effects_noint, is.na(cvb) & is.na(hvb))$gene # genes that have both NA so no effect after all
cat_effects_noint <- filter(cat_effects_noint, !(is.na(cvb) & is.na(hvb))) # remove rows where both NA
cat_effects_both_all <- filter(cat_effects_both_all, !(gene %in% nocateff))
# hm that removed quite a few like half why were they even in there? prob hvc
summary(cat_effects_noint) # 434 sig for hvb and 109 sig for cvb


### now also make one for all genes in the analysis
genelist <- as.list(unique(c(cat_effects$gene, cat_effects_noint$gene,linear_effects$gene,linear_effects_noint$gene)))
names(genelist) <- unique(c(cat_effects$gene, cat_effects_noint$gene,linear_effects$gene,linear_effects_noint$gene))
GeneFits_all <- lapply(genelist, function(x) suppressMessages(glmmRefit(fit_cat_noint_par, gene = x))) # takes a whiile
lingene <- GeneFits_all$g10974
emlingene <- emmeans(lingene, ~treatment)
emplingene <- pairs(emlingene)
summary(emplingene)
EMres_all <- lapply(GeneFits_all, function(fit) emmeans(fit, ~ treatment))
EMpairs_all <- lapply(EMres_all, function(emm) pairs(emm))
summary(EMpairs_all[[1]])
cvb_est_noint <- sapply(EMpairs_all, function(emm) summary(emm)$estimate[1])
hvb_est_noint <- sapply(EMpairs_all, function(emm) summary(emm)$estimate[2])
cvb_p_noint <- sapply(EMpairs_all, function(emm) summary(emm)$p.value[1])
hvb_p_noint <- sapply(EMpairs_all, function(emm) summary(emm)$p.value[2])
gene <- names(EMpairs_all)
cat_effects_all <- data.frame(gene=gene,cvb=cvb_est_noint,hvb=hvb_est_noint, p_cvb=cvb_p_noint,p_hvb=hvb_p_noint)


save(cat_effects_noint, linear_effects_noint, GeneFits_catNoint, GeneFitsLinNoint,GeneFits_all, cat_effects_all, file = "PostHoc_NoInt.RData")



### sanity check: are fits for interaction and non-interaction models similar?
table(cat_effects$gene %in% cat_effects_noint$gene)
table(cat_effects_noint$gene %in% cat_effects$gene) # do not overlap a lot

plot(rowMeans(cbind(cat_effects_both_all$cvb_week8,cat_effects_both_all$cvb_week1)),cat_effects_both_all$cvb)
abline(0,1)
plot(rowMeans(cbind(cat_effects_both_all$hvb_week8,cat_effects_both_all$hvb_week1)),cat_effects_both_all$hvb)
abline(0,1)
# yes they are similar


#### get numbers and overlaps -----------------------------------------------------------------

# load results
load("PostHoc_catInt.RData")
load("PostHoc_linInt.RData")
load("PostHoc_NoInt.RData")

# get genes
genes_lin1w <- filter(linear_effects, age == "1w" & !is.na(slope))$gene
genes_lin8w <- filter(linear_effects, age == "8w" & !is.na(slope))$gene
genes_linNoInt <- linear_effects_noint$gene
genes_hvb1w <- filter(cat_effects, age == "week1" & !is.na(hvb))$gene
genes_hvb8w <- filter(cat_effects, age == "week8" & !is.na(hvb))$gene
genes_hvbNoInt <- filter(cat_effects_noint, !is.na(hvb))$gene
genes_cvb1w <- filter(cat_effects, age == "week1" & !is.na(cvb))$gene
genes_cvb8w <- filter(cat_effects, age == "week8" & !is.na(cvb))$gene
genes_cvbNoInt <- filter(cat_effects_noint, !is.na(cvb))$gene

tempgenes3 <- unique(c(cat_effects$gene, cat_effects_noint$gene))


# categorise hot, cold and shared
make_colvec <- function(x){
  names(x) <- x
  if (x %in% genes_hvbNoInt & x %in% genes_cvbNoInt){
    return("hotandcold")
  }else if (x %in% genes_hvbNoInt){
    return("hot")
  }else if (x %in% genes_cvbNoInt){
    return("cold")
    #  }else if(x %in% lin_genes){
    #    return("linear")
  }else{
    return("out")
  }
}
colvec <- sapply(cat_effects_noint$gene, function(x) make_colvec(x)) # changed here from cat_effects_all to cat_effects_noint because we're not using linear anymore
table(colvec)
cat_effects_noint$colour <- colvec
# write.csv(cat_effects_noint, file = "cat_effects_noint.csv") # save for Mads


## up or downregulated
hotdown <- filter(cat_effects_noint, hvb > 0 & colour == "hot") 
hotup <- filter(cat_effects_noint, hvb < 0 & colour == "hot") 
colddown <- filter(cat_effects_noint, cvb > 0 & colour == "cold") 
coldup <- filter(cat_effects_noint, cvb < 0 & colour == "cold")

### Are the ones shared between hot and cold, same or opposite direction? (using also non-significant opposite genes)
hotdowncoldup <- filter(cat_effects_noint, hvb > 0 & cvb < 0 & colour == "hotandcold") 
hotupcolddown <- filter(cat_effects_noint, hvb < 0 & cvb > 0 & colour == "hotandcold") 
bothdown <- filter(cat_effects_noint, hvb > 0 & cvb > 0 & colour == "hotandcold") 
bothup <- filter(cat_effects_noint, hvb < 0 & cvb < 0 & colour == "hotandcold")


# and are same direction ones the ones that are also in linear?
table(c(bothdown$gene,bothup$gene) %in% unique(c(linear_effects$gene, linear_effects_noint$gene)))
table(c(hotupcolddown$gene,hotdowncoldup$gene) %in% unique(c(linear_effects$gene, linear_effects_noint$gene)))
# more than the smiley ones but still less than half...
# to count as shared they still have to be significant in both heat and cold
gplots::venn(list("cat opposite" = c(hotupcolddown$gene,hotdowncoldup$gene),
                  "hvb" = filter(cat_effects, !is.na(hvb))$gene,
                  "cvb" = filter(cat_effects, !is.na(cvb))$gene))


## how many are different in ages over all?
length(unique(c(cat_effects$gene,linear_effects$gene)))
nrow(linear_effects)/length(unique(c(linear_effects_noint$gene,linear_effects$gene)))*100

## compare intensity of age responses
# how many each age?
length(unique(c(genes_lin1w,genes_hvb1w,genes_cvb1w))) # 368
length(unique(c(genes_lin8w,genes_hvb8w,genes_cvb8w))) # 276
# effect intensity?
hist(filter(cat_effects, age == "week1")$hvb)
hist(filter(cat_effects, age == "week8")$hvb)
hist(filter(cat_effects, age == "week1")$cvb)
hist(filter(cat_effects, age == "week8")$cvb)
hist(filter(linear_effects, age == "1w")$slope)
hist(filter(linear_effects, age == "8w")$slope)
# higher in 1w a bit



### save gene lists for functional analyses
# ages: genes that are differently or same expressed between ages
ages_different <- unique(c(filter(cat_eff_wide, (hvb_week1 < 0 & hvb_week8 > 0) | (hvb_week1 > 0 & hvb_week8 < 0))$gene,
                           filter(cat_eff_wide, (cvb_week1 < 0 & cvb_week8 > 0) | (cvb_week1 > 0 & cvb_week8 < 0))$gene,
                           filter(lin_eff_wide, (week1 < 0 & week8 > 0) | (week1 > 0 & week8 < 0))$gene))
ages_different_cat <- unique(c(filter(cat_eff_wide, (hvb_week1 < 0 & hvb_week8 > 0) | (hvb_week1 > 0 & hvb_week8 < 0))$gene,
                           filter(cat_eff_wide, (cvb_week1 < 0 & cvb_week8 > 0) | (cvb_week1 > 0 & cvb_week8 < 0))$gene))
ages_same <- unique(c(filter(cat_eff_wide, (hvb_week1 < 0 & hvb_week8 < 0) | (hvb_week1 > 0 & hvb_week8 > 0))$gene,
                      filter(cat_eff_wide, (cvb_week1 < 0 & cvb_week8 < 0) | (cvb_week1 > 0 & cvb_week8 > 0))$gene,
                      filter(lin_eff_wide, (week1 < 0 & week8 < 0) | (week1 > 0 & week8 > 0))$gene))
ages_same_cat <- unique(c(filter(cat_eff_wide, (hvb_week1 < 0 & hvb_week8 < 0) | (hvb_week1 > 0 & hvb_week8 > 0))$gene,
                      filter(cat_eff_wide, (cvb_week1 < 0 & cvb_week8 < 0) | (cvb_week1 > 0 & cvb_week8 > 0))$gene))

save(hotup,hotdown,coldup,colddown,hotdowncoldup,hotupcolddown,bothup,bothdown,errgs_2cat2lin, file = "NewPHgenesForGO.RData") # should only be no interaction genes!

# also make list of genes for Fst analysis
genes_of_interest_noint <- cat_effects_noint$gene
table(duplicated(genes_of_interest_noint)) # good no overlaps
writeLines(genes_of_interest_noint, "../genes_of_interest_noint.txt")

#### venn plots -------------------------

gplots::venn(list("hvb" = filter(cat_effects, !is.na(hvb))$gene,
                  "cvb" = filter(cat_effects, !is.na(cvb))$gene))
gplots::venn(list("linear" = linear_effects$gene,
                  "hvb" = filter(cat_effects, !is.na(hvb))$gene,
                  "cvb" = filter(cat_effects, !is.na(cvb))$gene))
# gplots::venn(list("lin1w" = filter(linear_effects, age == "1w" & !is.na(slope))$gene,
#                   "lin8w" = filter(linear_effects, age == "8w" & !is.na(slope))$gene,
#                   "hvb1w" = filter(cat_effects, age == "week1" & !is.na(hvb))$gene,
#                   "hvb8w" = filter(cat_effects, age == "week8" & !is.na(hvb))$gene))
gplots::venn(list("cvb1w" = filter(cat_effects, age == "week1" & !is.na(cvb))$gene,
                  "cvb8w" = filter(cat_effects, age == "week8" & !is.na(cvb))$gene,
                  "hvb1w" = filter(cat_effects, age == "week1" & !is.na(hvb))$gene,
                  "hvb8w" = filter(cat_effects, age == "week8" & !is.na(hvb))$gene))
# from models without interaction: shared and distinct
# gplots::venn(list("linear" = linear_effects_noint$gene,
#                   "hvb" = filter(cat_effects_noint, !is.na(hvb))$gene,
#                   "cvb" = filter(cat_effects_noint, !is.na(cvb))$gene))
# compare mods with and without interaction
gplots::venn(list("lin_int" = linear_effects$gene,
                  "lin_noint" = linear_effects_noint$gene,
                  "cat_int" = cat_effects$gene,
                  "cat_noint" = cat_effects_noint$gene))

gplots::venn(list("cat_int" = cat_effects$gene,
                  "cat_noint" = cat_effects_noint$gene))

# quite different
# shared and distinct 1w
gplots::venn(list("lin1w" = filter(linear_effects, age == "1w" & !is.na(slope))$gene,
                  "hvb1w" = filter(cat_effects, age == "week1" & !is.na(hvb))$gene,
                  "cvb1w" = filter(cat_effects, age == "week1" & !is.na(cvb))$gene))
# shared and distinct 8
gplots::venn(list("lin8w" = filter(linear_effects, age == "8w" & !is.na(slope))$gene,
                  "hvb8w" = filter(cat_effects, age == "week8" & !is.na(hvb))$gene,
                  "cvb8w" = filter(cat_effects, age == "week8" & !is.na(cvb))$gene))
#### plot results --------------------------

### plot showing amount of shared and distinct response

cat_effects_noint$resptype <- ifelse(cat_effects_noint$colour %in% c("cold","hot"), "distinct","shared")
distshared_freq <- as.data.frame(table(cat_effects_noint$colour,cat_effects_noint$resptype))
distshared_df <- data.frame(cond = c(rep("cold",distshared_freq$Freq[1]+distshared_freq$Freq[6]), rep("hot",distshared_freq$Freq[2]+distshared_freq$Freq[6])),
                            resptype = c(rep("distinct",distshared_freq$Freq[1]),rep("shared",distshared_freq$Freq[6]),rep("distinct",distshared_freq$Freq[2]),rep("shared",distshared_freq$Freq[6])),
                            label = c(rep(paste0("n = ",distshared_freq$Freq[1]+distshared_freq$Freq[6]),distshared_freq$Freq[1]+distshared_freq$Freq[6]),
                                      rep(paste0("n = ",distshared_freq$Freq[2]+distshared_freq$Freq[6]),distshared_freq$Freq[2]+distshared_freq$Freq[6])))
ggplot(distshared_df, aes(x = cond, fill = resptype))+
  geom_bar(position = "fill")+
  scale_fill_manual(values = c(gray(0.4),gray(0.8)))+
  geom_text(aes(label = label, y = 1.05))+
  scale_fill_manual(values = c("#f98400", "#00a08a"))+
  xlab("")+ylab("proportion of genes")+
  theme_pubr(base_size = 20)+
  theme(legend.title = element_blank())
# ggsave("Bar_shareddist.pdf",path = "../../../figures/glmmSeq_PostHoc/")
table(distshared_df$resptype,distshared_df$cond) # check

## new version of this with both models where distinct = genes ever only in heat or ever only in cold across ages & models
# and shared = the rest (if it's ever DE in both across ages and models)
# run PH_ageefftable.R first

make_distshared_new <- function(row_gene){
  hvb <- unname(unlist(row_gene["hvb"]))
  cvb <- unname(unlist(row_gene["cvb"]))
  if(hvb == "not_sig" | cvb == "not_sig"){
    return("distinct")
  }else{
    return("shared")
  }
}
distshared <- apply(age_dirs,1,make_distshared_new)
table(distshared) # seems plausible
age_dirs$distshared <- distshared
table(age_dirs$distshared,age_dirs$hvb) # seems plausible

# reshape so we count hot and cold independently of age and treatment (same gene can appear in heat and in cold) 
age_dirs_long <- pivot_longer(age_dirs, cols = cvb:hvb, names_to = "cond", values_to = "age_eff")
age_dirs_long <- filter(age_dirs_long, age_eff != "not_sig") # delete the not significant rows
# make labels
nr_cold <- table(age_dirs_long$cond)[1]
nr_hot <- table(age_dirs_long$cond)[2]
age_dirs_long$label <- NA
age_dirs_long[1,"label"] <- nr_cold
age_dirs_long[2,"label"] <- nr_hot
# flip distinct and shared and gray shared out because we are talking about distinct here
age_dirs_long$distshared <- factor(age_dirs_long$distshared, levels = c("shared","distinct"))
# plot
ggplot(age_dirs_long, aes(x = cond, fill = distshared,alpha = distshared))+
  geom_bar(position = "fill")+
  geom_text(aes(label = label, y = 1.05), alpha = 1)+
  scale_fill_manual(values = c("#00a08a","#f98400"))+
  scale_alpha_manual(values = c(0.4,1))+
  scale_x_discrete(labels = c("cvb" = "cold", "hvb" = "heat"))+
  guides(fill = guide_legend(reverse = TRUE, override.aes = list(alpha = c(1,0.4))), alpha = "none")+
  xlab("")+ylab("proportion of genes")+
  theme_pubr(base_size = 20)+
  theme(legend.title = element_blank())
# voilà

# now from the distinct, how many in cold and how many in heat
table(filter(age_dirs_long, distshared == "distinct")$cond)
23/386*100
363/386*100
# and how many of all cold genes are distinct, and how many of all heat genes are distinct?
table(age_dirs_long$cond,age_dirs_long$distshared)
23/224*100
363/564*100
# how many of all genes are consistently involved in the response to both heat and cold across development?
table(filter(age_dirs_long, distshared == "shared" & age_eff == "same_slope")$cond)


## bar distinct with up and down
# make df with effect sizes for distinct genes from both models (then I don't have to deal with ages here at all)
distgenes <- filter(age_dirs, distshared == "distinct")$gene
cat_eff_dist <- bind_rows(subset(cat_effects, gene %in% distgenes, select = c(gene,cvb,hvb)),
                      filter(cat_effects_noint, gene %in% distgenes & !(gene %in% cat_effects$gene)))
rownames(cat_eff_dist) <- NULL
# in there are some genes that are upregulated in one age but downregulated in another, thus some genes may appear twice in the bar plot
# but will make genes that are expressed in different intensities across ages but same direction just once
same_dir_hot <- filter(age_dirs, distshared == "distinct" & hvb == "same_direction")$gene
same_dir_cold <- filter(age_dirs, distshared == "distinct" & cvb == "same_direction")$gene # there are none
which(duplicated(cat_eff_dist$gene)) # 19
cat_eff_dist <- cat_eff_dist[-intersect(which(cat_eff_dist$gene %in% same_dir_hot),which(duplicated(cat_eff_dist$gene))),] # remove double genes with same direction but only one! should remove 16
rownames(cat_eff_dist) <- NULL
which(duplicated(cat_eff_dist$gene)) # so note in legend that 3 genes are in there double (one cold, two hot)
cat_eff_dist$colour <- ifelse(is.na(cat_eff_dist$hvb), "cold","hot")
make_dir_dist <- function(row){
  colour <- unname(unlist(row["colour"]))
  cvb <- as.numeric(row["cvb"])
  hvb <- as.numeric(row["hvb"])
  if(!(is.na(hvb)) & colour == "hot" & (hvb > 0)){
    return("down_regulated")
  }else if(!(is.na(cvb)) & colour == "cold" & (cvb > 0)){
    return("down_regulated")
  }else{return("up_regulated")}
}
dist_dir <- apply(cat_eff_dist, 1, make_dir_dist)
cat_eff_dist$dist_dir <- dist_dir
table(cat_eff_dist$dist_dir, cat_eff_dist$colour)
ggplot(filter(cat_eff_dist, resptype == "distinct"), aes(x = colour, fill = dist_dir))+
  geom_bar(position = position_dodge())+
  scale_fill_manual(values = c("#f8c75d","#5bbcd6"))+
  xlab("")+ylab("Number of genes")+
  theme_pubr(base_size = 20)+
  theme(legend.title = element_blank())
# ggsave("Bar_updowndist.pdf",path = "../../../figures/glmmSeq_PostHoc/")

## But I probably actually want something else, where genes can appear twice once in cold and once in heat
cat_eff_noint_long <- pivot_longer(dplyr::select(cat_effects_noint, c(gene,cvb,hvb,resptype)), cols = cvb:hvb, names_to = "condition", values_to = "effect")
cat_eff_noint_long <- filter(cat_eff_noint_long, !is.na(effect)) # remove NA rows
cat_eff_noint_long$direction <- ifelse(cat_eff_noint_long$effect < 0, "up-regulated","down-regulated") # add up or down
# plot
ggplot(cat_eff_noint_long, aes(x = condition, fill = direction))+
  geom_bar(position = position_dodge())+
  scale_fill_manual(values = c("#f8c75d","#5bbcd6"))+
  scale_x_discrete(labels = c("cvb" = "In cold", "hvb" = "In heat"))+
  xlab("")+ylab("Number of genes")+
  theme_pubr(base_size = 20)+
  theme(legend.title = element_blank())


### plot with hot dist and shared and cold dist and shared showing behaviour in ages
# need to run PH_ageefftable script first to get the dfs
make_colvec2 <- function(x){
  names(x) <- x
  if (x %in% c(genes_hvbNoInt,genes_hvb1w,genes_hvb8w) & x %in% c(genes_cvbNoInt,genes_cvb1w,genes_cvb8w)){
    return("hotandcold")
  }else if (x %in% c(genes_hvbNoInt,genes_hvb1w,genes_hvb8w)){
    return("hot")
  }else if (x %in% c(genes_cvbNoInt,genes_cvb1w,genes_cvb8w)){
    return("cold")
  }else{
    return("wtf")
  }
}
colvec2 <- sapply(cat_eff2$gene, function(x) make_colvec2(x))
table(colvec2) # idk could be right I guess?
cat_eff2$colour <- colvec2
# check
filter(cat_effects_noint, colour == "hotandcold")$gene %in% filter(cat_eff2, colour == "hotandcold")$gene # good
filter(cat_effects_noint, colour == "hot")$gene %in% filter(cat_eff2, colour == "hot")$gene # that's ok because some are hotandcold now because of age
filter(cat_effects_noint, colour == "cold")$gene %in% filter(cat_eff2, colour == "cold")$gene # that's ok because some are hotandcold now because of age
# add response type per condition
hotcoldgenes_age <- unname(unlist(cat_eff2[which(duplicated(cat_eff2$gene)),"gene"]))
hotcoldgenes_age %in% filter(cat_eff2, colour == "hotandcold")$gene # good
filter(cat_eff2, colour == "hotandcold")$gene %in% hotcoldgenes_age # two false, g18643 and g17406, both just1week
# some kind of mistake just change it
cat_eff2[which(cat_eff2$gene %in% c("g18643","g17406") & cat_eff2$condition == "hvb"), "colour"] <- "hot"
# now add response type but check numbers first so no overlap!
nshared <- nrow(cat_eff2[which(cat_eff2$gene %in% hotcoldgenes_age & cat_eff2$condition == "hvb"),])
nshared <- nrow(cat_eff2[which(cat_eff2$gene %in% hotcoldgenes_age & cat_eff2$condition == "cvb"),])
nhot <- nrow(cat_eff2[which(!(cat_eff2$gene %in% hotcoldgenes_age) & cat_eff2$condition == "hvb"),])
ncold <- nrow(cat_eff2[which(!(cat_eff2$gene %in% hotcoldgenes_age) & cat_eff2$condition == "cvb"),])
nshared + nhot + ncold # correct!

cat_eff2[which(cat_eff2$gene %in% hotcoldgenes_age & cat_eff2$condition == "hvb"),"condresp"] <- "hot_shared"
cat_eff2[which(cat_eff2$gene %in% hotcoldgenes_age & cat_eff2$condition == "cvb"),"condresp"] <- "cold_shared"
cat_eff2[which(!(cat_eff2$gene %in% hotcoldgenes_age) & cat_eff2$condition == "hvb"),"condresp"] <- "hot_distinct"
cat_eff2[which(!(cat_eff2$gene %in% hotcoldgenes_age) & cat_eff2$condition == "cvb"),"condresp"] <- "cold_distinct"
# add column for label (adjust numbers manually if needed)
labsdf <- data.frame(condresp = c("cold_distinct","hot_distinct","cold_shared","hot_shared"), label = c(paste0("n = ",ncold), paste0("n = ",nhot),
                                                                                                        paste0("n = ", nshared),paste0("n = ", nshared)))
cat_eff2 <- left_join(cat_eff2, labsdf)
# make nicer order
cat_eff2$ageeff <- factor(cat_eff2$age_eff, c("same_slope","same_direction","opposite_direction","just_1week","just_8week"))
cat_eff2$condresp <- factor(cat_eff2$condresp, c("cold_distinct","hot_distinct","cold_shared","hot_shared"))
# also save it for targeted fst analysis
# save(cat_eff2, file = "cat_eff_age_resp.RData")


ggplot(cat_eff2, aes(x = condresp, fill = ageeff))+
  geom_bar(position = "fill")+
  scale_fill_manual(values = c("khaki2","khaki4","brown4",unname(agecols[1]),unname(agecols[2])))+
  geom_text(aes(label = label, y = 1.05))+
  xlab("")+ylab("proportion of genes")+
  theme_pubr(base_size = 18)+
  theme(legend.title = element_blank())


### new bar plot with modulated by age yes/no and colour shared/distinct
# each gene only once
age_dirs$age_mod <- ifelse((age_dirs$hvb == "same_slope" | age_dirs$cvb == "same_slope"), "identical_response", "age_modulated")
table(age_dirs$age_mod) # adds up now, 458-96 = 362!
diffgenes <- cat_effects_noint[which(!(cat_effects_noint$gene %in% filter(age_dirs, age_mod == "identical_response")$gene)), "gene"] # I mean it is different models, I guess they just have slightly different results
# I am guessing they are mostly just different slope but let's check
diffgenes <- filter(age_dirs, gene %in% diffgenes) # no it's all types of things
table(diffgenes$cvb, diffgenes$hvb)
# anyway, add response type hotandcold if there's anything going on in both temps in any age!
cat_eff2$resptype <- ifelse(cat_eff2$colour == "hotandcold","shared","distinct")
# test that the same gene always has the same resptype, thus, distinct genes should only occur once!
unname(unlist(cat_eff2[which(duplicated(cat_eff2$gene)),"gene"])) %in% filter(cat_eff2, resptype == "shared")$gene # good
table(filter(cat_eff2, resptype == "shared")$gene %in% unname(unlist(cat_eff2[which(duplicated(cat_eff2$gene)),"gene"]))) # good!
table(filter(cat_eff2, resptype == "distinct")$gene %in% unname(unlist(cat_eff2[which(duplicated(cat_eff2$gene)),"gene"]))) # good!

age_dirs2 <- inner_join(age_dirs, subset(cat_eff2, select = c("gene","resptype")))
age_dirs2 <- distinct(age_dirs2) # it duplicates the rows, so remove again
table(age_dirs2$age_mod) #
table(age_dirs2$resptype)
# plot
ggplot(age_dirs2, aes(x = age_mod, fill = resptype))+
  geom_bar(position = position_dodge())+
  scale_fill_manual(values = c("#f98400", "#00a08a"))+
  ylab("Number of genes")+ xlab("")+
  theme_pubr(base_size = 20)+
  theme(legend.title = element_blank())
# plot other way around
ggplot(age_dirs2, aes(x = resptype, fill = age_mod))+
  geom_bar(position = position_dodge())+
  scale_fill_manual(values = c("#f98400", "#00a08a"))+
  ylab("Number of genes")+ xlab("")+
  theme_pubr(base_size = 20)+
  theme(legend.title = element_blank())


# NEW: hot, cold and hotandcold separately
# make facet labeller
colour_labs <- c("Distinct: hot", "Distinct: cold", "Shared")
names(colour_labs) <- c("hot","cold","hotandcold")
ggplot(filter(cat_eff2, condition == "cvb" & ageeff != "same_slope"), aes (x = ageeff, fill = colour))+
  facet_wrap(facets = vars(colour), nrow = 2, labeller = labeller(colour = colour_labs))+
  geom_bar()+
  xlab("")+ylab("number of genes")+
  theme_pubr(base_size = 15)+
  scale_fill_manual(values = unname(tempcols[c(2,1)]))+
  scale_x_discrete(labels = c( "same direction","opposite direction","just 1-week-old", "just 8-week-old"))+
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
ggplot(filter(cat_eff2, condition == "hvb" & ageeff != "same_slope"), aes (x = ageeff, fill = colour))+
  facet_wrap(facets = vars(colour), nrow = 2, labeller = labeller(colour = colour_labs))+
  geom_bar()+
  xlab("")+ylab("number of genes")+
  scale_fill_manual(values = unname(tempcols[c(3,1)]))+
  scale_x_discrete(labels = c( "same direction","opposite direction","just 1-week-old", "just 8-week-old"))+
  theme_pubr(base_size = 15)+
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

# save this info for targeted fst analysis



### plot of all genes for paper
# what happens when we're a bit stricter about linearity?
fitfull_gr <- glmmQvals(fitfull_gr, cutoff = 0.001)
fit_lin_noint <- glmmQvals(fit_lin_noint, cutoff = 0.005)
statsct_agesig005 <- as.data.frame(statsct_agesig[1:219,])
statsct_agesig001 <- as.data.frame(statsct_agesig[1:205,])
statslin_noint_sig005 <- as.data.frame(statslin_noint_sig[1:257,])
statslin_noint_sig001 <- as.data.frame(statslin_noint_sig[1:206,])

lin_genes005 <- lin_genes[lin_genes %in% c(statsct_agesig005$gene, statslin_noint_sig005$gene)]
lin_genes001 <- lin_genes[lin_genes %in% c(statsct_agesig001$gene, statslin_noint_sig001$gene)]

inbyage <- function(x){
  names(x) <- x
  if (x %in% statscat_withint_sig$gene){ # order is important! first cat then lin, and first age dependent
    return("age_dependent")
  }else if (x %in% statsFull_cat_sig$gene){
    return("ages_same")
  }else if (x %in% statsct_agesig$gene){
    return("age_dependent")
  }else if(x %in% statslin_noint_sig$gene){
    return("ages_same")
  }else{
    return("unclear")
  }
}

agedep <- sapply(cat_effects_all$gene, function(x) inbyage(x))
table(agedep) # looking good
cat_effects_all$agedep <- agedep
# also need colour thing

colvecall <- sapply(cat_effects_all$gene, function(x) make_colvec2(x))
table(colvecall)
cat_effects_all$colour <- colvecall

cat_eff_wide <- pivot_wider(dplyr::select(cat_effects, gene:hvb), names_from = age, values_from = c(cvb,hvb))
cat_eff_wide[is.na(cat_eff_wide$cvb_week1),"cvb_week1"] <- 0
cat_eff_wide[is.na(cat_eff_wide$cvb_week8),"cvb_week8"] <- 0
cat_eff_wide[is.na(cat_eff_wide$hvb_week1),"hvb_week1"] <- 0
cat_eff_wide[is.na(cat_eff_wide$hvb_week8),"hvb_week8"] <- 0
# add info where expressed

cat_eff_wide <- left_join(cat_eff_wide, subset(cat_effects_all, select = c(gene, colour)))

## could also include no interaction with equal effects in hot and cold? But only maybe
cat_effects_noint_ageish <- dplyr::rename(cat_effects_noint, cvb_week1 = cvb, hvb_week1 = hvb)
cat_effects_noint_ageish$cvb_week8 <- cat_effects_noint_ageish$cvb_week1
cat_effects_noint_ageish$hvb_week8 <- cat_effects_noint_ageish$hvb_week1

# cat_eff_wide <- rbind(cat_eff_wide, subset(cat_effects_noint_ageish, !(gene %in% cat_eff_wide$gene), select = c(gene, cvb_week1,cvb_week8,hvb_week1,hvb_week8,colour)))

# turn signs around so it fits with the text
cat_eff_wide$cvb_week1 <- cat_eff_wide$cvb_week1*-1
cat_eff_wide$cvb_week8 <- cat_eff_wide$cvb_week8*-1
cat_eff_wide$hvb_week1 <- cat_eff_wide$hvb_week1*-1
cat_eff_wide$hvb_week8 <- cat_eff_wide$hvb_week8*-1

ggplot(filter(cat_eff_wide, !(hvb_week1 == 0 & hvb_week8 == 0)), aes(x = hvb_week1, y = hvb_week8, colour = colour))+
  geom_jitter(height=0.02, width = 0.02, alpha = 0.6)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_colour_manual(values = unname(tempcols[c(3,1)]), labels = c("Distinct: hot", "Shared"))+
  scale_x_continuous(name="Expression in heat compared to benign \n (1-week-olds)", breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.85,1.85))+
  scale_y_continuous(name="Expression in heat compared to benign \n (8-week-olds)", breaks = c(-1,-0.5,0,0.5,1), limits = c(-1.2,1.2))+
  theme_minimal(base_size = 20)+
  theme(legend.title = element_blank(), legend.position = "top")

ggplot(filter(cat_eff_wide, !(cvb_week1 == 0 & cvb_week8 == 0)), aes(x = cvb_week1, y = cvb_week8, colour = colour))+
  geom_jitter(height=0.02, width = 0.02, alpha = 0.6)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_colour_manual(values = unname(tempcols[c(2,1)]), labels = c("Distinct: cold", "Shared"))+
  scale_x_continuous(name="Expression in heat compared to benign \n (1-week-olds)", breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.85,1.85))+
  scale_y_continuous(name="Expression in heat compared to benign \n (8-week-olds)", breaks = c(-1,-0.5,0,0.5,1), limits = c(-1.2,1.2))+
  theme_minimal(base_size = 20)+
  theme(legend.title = element_blank(), legend.position = "top") # quite a lot different directions!

# in one plot difficult because genes are there twice
# it looks like there are proportionally more genes that change direction with age in cold than in heat is that right?
length(which(cat_effects_both_all$cvb_week8 < 0 & cat_effects_both_all$cvb_week1 > 0)) # 26
length(which(cat_effects_both_all$cvb_week8 > 0 & cat_effects_both_all$cvb_week1 < 0)) # 10
length(which(cat_effects_both_all$cvb_week8 > 0 & cat_effects_both_all$cvb_week1 > 0)) # 28
length(which(cat_effects_both_all$cvb_week8 < 0 & cat_effects_both_all$cvb_week1 < 0)) # 10
# that does not add up there should be 114 in total according to below
# ah no it's because genes that are in both hot and cold are classified as both sig if they are sig either in hot or cold!


ggplot(cat_effects_all, aes(x = hvb, y = cvb, fill = colour))+
  geom_point(shape = 21, alpha = 0.7)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_fill_manual(values = c(linear = "green3",hot=tempcols[["hot"]],hotandcold=tempcols[["benign"]],cold=tempcols[["cold"]]))+
  scale_x_continuous(name="Gene expression in benign compared to heat", breaks = c(-1,-0.5,0,0.5,1), limits = c(-1,1.7))+
  scale_y_continuous(name="Gene expression in benign compated to cold", breaks = c(-0.5,0,0.5,1), limits = c(-0.55,0.9))+
  theme_minimal(base_size = 20)
# plots with separate hot cold and both
ggplot(filter(cat_effects_all, colour %in% "hot"), aes(x = hvb, y = cvb, fill = colour))+
  geom_point(shape = 21, alpha = 0.7)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_fill_manual(values = c(linear = "green3",hot=tempcols[["hot"]],hotandcold=tempcols[["benign"]],cold=tempcols[["cold"]]))+
  scale_x_continuous(name="Gene expression in benign compared to heat", breaks = c(-1,-0.5,0,0.5,1), limits = c(-1,1.7))+
  scale_y_continuous(name="Gene expression in benign compated to cold", breaks = c(-0.5,0,0.5,1), limits = c(-0.55,0.9))+
  theme_minimal(base_size = 20)
# kind of defies the message that a lot are opposite because of all the cat heat ones
table(filter(cat_effects_all, colour == "linear")$hvb < 0 & filter(cat_effects_all, colour == "linear")$cvb > 0) # 118 in upper left
table(filter(cat_effects_all, colour == "linear")$hvb > 0 & filter(cat_effects_all, colour == "linear")$cvb < 0) # 62 in lower right
table(filter(cat_effects_all, colour == "linear")$hvb > 0 & filter(cat_effects_all, colour == "linear")$cvb > 0) # 120 in upper right
table(filter(cat_effects_all, colour == "linear")$hvb < 0 & filter(cat_effects_all, colour == "linear")$cvb < 0) # 25 in lower left
# upper right and lower left should not exist for linear!!
# but maybe that's because I am including genes that have different slopes depending on age which will muddle it all
# so now same plot with only genes from the models without intercept
ggplot(filter(cat_effects_all, gene %in% c(statslin_noint_sig$gene,statsFull_cat_sig$gene)), aes(x = hvb, y = cvb, fill = colour))+
  geom_point(shape = 21, alpha = 0.7)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_fill_manual(values = c(linear = "green3",hot=tempcols[["hot"]],hotandcold=tempcols[["benign"]],cold=tempcols[["cold"]],unclear = "gray"))+
  scale_x_continuous(name="Gene expression in benign minus heat", breaks = c(-1,-0.5,0,0.5,1), limits = c(-1,1.7))+
  scale_y_continuous(name="Gene expression in benign minus cold", breaks = c(-0.5,0,0.5,1), limits = c(-0.55,0.9))+
  theme_minimal(base_size = 20)
# sadly that did not solve the problem
# let's just look at linear ones
ggplot(filter(cat_effects_all, colour %in% c("linear")), aes(x = hvb, y = cvb, fill = agedep))+
  geom_point(shape = 21,alpha = 0.7)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_fill_manual(values = c("green2","green4"))+
  scale_x_continuous(name="Gene expression in benign compared to heat", breaks = c(-1,-0.5,0,0.5,1), limits = c(-1,1.7))+
  scale_y_continuous(name="Gene expression in benign compated to cold", breaks = c(-0.5,0,0.5,1), limits = c(-0.55,0.9))+
  ggtitle("linear genes at p < 0.001")+
  theme_minimal(base_size = 20)
# take a look at the genes that don't behave as expected
weirdg_effs <- filter(cat_effects_all, agedep == "ages_same" & colour == "linear" & hvb > 0 & cvb > 0)
rownames(weirdg_effs) <- NULL
summary(weirdg_effs) # most of them are not significant, esp not for cold
# look at some
summary(GeneFitsLinNoint$g3872) # negative slope although clearly up in heat according to cat_effects_all table
summary(GeneFits_all$g3872) # here benign estimate is highest and cold close, and heat lowest -> in line with what linear models says more or less
EMres_all$g3872 # again bgn and cold higher than hot
EMpairs_all$g3872 # but now here benign - hot is positive! What happened? I read it wrong, of course b-h is pos when bgn is higher duuhhh
# is the slope always negative?
table((weirdg_effs$hvb - weirdg_effs$cvb)>0) # ok for very few of them the difference b-c is bigger than b-h, so their slope should be positive!
weirdg_effs$hvb - weirdg_effs$cvb # eg nr 31 aka g16950
summary(GeneFitsLinNoint$g14301) # yep positive slope, so one cannot even count on that
summary(GeneFits_all$g14301) # hot and bgn basically the same
EMres_all$g14301 # again bgn and cold higher than hot
EMpairs_all$g14301 # but now here benign - hot is positive! What happened? I read it wrong, of course b-h is pos when bgn is higher duuhhh
# for almost all of these the difference from benign is not significant (at 0.01) so maybe we should just put them to 0 to exemplify that
cat_effects_all0 <- cat_effects_all
cat_effects_all0[cat_effects_all0$p_hvb > 0.01, "hvb"] <- 0
cat_effects_all0[cat_effects_all0$p_cvb > 0.01, "cvb"] <- 0
ggplot(cat_effects_all0, aes(x = hvb, y = cvb, fill = colour))+
  geom_jitter(shape = 21, alpha = 0.7, width = 0.03, height = 0.03)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_fill_manual(values = c(linear = "green3",hot=tempcols[["hot"]],hotandcold=tempcols[["benign"]],cold=tempcols[["cold"]]))+
  scale_x_continuous(name="Gene expression in heat compared to benign", breaks = c(-1,-0.5,0,0.5,1))+
  scale_y_continuous(name="Gene expression in cold compared to benign", breaks = c(-0.5,0,0.5,1))+
  theme_minimal(base_size = 20)
# isn't really helpful maybe
linlike_effs <- filter(cat_effects_all, agedep == "ages_same" & colour == "linear" & hvb < 0 & cvb > 0)
rownames(linlike_effs) <- NULL
summary(linlike_effs) # they are all not significant
# check that they have positive slope as they should
summary(GeneFitsLinNoint$g10463)$coefficients[2,1] 
sloopes <- sapply(GeneFitsLinNoint, function(fit) summary(fit)$coefficients[2,1])
sloopes <- sloopes[names(sloopes) %in% linlike_effs$gene]
summary(sloopes) # yep all positive, and also tested the opposite good all negative




# I don't know how I was thinking above, but from this plot it's clear which of the shared cat ones are same and opposite direction
# how many?
table(filter(cat_effects_all, colour == "hotandcold")$hvb < 0 & filter(cat_effects_all, colour == "hotandcold")$cvb > 0) # 29 in upper left
table(filter(cat_effects_all, colour == "hotandcold")$hvb > 0 & filter(cat_effects_all, colour == "hotandcold")$cvb < 0) # 43 in lower right
# take these out as well for GO just for completeness sake
catshared_hotupcolddown <- cat_effects_all[which(cat_effects_all$colour == "hotandcold" & (cat_effects_all$hvb >0 & cat_effects_all$cvb < 0)),"gene"]
catshared_hotupcoldup <- cat_effects_all[which(cat_effects_all$colour == "hotandcold" & (cat_effects_all$hvb >0 & cat_effects_all$cvb > 0)),"gene"]
catshared_hotdowncolddown <- cat_effects_all[which(cat_effects_all$colour == "hotandcold" & (cat_effects_all$hvb < 0 & cat_effects_all$cvb < 0)),"gene"]
catshared_hotdowncoldup <- cat_effects_all[which(cat_effects_all$colour == "hotandcold" & (cat_effects_all$hvb <0 & cat_effects_all$cvb > 0)),"gene"]



## what are hot doing in cold and vice versa?
cat_effects_noint[is.na(cat_effects_noint$cvb),"cvb"] <- 0
cat_effects_noint[is.na(cat_effects_noint$hvb),"hvb"] <- 0
ggplot(cat_effects_noint, aes(x = hvb, y = cvb))+
  geom_point(shape = 1)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_x_continuous(name="Gene expression in heat compared to benign", breaks = c(-1,-0.5,0,0.5,1))+
  scale_y_continuous(name="Gene expression in cold compared to benign", breaks = c(-0.5,0,0.5,1))+
  theme_minimal(base_size = 15)
# some genes are in there twice that's not great


## check how many same direction and how many up down in each age (corresponds with fold change plot?)
# linear
lin_eff_wide <- pivot_wider(dplyr::select(linear_effects, c(gene, age, slope)), names_from = age, values_from = slope)
lin_eff_wide <- dplyr::rename(lin_eff_wide, week1 = "1w", week8 = "8w")
table(lin_eff_wide$week8<0,lin_eff_wide$week1>0) # (11+56)/(11+56+8+39)*100 = 58.8 % are same direction
lin_eff_wide[which(is.na(lin_eff_wide$week1)),"week1"] <- 0 # for the sake of plotting
lin_eff_wide[which(is.na(lin_eff_wide$week8)),"week8"] <- 0
ggplot(lin_eff_wide, aes(x = week1, y = week8))+
  geom_point()+
  theme_minimal()
# resembles other plot in some ways! But less genes
