## chamber experiment chick cloacal temperature ##

#### load packages ####
library("readxl")
library("tidyverse")
library("ggpubr")
library("RColorBrewer")
library("ggExtra")
library(lme4)
library(viridis)

#### preliminaries ####
setwd("/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/data/working")
figurepath="/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/figures/Physiology"
### colours
tempcols = c(benign = "#6800a4", cold = "#377EB8", hot = "#E41A1C") # or should I do purple for benign? #4DAF4A
agecols = c("#1B9E77" ,"#D95F02")
sexcols = c("#7570B3", "#E7298A")
shadecols = c("orange","gray")
shadecols2 = c("black","gray")

# in here there's everything, don't need the other things any more
load("ChamberTempSI_new.RData")
# relevel treatment so benign in the middle
sampleInfo$treatment <- factor(sampleInfo$treatment, levels = c("cold","benign","hot")) # maybe this messes up some things watch out!

# one row per chick
AllChicks <- pivot_wider(datTall, id_cols = !c(Date, SampleNr, Samplingtime, chamberT), names_from = treatment, values_from = cloacat) # should have 34 rows
SIwide <- pivot_wider(sampleInfo, id_cols = !c(Date, SampleNr, Samplingtime, chamberT), names_from = treatment, values_from = cloacat) # should have 20 rows

# often I called it datRNAseq, just the same as sampleInfo I think
datRNAseq <- sampleInfo


#### plotting cloacal temp data ####

ggplot(datTalive, aes(x=treatment, y=cloacat, fill= Age))+
  geom_boxplot()+
  geom_jitter(position = position_jitterdodge(jitter.height = 0), aes(shape = Age))+
  facet_wrap(vars(as.factor(Date)))+
  scale_fill_manual(values=sexcols)+
  theme_pubr()
  
#per group: group 1 is colder in benign than cold (that was day 2), no consistent age difference
#per date: first day best differenciated
# using only survivors cold variation on first day 1 week and benign variation on 2nd day 8 week smaller
# don't see any sex differences whatsoever

# plot only with the ones that survived
ggplot(datTalive, aes(x=treatment, y=cloacat, colour=Age, group=chickno))+
  #geom_boxplot()+
  facet_wrap(vars(group))+
  geom_point()+
  geom_line()+
  theme_pubr()

# include info wether chick died the next day
datTall$deadnextday =c(rep(NA))
for (i in 1:nrow(datTall)){
    if(is.na(datTall$mort_dat[i])){ 
    }else if(datTall$Date[i] == as.POSIXct("2023-10-25", tz = "UTC") & datTall$mort_dat[i] == as.POSIXct("2023-10-26", tz = "UTC") |
       (datTall$Date[i] == as.POSIXct("2023-10-26", tz = "UTC") & datTall$mort_dat[i] == as.POSIXct("2023-10-27", tz = "UTC"))){
      datTall$deadnextday[i] = "yes"
    }else{
      datTall$deadnextday[i] = "no_or_NA"
  }
}
datTall$deadnextday <- as.factor(datTall$deadnextday)
ggplot(datTall, aes(x=treatment, y=cloacat, fill= Age))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(position = position_jitterdodge(jitter.height = 0), aes(color = deadnextday))+
  scale_fill_manual(values=agecols)+
  scale_color_manual(values = c("black","red"))+
  theme_pubr()


### how hot was it actually? eg mean of last two hours or so?


ggplot(sampleInfo, aes(x=chamberT, y=cloacat, color=treatment, shape = Age))+
  geom_jitter(height = 0.1, width = 0.1, size = 3)+
  scale_x_continuous(limits = c(10,47))+
  scale_color_manual(values= tempcols)+
  xlab("Ambient temperature [°C]")+
  ylab("Cloacal temperature [°C]")+
  theme_pubr(base_size = 25)
# ggsave("cloacaTvchamberT_120to-30_SIdata.pdf", path = figurepath)



# import temp data
hdrs = c("Reading","Date.Time","temp","humidity","dewpoint","TrueTime","ManTime")
tempdata = vector("list", length = 9)
names(tempdata) = c("hot1", "hot2", "hot3","cold1","cold2","cold3","benign1","benign2","benign3")
tempdata$hot1 = read_xlsx("Temperatures/Day1/1HotDay1.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$hot3 = read_xlsx("Temperatures/Day3/1HotDay3.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$hot2 = read_xlsx("Temperatures/Day2/1HotDay2.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$cold1 = read_xlsx("Temperatures/Day1/2ColdDay1.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$cold2 = read_xlsx("Temperatures/Day2/2ColdDay2.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$cold3 = read_xlsx("Temperatures/Day3/2ColdDay3.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$benign1 = read_xlsx("Temperatures/Day1/3BenignDay1.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$benign2 = read_xlsx("Temperatures/Day2/3BenignDay2.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$benign3 = read_xlsx("Temperatures/Day3/3BenignDay3.xlsx", skip = 11, col_names = hdrs)[,1:7]


## plot all temps on all days together so I can really see
pday1=ggplot()+
  geom_line(data=tempdata$hot1, aes(x=ManTime,y=temp, color = tempcols[2]))+
  geom_line(data=tempdata$cold1,aes(x=ManTime,y=temp, color = tempcols[3]))+
  geom_line(data=tempdata$benign1,aes(x=ManTime,y=temp, color = tempcols[1]))+
  scale_x_datetime(limits = c(tempdata$hot1$ManTime[1],tempdata$benign3$ManTime[431]), date_breaks = "1 hour",
                   date_labels = "%H:%M", name = "Time")+
  scale_y_continuous(limits = c(10,48))+
  scale_color_discrete(name="Chamber", labels = c("Hot","Benign","Cold"))+
  theme_pubr()
pday2=ggplot()+
  geom_line(data=tempdata$hot2, aes(x=ManTime,y=temp, color = tempcols[2]))+
  geom_line(data=tempdata$cold2,aes(x=ManTime,y=temp, color = tempcols[3]))+
  geom_line(data=tempdata$benign2,aes(x=ManTime,y=temp, color = tempcols[1]))+
  scale_x_datetime(limits = c(tempdata$hot1$ManTime[1],tempdata$benign3$ManTime[431]), date_breaks = "1 hour",
                   date_labels = "%H:%M", name = "Time")+
  scale_y_continuous(limits = c(10,48))+
  scale_color_discrete(name="Chamber", labels = c("Hot","Benign","Cold"))+
  theme_pubr()
pday3=ggplot()+
  geom_line(data=tempdata$hot3, aes(x=ManTime,y=temp, color = tempcols[2]))+
  geom_line(data=tempdata$cold3,aes(x=ManTime,y=temp, color = tempcols[3]))+
  geom_line(data=tempdata$benign3,aes(x=ManTime,y=temp, color = tempcols[1]))+
  scale_x_datetime(limits = c(tempdata$hot1$ManTime[1],tempdata$benign3$ManTime[431]), date_breaks = "1 hour",
                   date_labels = "%H:%M", name = "Time")+
  scale_y_continuous(limits = c(10,48))+
  scale_color_discrete(name="Chamber", labels = c("Hot","Benign","Cold"))+
  theme_pubr()
ggarrange(pday1,pday2,pday3, nrow=3)

ggplot()+
  geom_line(data=tempdata$hot1, aes(x=ManTime,y=temp, linetype = "solid"),color = tempcols[3])+
  geom_line(data=tempdata$cold1,aes(x=ManTime,y=temp, linetype = "solid"), color = tempcols[2])+
  geom_line(data=tempdata$benign1,aes(x=ManTime,y=temp, linetype = "solid"), color = tempcols[1])+
  geom_line(data=tempdata$hot2, aes(x=ManTime,y=temp, linetype = "dashed"), color = tempcols[3])+
  geom_line(data=tempdata$cold2,aes(x=ManTime,y=temp, linetype = "dashed"), color = tempcols[2])+
  geom_line(data=tempdata$benign2,aes(x=ManTime,y=temp, linetype = "dashed"), color = tempcols[1])+
  geom_line(data=tempdata$hot3, aes(x=ManTime,y=temp, linetype = "twodash"), color = tempcols[3])+
  geom_line(data=tempdata$cold3,aes(x=ManTime,y=temp, linetype = "twodash"), color = tempcols[2])+
  geom_line(data=tempdata$benign3,aes(x=ManTime,y=temp, linetype = "twodash"), color = tempcols[1])+
  scale_x_datetime(limits = c(tempdata$hot1$ManTime[1],tempdata$benign3$ManTime[431]), date_breaks = "1 hour",
                   date_labels = "%H:%M", name = "Time")+
  scale_y_continuous(limits = c(10,48))+
  scale_color_discrete(name="Chamber", labels = c("Hot","Benign","Cold"))+
  scale_linetype_discrete(name = "Day", labels = c("1","2","3"))+
  ylab("Air temperature [°C]")+
  theme_pubr(base_size = 20)
#ggsave("Chamber_temp_all.pdf", path = figurepath)

## graph with average temps of the 3 days, not to replace but to put into paper
# Should quickly show what temps chicks were exposed to over time, more to illustrate than show data (that's for SI)
# -> average?
len <- length(tempdata$hot1$temp[8:367]) # from 09:16 to 15:15 -> 6 hours
tempdatframe <- data.frame(treatment = c(rep("hot",len*3),rep("cold",len*3),rep("benign",len*3)),
                           day = rep(c(rep(1,len),rep(2,len),rep(3,len)),3),
                           temp = c(tempdata$hot1$temp[32:(31+len)],tempdata$hot2$temp[8:(7+len)],tempdata$hot3$temp[18:(17+len)],
                                    tempdata$cold1$temp[32:(31+len)],tempdata$cold2$temp[8:(7+len)],tempdata$cold3$temp[18:(17+len)],
                                    tempdata$benign1$temp[32:(31+len)],tempdata$benign2$temp[8:(7+len)],tempdata$benign3$temp[18:(17+len)]),
                           ManTime = rep(tempdata$hot3$ManTime[18:(17+len)],9))
meantemp <- dplyr::summarise(group_by(tempdatframe, ManTime, treatment), mean_temp=mean(temp))

ggplot()+
  geom_line(data=filter(meantemp, treatment == "hot"), aes(x=ManTime,y=mean_temp, color = tempcols[2]))+
  geom_line(data=filter(meantemp, treatment == "cold"), aes(x=ManTime,y=mean_temp, color = tempcols[3]))+
  geom_line(data=filter(meantemp, treatment == "benign"), aes(x=ManTime,y=mean_temp, color = tempcols[1]))+
  scale_x_datetime(date_breaks = "1 hour", date_labels = "%H:%M", name = "Time")+
  scale_y_continuous(limits = c(10,48))+
  scale_color_discrete(name="Chamber", labels = c("Hot","Benign","Cold"))+
  theme_pubr()
#ggsave("Mean_chamber_temp.pdf", path = figurepath)


# which time of exposure correlates best with cloacal temp? Doesn't matter much for cold and benign, but in hot when plotting ~2.5
# hour mean cloacal temp seems to increase exponentially, while when plotting the last hour it looks more like its plateauing.
# btw the variation within days and treatments looks very nice

# cold didn't have a big impact

#### Cloaca data statistics ------------------------------------------------------------------------
### Only with sequenced chicks
# Question: Does treatment significantly influence cloaca temp, and does age significantly influence cloaca temp?
# response: Cloaca temp, assume normal distribution
# regression: cloacat ~ treatment*Age + (1|group) + (1|chickno)
# and then two-way repeated measures ANOVA -> is it different?
# and then post-hoc test -> where is the difference?
library(lme4)
library(GGally)
library(DHARMa)
library(MuMIn)
library(lmerTest)
library(rstatix)
library(emmeans)
library(visreg)

SIwide = pivot_wider(dplyr::select(sampleInfo,c(treatment,chickno, cloacat,group,Age,sex,Mass,massScaled)),
                     names_from = treatment, values_from = cloacat) # wide with one row per chick
SIwide = mutate(SIwide,TchangeHot = hot-benign)
SIwide = mutate(SIwide,TchangeCold = cold-benign)
SIwide <- mutate(SIwide, TchangeAv = (abs(TchangeHot)+abs(TchangeCold))/2)
SIwide$chickno <- as.factor(SIwide$chickno)
# add to sampleInfo
sampleInfo <- left_join(sampleInfo, select(SIwide, c(chickno, TchangeHot, TchangeCold, TchangeAv)))


datRNAseq <- sampleInfo # use only chicks that I am actually using

m1 <- lmer(cloacat ~ treatment*Age + treatment*sex + (1|group) + (1|chickno) + (1|Date), data = datRNAseq)
simulationOutput <- simulateResiduals(fittedModel = m1, n = 500) # simulate data from our model n times
# This method checks if your model is useful by seeing if it can even produce data that looks like the data you used to fit your model (i.e. your observations). 
plot(simulationOutput) # looks ok

options(na.action = "na.fail") # needed for dredge() function to prevent illegal model comparisons
dredgeOut<-dredge(m1, extra = "R^2") # fit and compare a model set representing all possible predictor combinations
# gives singularity warning
dredgeOut # seems like only treatment has an effect, age is at the border (but w/o interaction)

datRNAseq8 <- filter(datRNAseq, Age == "8week")

# so lets do it with only treatment, m2 is the right one, or rather m5
m2 <- lmerTest::lmer(cloacat ~ treatment + (1|chickno) + (1|group)+ (1|Date), data = datRNAseq)
m3 <- lmerTest::lmer(cloacat ~ treatment + (1|chickno) + (0 + treatment|group)+ (0 + treatment|Date), data = datRNAseq)
m4 <- lmerTest::lmer(cloacat ~ treatment + (1|chickno) + (1 + treatment|group)+ (1 + treatment|Date), data = datRNAseq)
m5 <- lmerTest::lmer(cloacat ~ treatment + (1|chickno) + (1|Date), data = datRNAseq)
m6 <- lmerTest::lmer(cloacat ~ treatment + Age + sex + (1|chickno) + (1|Date), data = datRNAseq)
simulationOutput <- simulateResiduals(fittedModel = m6, n = 500) # simulate data from our model n times
# This method checks if your model is useful by seeing if it can even produce data that looks like the data you used to fit your model (i.e. your observations). 
plot(simulationOutput) # looks ok 
AIC(m2,m3,m4,m5,m6) # m2 and m5 have almost the same AIC, lower than m3 and m4 and m6

summ5 <- summary(m5)# here there are already p values saying all treatments are significantly different from each other
summ_co <- summ5$coefficients # extract model coefficients
# extract variance explained by random effects
v2 <- VarCorr(m5)
v2df <- as.data.frame(v2)
v2df <- subset(v2df, select = -var2)
v2df <- reshape::rename(v2df, c(grp = "Random effect", var1 = "Var1", vcov = "Variance", sdcor = "Std. Dev."))
# chickno explains more variance than Date but less than residual variance

# anova
aov_res1 <- anova(m1)
aov_res1 # shows again that sex and age are not significant but treatment is
aov_res2 <- anova(m2)
aov_res2
aov_res5 <- anova(m5)
aov_res5
aov_res6 <- anova(m6)
aov_res6
ranova(m5)
drop1(m5)
step(m4)
get_model(step(m4)) # here it says m5 is the best
# but with random effects the important thing is whether they explain any variance, which group does, so I'll leave it in (also no singularity error)
# or maybe not after all because it should be nested in chickno

emmeans(m5, pairwise ~ treatment, type = "response") # this is what I wanted !!!!
plot(emmeans(m5, pairwise ~ treatment, type = "response"), comparisons = TRUE)+
  coord_flip()+
  theme_pubr()

# visualise
visreg(m5, scale = "response", xvar = "treatment", gg = TRUE)+
  theme_pubr()
  
ggplot(sampleInfo, aes(x=treatment, y=cloacat, fill=Age))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(shape = sex, group = interaction(treatment, Age)), size = 2, 
              position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.75))+
  geom_signif(comparisons = list(c("hot","cold")),map_signif_level=TRUE)+
  geom_signif(comparisons = list(c("hot","benign")),map_signif_level=TRUE, y_position = c(42,1))+
  ylab("cloaca temperature [°C]")+
  scale_fill_manual(values=agecols)+
  theme_pubr(base_size = 22)
# ggsave("cloacatVtreatment_box_stats.pdf", path = figurepath)

# no age and no sex for beginning of paper as not important yet
ggplot(sampleInfo, aes(x=treatment, y=cloacat, colour=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  geom_signif(comparisons = list(c("hot","cold")),map_signif_level=TRUE,color = "black")+
  geom_signif(comparisons = list(c("hot","benign")),map_signif_level=TRUE, y_position = c(42,1),color = "black")+
  ylab("cloaca temperature [°C]")+xlab("")+
  scale_colour_manual(values=tempcols)+
  theme_pubr(base_size = 22)
#ggsave("cloacatVtreatment_box_stats2.pdf", path = figurepath)

## for presentation hot and cold graphs separately
ggplot(filter(sampleInfo, treatment != "hot"), aes(x=treatment, y=cloacat, fill=Age))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(shape = sex, group = interaction(treatment, Age)), size = 2, 
              position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.75))+
  #geom_signif(comparisons = list(c("hot","benign")),map_signif_level=TRUE, y_position = c(42,1))+
  ylab("cloaca temperature [°C]")+
  scale_fill_manual(values=agecols)+
  theme_pubr(base_size = 22)
ggsave("cloacatVtreatment_box_justcold.pdf", path = figurepath)

# or even without age
ggplot(filter(sampleInfo, treatment != "hot"), aes(x=treatment, y=cloacat, fill = treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(shape = sex), size = 3, 
              position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.75))+
  #geom_signif(comparisons = list(c("hot","benign")),map_signif_level=TRUE, y_position = c(42,1))+
  scale_fill_manual(values = tempcols)+
  ylab("cloaca temperature [°C]")+
  theme_pubr(base_size = 25)
ggsave("cloacatVtreatment_box_coldnoage.pdf", path = figurepath)

### anova for repeated measures (without regression)
# Are residuals normally distributed at each time point?
ggqqplot(datRNAseq, "cloacat") # hm not sure
shapiro.test(filter(datRNAseq, as.factor(Date) == "2023-10-25")$cloacat) # p-value =  0.08102 -> normal
shapiro.test(filter(datRNAseq, as.factor(Date) == "2023-10-26")$cloacat) # p-value =  0.04358 -> njäääh
shapiro.test(filter(datRNAseq, as.factor(Date) == "2023-10-27")$cloacat) # p-value =  0.3642 -> normal
datRNAseq %>%
  group_by(treatment, Date) %>%
  shapiro_test(cloacat)
# there is no alternative anyways
aov1 <- aov(cloacat ~ treatment + Age + sex + Error(chickno), data = datRNAseq)
summary(aov1) # age and sex are not significant (p: 0.202 and 0.392) but treatment is (<2e-16)
aov2 <- aov(cloacat ~ treatment + Age + sex + Error(chickno,Date), data = datRNAseq)
summary(aov2) # same result

# post hoc test: which treatments? -> run one-way model at each level of variable


### Temp deviation in heat vs cold
# significant?
SIwide <- filter(SIwide, !is.na(TchangeCold) & !is.na(TchangeHot))
SIwide$massScaled <- as.numeric(SIwide$massScaled)
fit <- lm(TchangeCold ~ TchangeHot, data = SIwide)
summary(fit) # significant
cor.test(SIwide$TchangeCold,SIwide$massScaled, method = "spearman") #nothing
# account for body mass
fit_bm <- lm(TchangeCold ~ TchangeHot + massScaled, data = SIwide) # can only be included as fixed effect obv
fit_bm_int <- lm(TchangeCold ~ TchangeHot*massScaled, data = SIwide) # interaction conceivable
# test other models that might fit even better
fit_bm_as <- lm(TchangeCold ~ TchangeHot + massScaled + Age + sex, data = SIwide) # can only be included as fixed effect obv
summary(fit_bm)
summary(fit_bm_int) # mass or interaction are not significant
summary(fit_bm_as)
simulationOutput <- simulateResiduals(fittedModel = fit, n = 500) # simulate data from our model n times
plot(simulationOutput, asFactor=FALSE) # looking OK
# so what is best model now?
AIC(fit,fit_bm,fit_bm_int,fit_bm_as) # with mass slightly better although not significant
anova(fit_bm,fit_bm_int) # p > 0.05 -> no interaction is better confirmed
# so let's go ahead with additive mass
visreg(fit = fit_bm, xvar = "TchangeHot",scale = "response", gg = TRUE, line=list(col="black"), rug = FALSE)+
  geom_point(data = SIwide, mapping = aes(x = TchangeHot, y = TchangeCold, shape = Age))+
  #scale_color_gradient(colours = mako(n=100, begin = 0.2, end = 0.8))+
  scale_shape_discrete(labels = c("1-week-old","8-week-old"))+
  xlab("Cloaca temperature change in heat [°C]")+
  ylab("Cloaca temperature change in cold [°C]")+
  #labs(size = "Body mass [kg]")+
  theme_pubr(base_size = 19)+
  theme(legend.title = element_blank())


ggplot(SIwide, aes(x = cold, y = benign))+
  geom_point()+
  theme_pubr()
fit <- lm(cold ~ benign, data = SIwide)
summary(fit) # significant
simulationOutput <- simulateResiduals(fittedModel = fit, n = 500) # simulate data from our model n times
plot(simulationOutput, asFactor=FALSE) # looking ok
visreg(fit = fit, xvar = "benign",scale = "response", gg = TRUE)+
  geom_point(data = SIwide, mapping = aes(x = benign, y = cold))+
  theme_pubr()


#### incorporating morphological measurements ---------------------------------------------------------------------

# overview of chick size (with already summarised data)
ggplot(datTall, aes(y = Mass, x = Age, fill = sex))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(position = position_jitterdodge())+
  scale_fill_manual(values = sexcols)+
  theme_pubr()

dplyr::summarise(group_by(datRNAseq, Age), sd(Mass))


## plot chick size vs. deviation from median temp per treatment
datTall1 = mutate(group_by(filter(datTall, Age=="1week"),treatment), devmed = cloacat-mean(cloacat)) #only 1 week olds
ggplot(subset(datTall1, !is.na(Mass)), aes(x=Mass,y=devmed, color=treatment))+
  geom_point(size=2)+
  scale_color_manual(values= c("#4DAF4A","#377EB8" ,"#E41A1C"))+
  scale_y_continuous(name = "cloaca temp - median cloaca temp per treatment",limits = c(-2.2,2))+
  theme_pubr()
#ggsave("dev_med_cloacat_v_mass_1wk.pdf", path = figurepath)

#same for 8 week olds
datTall8 = mutate(group_by(filter(datTall, Age=="8week"),treatment), devmed = cloacat-mean(cloacat))
ggplot(subset(datTall8, !is.na(Mass)), aes(x=Mass,y=devmed, color=treatment))+
  geom_point(size=2)+
  scale_y_continuous(name = "cloaca temp - median cloaca temp per treatment")+
  scale_color_manual(values= c("#4DAF4A","#377EB8" ,"#E41A1C"))+
  theme_pubr()
#ggsave("dev_med_cloacat_v_mass_8wk.pdf", path = figurepath)
# the cold ones are on top here (esp when plotting mean) which I don't quite get, shouldn't each treatment have same amount above
# and below 0 -> something wrong? I think graph doesn't show me what I think it does

ggplot(filter(datTall, Age%in%c("1week")), aes(x=Mass, y=TchangeAv))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values = unname(tempcols[2:3]))+
  theme_pubr()
# same with mass:neck length ratio
ggplot(filter(datTall, Age%in%c("8week")), aes(x=Mass/NeckLength, y=TchangeAv))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values = unname(tempcols[2:3]))+
  theme_pubr()

ggplot(filter(SIwide, !is.na(TchangeCold) & Age == "1week"), aes(x=Mass, y = TchangeCold))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_pubr()

ggplot(filter(SIwide, Age == "1week"), aes(x=Mass, y = TchangeCold))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_pubr()

## statistical tests
fit1 <- lm(TchangeHot ~ Mass, data = SIwide) # intercept is significant but mass is not
fit2 <- lm(TchangeCold ~ Mass, data = SIwide) # not significant
fit3 <- lm(TchangeAv ~ Mass, data = SIwide) # intercept is significant but mass is not
summary(fit1)

simulationOutput <- simulateResiduals(fittedModel = fit1, n = 500) # simulate data from our model n times
# This method checks if your model is useful by seeing if it can even produce data that looks like the data you used to fit your model (i.e. your observations). 
plot(simulationOutput) # looks ok

options(na.action = "na.fail") # needed for dredge() function to prevent illegal model comparisons
dredgeOut<-dredge(fit1, extra = "R^2") # fit and compare a model set representing all possible predictor combinations
dredgeOut # interaction shouldn't be included but without it doesn't make sense, mass could have an effect

summary(fit2)

simulationOutput <- simulateResiduals(fittedModel = fit2, n = 500) # simulate data from our model n times
# This method checks if your model is useful by seeing if it can even produce data that looks like the data you used to fit your model (i.e. your observations). 
plot(simulationOutput) # looks ok

options(na.action = "na.fail") # needed for dredge() function to prevent illegal model comparisons
(dredgeOut<-dredge(fit2)) # fit and compare a model set representing all possible predictor combinations
dredgeOut # better w/o sex or age

summary(fit3) # even without group, slope is not significant for neither cold or hot, only 1 week or all

aov2 <- aov(TchangeCold ~ NeckLength + Error(group), data = SIwide)
summary(aov2)

fit5 <- lmer(cloacat ~ treatment*Mass + treatment*Age + (1|group) + (1|chickno), data = datRNAseq)
simulationOutput <- simulateResiduals(fittedModel = fit5, n = 500) # simulate data from our model n times
plot(simulationOutput) # ok

options(na.action = "na.fail") # needed for dredge() function to prevent illegal model comparisons
(dredgeOut<-dredge(fit5)) # fit and compare a model set representing all possible predictor combinations
dredgeOut # it's still just treatment


# perhaps find some papers on correlation of body and ambient temp in birds and compare (eg Maloney 2008, or see Gunderson 2024)

## size correlations
ggplot(sampleInfo, aes(x=NeckLength, y= Mass, colour = Age))+
  geom_point(size=3.5, aes(shape=Age))+
  scale_color_manual(values= agecols)+
  labs(x = "Neck length [cm]", y = "Weight [kg]", title = "Weight and neck length of the chicks")+
  theme_pubr(base_size = 24)
ggsave("FrontVBackGirth.pdf", path = figurepath)

p1=ggplot(datTall, aes(x=NeckLength, y= Mass))+
  geom_point(size=2, aes(shape=Age))+
  scale_color_manual(values= tempcols)+
  theme_pubr()
p2=ggplot(SIwide, aes(x=NeckLength, y= FrontGirth))+
  geom_point(size=2, aes(shape=Age))+
  scale_color_manual(values= tempcols)+
  theme_pubr()
p3=ggplot(SIwide, aes(x=BackGirth, y= Mass))+
  geom_point(size=2, aes(shape=Age))+
  scale_color_manual(values= tempcols)+
  theme_pubr()
p4=ggplot(SIwide, aes(x=FrontGirth, y= Mass))+
  geom_point(size=2, aes(shape=Age))+
  scale_color_manual(values= tempcols)+
  theme_pubr()
ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2)
#best between girths and age and back girth, neck seems to be a bit plateauing


### control sex effects
# on morphology
psex1 = ggplot(filter(SIwide, Age == "1week"), aes(x = sex, y = Mass))+
  geom_boxplot()+
  geom_jitter(height = 0)+
  theme_pubr()
psex8 = ggplot(filter(SIwide, Age == "8week"), aes(x = sex, y = Mass))+
  geom_boxplot()+
  geom_jitter(height = 0)+
  theme_pubr()
ggarrange(psex1,psex8, nrow = 1)
# 8 weeks: most females heavier with a few very small ones so more variation, also same with girths but not neck length
# 1 week: basically the same

## on cloaca T
ggplot(datTalive, aes(x = treatment, y = cloacat, fill=sex))+
  geom_boxplot()+
  geom_jitter(position = position_jitterdodge(jitter.height = 0, jitter.width = 0.1))+
  scale_fill_manual(values = sexcols)+
  theme_pubr()


