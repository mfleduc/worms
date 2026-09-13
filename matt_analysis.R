library("DescTools")
library("ggplot2")
library("geomtextpath")
library("glmmTMB")
library("DHARMa")
library("ggdist")
library("usethis")
library("effects")
library("car")
library("emmeans")
library("png")
library("ggpubr")
library("tibble")
library("grid")
library("tidytable")
library("parameters")
#install.packages("parameters")

### Table of Contents
# 1 - File Management & Data
# 2 - Visualizations
# 3 - Summary Statistics
# 4 - Model Comparisons
# 5 - Best Fit Model & Tukey's Test

### Analysis
# 1. set up data
# 2. visualizations
# 3. determine best model
# 4. Use model
# 5. Conduct Tukey's

##################################################
# 1 - FILE MANAGEMENT & DATA
##################################################

# Read Files for Gardner
wfish <- read.csv(header = TRUE, "dec16revis.csv")
wfish

# Renaming/Ordering/Formatting Fish Names
wfish$SPECIES[wfish$SPECIES=="Red"] <- "Sockeye"
wfish$SPECIES <- factor(wfish$SPECIES, levels=c("Pink", "Sockeye", "Chum"))
wfish$year <- c(rep(as.numeric(2024)))

wfishes <- lm(wfish$W.G.INDEX ~ unlist(wfish$SPECIES))

# Read + Format Files for Mastick
mast <- read.csv("Can Data_Jan23.csv")
mast <- mast[mast$box.no.!="Practice",]
mast <- mast[mast$can.size!="small" & mast$can.size!="med",]
mast$salmon.species[mast$salmon.species=="coho"] <- "Coho"
mast$salmon.species[mast$salmon.species=="pink"] <- "Pink"
mast$salmon.species[mast$salmon.species=="red"] <- "Sockeye"
mast$salmon.species[mast$salmon.species=="chum"] <- "Chum"
mast$salmon.species <- factor(mast$salmon.species, levels=c("Coho", "Pink", "Sockeye", "Chum"))
mast$wgindex <- as.numeric(mast$wgindex)
# Testing and removing outlier chum
zc <- (0.2750127 - mean(mast$wgindex))/sd(mast$wgindex)
zc
mast <- mast[mast$wgindex<="0.2750127",]

mast$whg <- as.numeric(mast$whg)
mast$year <- as.numeric(mast$year)

wfish$set <- c("gardner")
mast$set <- c("mastick")

combo <- data.frame(
  species <- c(mast$salmon.species, wfish$SPECIES),
  whg <- c(mast$whg, wfish$WHG),
  year <- c(mast$year, wfish$year),
  set <- c(mast$set, wfish$set)
  
)
names(combo)[names(combo) == 'species....c.mast.salmon.species..wfish.SPECIES.'] <- 'species'
names(combo)[names(combo) == 'year....c.mast.year..wfish.year.'] <- 'year'
names(combo)[names(combo) == 'whg....c.mast.whg..wfish.WHG.'] <- 'whg'
names(combo)[names(combo) == 'set....c.mast.set..wfish.set.'] <- 'set'
combo$whg <- as.numeric(combo$whg)

comboModern = combo[combo$year >= "2000",]

comboHist = combo[combo$year < "2000",]

##################################################
# 2 - VISUALIZATIONS
##################################################

# Regional Distribution Intensity Boxplots
boxplot(wfish$WHG ~ wfish$REGION,
        main="Regional Comparison of 2024 Samples",
        xlab="", ylab="Nematodes per 100g",
        col=c("slateblue" , "royalblue", "purple"),
        border=c("slateblue4", "royalblue4", "purple4"))

ggplot(wfish, aes(x=REGION, y=WHG, fill=REGION)) +
  #background_image(akm1) +
  #geom_point(color = "red", size = 5) +
  geom_boxplot() +
  scale_fill_manual(values=c("slateblue" , "royalblue", "purple") )
#annotation_custom(rasterGrob(im2,  width = unit(1,"npc"),  height = unit(1,"npc")), -Inf, Inf, -Inf, Inf))

# Half-eye for Regional Distribution
ggplot (wfish, aes(y = REGION, x = WHG)) +
  stat_halfeye()

# Boxplot of Species for 2024
boxplot(wfish$WHG ~ wfish$SPECIES,
        main="Species Comparison of 2024 Samples",
        xlab="Species", ylab="Nematodes per 100g",
        #xlab="Species",
        #ylab="Infection Index (nematodes/gram)",
        col=c("pink" , "tomato", "#a7cdd6"),
        border=c("#f06790", "red", "#5894a3")
)

# Boxplot of Mastick data 
boxplot(mast$whg ~ mast$salmon.species,
        main="Species Comparison of Samples from 1982-2020\n(Mastick et al. 2024)",
        xlab="Species", ylab="Nematodes per 100g",
        col=c("gray", "pink", "tomato", "#a7cdd6"),
        border=c("darkgray", "#f06790", "red", "#5894a3"))


# Mastick WHG over time
ggplot(mast, aes(x = year, y = whg, color = salmon.species)) +
  geom_point() +
  geom_labelsmooth(aes(label = salmon.species), fill = "white",
                   method = "lm", formula = y ~ x,
                   size = 3, linewidth = 1, boxlinewidth = 0.4) +
  xlab ("Year") +
  xlim (1979, 2020) +
  theme_bw() + guides(color = 'none')

# VISUALIZATION - Combined sets of WHG over time
# wfish$year <- c(as.numeric("2024"))

ggplot(combo, aes(x = year, y = whg, color = species)) +
  geom_point() +
  geom_labelsmooth(aes(label = species), fill = "white",
                   method = "lm", formula = y ~ x,
                   size = 3, linewidth = 1, boxlinewidth = 0.4) +
  xlab ("Year") +
  xlim (1979, 2024) +
  theme_bw() + guides(color = 'none')

# VISUALIZATION - Distributions of WHG per data set
# Gardner 2024
ggplot (wfish, aes(y = SPECIES, x = WHG)) +
  stat_halfeye()

# Mastick 1979-2020
ggplot (mast, aes(y = salmon.species, x = whg)) +
  stat_halfeye()

# Combined 1979-2024
ggplot (combo, aes(y = species, x = whg)) +
  stat_halfeye()
# This tells us the data is not normally distributed
# and we need to use a zero-inflated gamma model.

# VISUALIZATION - Data sets on box plot together
ggplot(combo, aes(x=species, 
                  y=whg, 
                  fill=set,
                  #colour = 'red',
)) + 
  scale_fill_manual(values = c("tomato", "#a7cdd6")) +
  geom_boxplot() +
  ylab ("Nematodes per 100g") +
  xlab ("Species") 

### FIGURE X - Combined data set box plot

boxplot(combo$whg ~ combo$species,
        #subset=(threemast$salmon.species!="Coho"),
        main="Species Comparison of Combined Samples",
        xlab="Species", ylab="Nematodes per 100g",
        col=c("gray", "pink", "tomato", "#a7cdd6"),
        border=c("darkgray", "#f06790", "red", "#5894a3"))

##################################################
# 3 - SUMMARY STATS
##################################################

# SUMMARY - Mastick
as.numeric(mast$whg)
MPinks<-mast[mast$salmon.species=="Pink",]
MReds<-mast[mast$salmon.species=="Sockeye",]
MChum<-mast[mast$salmon.species=="Chum",]
MCoho<-mast[mast$salmon.species=="Coho",]

MPinks$whg<-as.numeric(MPinks$whg)
MReds$whg<-as.numeric(MReds$whg)
MChum$whg<-as.numeric(MChum$whg)
MCoho$whg<-as.numeric(MCoho$whg)

MStatsSum <- data.frame(
  Species = c("Pink", "Sockeye", "Chum", "Coho"),
  SampleSize = c(nrow(MPinks), nrow(MReds), nrow(MChum), nrow(MCoho)),
  Mean = c(mean(MPinks$whg), mean(MReds$whg), mean(MChum$whg), mean(MCoho$whg)),
  StDev = c(sd(MPinks$whg), sd(MReds$whg), sd(MChum$whg), sd(MCoho$whg))
)
MStatsSum

# SUMMARY - Gardner
as.numeric(wfish$WHG)
GPinks<-wfish[wfish$SPECIES=="Pink",]
GReds<-wfish[wfish$SPECIES=="Sockeye",]
GChums<-wfish[wfish$SPECIES=="Chum",]

GStatsSum <- data.frame(
  Species = c("Pink", "Sockeye", "Chum"),
  SampleSize = c(nrow(GPinks), nrow(GReds), nrow(GChums)),
  Mean = c(mean(GPinks$WHG), mean(GReds$WHG), mean(GChums$WHG)),
  StDev = c(sd(GPinks$WHG), sd(GReds$WHG), sd(GChums$WHG))
)
GStatsSum


##################################################
# 4 - MODEL COMPARISONS
##################################################
mask = comboModern[is.finite(comboModern$whg),]
comboModern=comboModern[is.finite(comboModern$whg),]
mask_coho = comboModern$species!="Coho"
comboModern = comboModern[mask_coho,]

# Conway-Maxwell Poisson, compois(link = "log")
# modConway = glmmTMB(  
#   whg ~ species + (1|year), 
#   data = combo,
#   family=
#     compois(link = "log"))
# print(summary(modConway),show.residuals=TRUE)

# tweedie
modTwee = glmmTMB(  
  whg ~ species + (1|year), 
  data = comboModern,
  family= tweedie(link = "log"))

# ziGamma log
modZiGam = glmmTMB(  
  whg ~ species + (1|year), 
  data = comboModern,
  ziformula = ~.,
  family=ziGamma("log"))

# ziGamma inverse
modZiGamIn = glmmTMB(  
  whg ~ species + (1|year), 
  data = comboModern,
  ziformula = ~.,
  family=ziGamma("inverse"))

# best fit
models <- list(modTwee, modZiGam,modZiGamIn)
model_names <- c("Tweedie", "ziGamma Log", "ZI-Gaminv")

residual_summaries <- map2_df(models, model_names, function(model, name) {
  sim <- simulateResiduals(fittedModel = model, plot = FALSE, n=2500)
  data.frame(
    model = name,
    dispersion = testDispersion(sim)$p.value,
    uniformity = testUniformity(sim)$p.value,
    outliers = testOutliers(sim)$p.value,
    AIC = AIC(model)
  )
})

residual_summaries %>%
  arrange(AIC) %>%
  print()

modFinal <- modTwee

simresults = simulateResiduals(fittedModel = modFinal, plot = TRUE, n=2500)
summary(modFinal)

##################################################
# 5 - MODEL CHOICE & TUKEY'S TEST
##################################################

## Regional comparison using 2024: Likelihood ratio test
## Find that region is not a significant contributor

gModTwee_regional = glmmTMB(  
  WHG ~ REGION, 
  data = wfish,
  family= tweedie(link = "log"))
gModTwee_null <- glmmTMB(
  WHG ~ 1,
  data = wfish,
  family = tweedie(link = "log")
)
anova(gModTwee_null, gModTwee_regional)

print(summary(gModTwee_regional),show.residuals=TRUE)

## Species comparison using 2000-2020 + 2024:

cModAnovaSpecies = glmmTMB(
  whg ~ -1+species + (1|year),
  data = comboModern,
  family=tweedie("log"),
)
print(summary(cModAnovaSpecies),show.residuals=TRUE)

# anova_table = glmmTMB:::Anova.glmmTMB(cModAnovaSpecies,type=3)
# anova_table

# Estimate ratio of means in physical space rather than log space, 
# then look at their ratios. CI doesn't include 1 -> Reject null
# Internally, everything is done on the log scale and then transformed 
# into physical space for our analysis
emm_cond = emmeans(cModAnovaSpecies, ~species, type="response")
tukey_cond = pairs(emm_cond, adjust = "tukey")
tukey_cond
confint(tukey_cond)



## Species comparison using 2024: 

gModTwee = glmmTMB(  
  WHG ~ -1+SPECIES, 
  data = wfish,
  family= tweedie(link = "log"))

print(summary(gModTwee),show.residuals=TRUE)

######### Probability of infection
data_p0_testing = wfish[9:dim(wfish)[1],]
infection_yes_no = (data_p0_testing$WHG>1e-6)
data_p0_testing$infected = infection_yes_no
contingency_table = array(0,dim=c(2,2))
species_names = unique(data_p0_testing$SPECIES)
for(sn in 1:length(species_names)){
  contingency_table[sn,1] = sum( infection_yes_no&(data_p0_testing$SPECIES==species_names[sn]) )
  contingency_table[sn,2] = sum( (!infection_yes_no)&(data_p0_testing$SPECIES==species_names[sn]) )
}
chisq_res = chisq.test(contingency_table)#[2:3,])
fisher_res = fisher.test(contingency_table)#[2:3,])
barnard_res = BarnardTest(contingency_table)#[2:3,])
boschloo_res = BarnardTest(contingency_table)#[2:3,], method="boschloo")
## Could try logistic regression with can size as a predictor
## Something like
## log(p) = \beta_species+\beta_year+\beta_{size}log(cansize)
## I worry that these datasets are too small to get good results. All of the ones above are 
## hovering on the edge of significance, which, while suggestive, mostly 
## points toward sample size limitations that we can solve by including the "whole" dataset
m_offset <- glmmTMB(
  infected ~ -1+SPECIES + offset(log(WEIGHT)),
  family = binomial(link = "cloglog"),
  data =data_p0_testing
)
## Note: slight variations in mass examined, so let's account for it explicitly.

emmeans_binom = emmeans(m_offset, ~SPECIES, type="response")
tukey_binom = pairs( emmeans_binom, adjust="tukey")
emm_rg <- regrid(emmeans_binom)
tukey_rg = pairs(emm_rg)
confint(tukey_rg)
## NOTE: Use the p-values from tukey_binom and confidence intervals from tukey_rg. 
## p-values calculated where the effects are linear are more accurate
## but regridding to get them out of the cloglog space is more interpretable. 

# m_free <- glmmTMB(
#   infected ~ -1+SPECIES + log(WEIGHT),
#   family = binomial(link = "cloglog"),
#   data = data_p0_testing
# )



