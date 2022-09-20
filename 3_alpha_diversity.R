library(ggplot2)
library(ggpubr)
library(MASS)
library(multcomp)
library(rstatix)
library(emmeans)
library(export)

install.packages("ggpubr")

tmp_alphadiv = alpha_div_ITS_exp1

# Make Fraction and Dose factors with ordered levels
tmp_alphadiv$dose <- factor(tmp_alphadiv$dose, levels = c("Control", "1x", "10x", "100x"))
tmp_alphadiv$fraction <- factor(tmp_alphadiv$fraction, levels = c("Unfiltered", "F1", "F2", "F3", "F4"))

# set which diversity index you want to analyze
# "chao1", "dominance", "equitability", "observed_species", "shannon", "simpson_reciprocal"
tmp_index = "simpson_reciprocal"

# subset to include only your index of interest while keeping mapping data
tmp_sub <- tmp_alphadiv[,c(colnames(tmp_alphadiv[,1:6]),tmp_index)]

# set index name to "value" for ease of executing code
colnames(tmp_sub)[7] <- "value"


# subset data by fraction ----
tmp_U <- subset(tmp_sub, fraction == "Unfiltered")
tmp_F1 <- subset(tmp_sub, fraction == "F1")
tmp_F2 <- subset(tmp_sub, fraction == "F2")
tmp_F3 <- subset(tmp_sub, fraction == "F3")
tmp_F4 <- subset(tmp_sub, fraction == "F4")

# set which subset you want to analyze
tmp_data = tmp_F4

# ANOVA ----
#tmp_aov <- aov(value ~ fraction * dose, data = tmp_sub)

tmp_aov <- aov(value ~ dose, data = tmp_data)
summary(tmp_aov)
#plot(tmp_aov)

# check normality of residuals; with shapiro test p value over 0.05 means we can assume normality
ggqqplot(residuals(tmp_aov))
shapiro_test(residuals(tmp_aov))
#dominance F3 p = 0.005
#observed_species F3 p = 0.035
#PD_whole_tree F3 p = 0.0386
#simpson_reciprocal F1 p = 0.0449
#simpson_reciprocal F2 p = 0.00219

# check for homogeneity of variance
plot(tmp_aov, 1) 
tmp_aov %>% levene_test(value ~ dose) #if p value greater than 0.05, we can assume homogeneity of variance
#dominance F1 p = 0.01
#observed_species F3 p = 0.0201
#equitability F3 p = 0.026
#shannon F2 p = 0.0229
#shannon F3 p = 0.0113
#simpson_reciprocal U p = 0.0396
#simpson_reciprocal F3 p = 0.0220

# post hoc Tukey test
tmp_ht <- glht(tmp_aov, linfct = mcp(dose = "Tukey"))
summary(tmp_ht) # get p values
tmp_cld <- cld(tmp_ht, decreasing = FALSE)
tmp_cld$mcletters$Letters # get compact letters

# check for outliers
tmp_data %>% 
  group_by(dose) %>%
  identify_outliers(value)

# Plot ----
tmp_plot <-ggplot(tmp_sub, aes(x=dose, y=value)) + 
  geom_boxplot() +
  #geom_text(aes(x = Pesticide, y = tops, label = cld), inherit.aes = FALSE, data = anno.cld) +
  ggtitle(tmp_index) +
  xlab("Pesticide Dose") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(tmp_sub$fraction~., ncol = 5)

tmp_plot

graph2ppt(x= tmp_plot ,file="out/simpson_reciprocal.pptx")

# GLM ----
# glm with negative binomial for observed species, gamma for chao1, shannon, simpson reciprocal, PD whole tree
tmp_glm <- glm(value ~ dose, data = tmp_data, family = "Gamma") 
summary(tmp_glm)
plot(tmp_glm)
ggqqplot(residuals(tmp_glm))

hist(tmp_sub$value)
hist(residuals(tmp_glm))

# tmp_glm <- glm.nb(value ~ dose, data = tmp_data) # for observed species
# summary(tmp_glm)

tmp_ht <- glht(tmp_glm, linfct = mcp(dose = "Tukey"))
summary(tmp_ht) # get p values
tmp_cld <- cld(tmp_ht, decreasing = FALSE)
tmp_cld$mcletters$Letters # get compact letters

library(FSA)
library(rcompanion)

kruskal.test(value ~ dose, data = tmp_data) #chi-squared = 18.195, df = 3, p-value = 0.000401
d <- dunnTest(tmp_data$value, tmp_data$dose, method = "bh")
cld <- cldList(P.adj ~ Comparison, data = d$res, threshold = 0.05, remove.zero = FALSE)
cld

# remove temporary objects from working environment
rm(list = names(.GlobalEnv)[grep("tmp",names(.GlobalEnv))])
