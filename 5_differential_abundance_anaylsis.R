library(phyloseq)
library(readr)
library(dplyr)
library(phytools)
library(emmeans)
library(tibble)
library(lme4)
library(reshape2)

setwd("/Users/carameyer/Desktop/PhD/exp1repo/ITS")
rm(list=ls())
load("/Users/carameyer/Desktop/PhD/exp1repo/ITS_env.RData")

# select data from the chosen barcode
ITS_ps0 = ps_ITS_most_abund.1

# subset data by fraction
ITS_ps_U = subset_samples(ITS_ps0,ITS_ps0@sam_data$fraction  == "Unfiltered")
ITS_ps_F1 = subset_samples(ITS_ps0,ITS_ps0@sam_data$fraction  == "F1")
ITS_ps_F2 = subset_samples(ITS_ps0,ITS_ps0@sam_data$fraction  == "F2")
ITS_ps_F3 = subset_samples(ITS_ps0,ITS_ps0@sam_data$fraction  == "F3")
ITS_ps_F4 = subset_samples(ITS_ps0,ITS_ps0@sam_data$fraction  == "F4")

# select which fraction to analyze 
##### ALSO CHANGE THE NAME OF THE DATA FRAME AND MODEL AT END OF LOOP
tmp_ps = ITS_ps_F4

# treatments
a = tibble("id"= tmp_ps@sam_data$id,
           "dose"= as.character(tmp_ps@sam_data$dose))
a = factor(a$dose, levels = c("Control", "1x", "10x", "100x"))
# offset
o = log(sample_sums(tmp_ps))
# random effect
z <- as.factor(tmp_ps@sam_data$id)


# r multiple comparaison loop with model

glmITS.sum.global = data.frame()
glmITS.pairwise.global = data.frame()


for (i in 1:ntaxa(tmp_ps)) {
  # select one OTU
  OTU = taxa_names(tmp_ps)[i]
  # response variable
  y = as.vector(tmp_ps@otu_table[OTU,]@.Data)
  
  tryCatch({
    ### model
    glmITS <- glmer(y ~ -1 + a + (1 | z),
                    family='poisson', offset = o)
    
    glmITS.sum = summary(glmITS)$coefficients
    glmITS.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmITS.sum),
                        as_tibble(glmITS.sum))
    glmITS.sum
    glmITS.sum.global = rbind(glmITS.sum.global,glmITS.sum)
    ### multiple comparaison
    glmITS.pairwise = emmeans(glmITS,pairwise~a,adjust="tukey")
    # select p value
    glmITS.pairwise.sum = summary(glmITS.pairwise)
    glmITS.pairwise.sum = glmITS.pairwise.sum[["contrasts"]]
    # extract summary
    tmp_df = glmITS.pairwise.sum
    # keep only comparisons of interest
    tmp = unlist(strsplit(as.character(tmp_df$contrast)," - "))
    tmp_df[,"a"] <- tmp[seq(1,length(tmp),by=2)]
    tmp_df[,"b"] <- tmp[seq(2,length(tmp),by=2)]
    tmp_df = tmp_df[tmp_df$a == "Control" | tmp_df$b == "Control",]
    tmp_df = cbind("OTU"=OTU,tmp_df)
    # extract results in data frame
    glmITS.pairwise.global = rbind(glmITS.pairwise.global,tmp_df)
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  rm(OTU,y,glmITS,glmITS.sum)
}
glmITS.model.global.F4 = glmITS.sum.global
glmITS.pairwise.global.F4 = glmITS.pairwise.global

# adjust p values for all OTUs
glmITS.pairwise.global.U[,"p.adjust"] <- p.adjust(glmITS.pairwise.global.U$p.value,"fdr")
glmITS.pairwise.global.F1[,"p.adjust"] <- p.adjust(glmITS.pairwise.global.F1$p.value,"fdr")
glmITS.pairwise.global.F2[,"p.adjust"] <- p.adjust(glmITS.pairwise.global.F2$p.value,"fdr")
glmITS.pairwise.global.F3[,"p.adjust"] <- p.adjust(glmITS.pairwise.global.F3$p.value,"fdr")
glmITS.pairwise.global.F4[,"p.adjust"] <- p.adjust(glmITS.pairwise.global.F4$p.value,"fdr")

save.image("/Users/carameyer/Desktop/PhD/exp1repo/ITS_env.RData")