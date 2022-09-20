library(export)
library(UpSetR)
library(dplyr)
library(tidyr)
library(grid)
library(ggplot2)

setwd("/Users/carameyer/Desktop/PhD/exp1repo/ITS")
#rm(list=ls())
#load("/Users/carameyer/Desktop/PhD/exp1repo/ITS_env.RData")

# Find significant effects of pesticide in fractions ----
# subset significant changes from differential abundance analysis; add column identifying fraction
ITS.sig.otu.U <- subset(glmITS.pairwise.global.U, p.adjust < 0.05) %>%
  mutate(fraction = "Unfiltered") %>%
  mutate(comparison = paste(fraction, contrast))
ITS.sig.otu.F1 <- subset(glmITS.pairwise.global.F1, p.adjust < 0.05) %>%
  mutate(fraction = "F1") %>%
  mutate(comparison = paste(fraction, contrast))
ITS.sig.otu.F2 <- subset(glmITS.pairwise.global.F2, p.adjust < 0.05) %>%
  mutate(fraction = "F2") %>%
  mutate(comparison = paste(fraction, contrast))
ITS.sig.otu.F3 <- subset(glmITS.pairwise.global.F3, p.adjust < 0.05) %>%
  mutate(fraction = "F3") %>%
  mutate(comparison = paste(fraction, contrast))
ITS.sig.otu.F4 <- subset(glmITS.pairwise.global.F4, p.adjust < 0.05) %>%
  mutate(fraction = "F4") %>%
  mutate(comparison = paste(fraction, contrast))

ITS.sig.otu.all <- rbind(ITS.sig.otu.U, ITS.sig.otu.F1, ITS.sig.otu.F2, ITS.sig.otu.F3, ITS.sig.otu.F4)
ITS.sig_otus <- unique(ITS.sig.otu.all$OTU) #49

# Since it isn't very many OTUs, we'll look at the OTU table
ITS_sig_otu_table <- as.data.frame(otu_table(ps_ITS_most_abund.1)) %>%
  subset(rownames(.) %in% ITS.sig_otus)

# count the number of OTUs per comparison 
ITS.OTU.comp.all <- ITS.sig.otu.all %>%
  group_by(comparison) %>%
  summarise(OTU = list(unique(OTU)))
head(ITS.OTU.comp.all)
# Significant estimates per OTU in each comparison 
ITS.sig.all <- ITS.sig.otu.all %>%
  group_by(comparison, OTU) %>%
  summarise(estimates = estimate)
head(ITS.sig.all) 

ITS.sig.all.df <- tibble(ITS.sig.all)

#transform from long to wide tibble
ITS.sig.all.df<- ITS.sig.all.df %>%
  spread(comparison, estimates)
ITS.sig.all.df[is.na(ITS.sig.all.df)] <- 0
head(ITS.sig.all.df)
dim(ITS.sig.all.df)

#write.table(x=ITS.sig.E, file="ITS.sig.Estimates.csv")

##how many OTUs significant per comparison
nb.ITS.sig.all <- colSums(ITS.sig.all.df != 0) %>% as.data.frame()
nb.ITS.sig.all
write.table(x=nb.ITS.sig.all, file="out/nb.ITS.sigificant.csv", sep = ";")

###Add % -> % significant OTUs of most abundant OTUs per comparison
nb.ITS.sig.all.perc <- (colSums(ITS.sig.all.df != 0)/ nrow(tax_table(ps_ITS_most_abund.1)))*100 
nb.ITS.sig.all.perc <- as.data.frame(nb.ITS.sig.all.perc)
nb.ITS.sig.all.perc
write.table(x=nb.ITS.sig.all.perc, file="out/nb.ITS.sigificant.all.perc.csv", sep = ";")

# Make into a data frame for saving UpSet results later
ITS.affected.otus <- data.frame("Affected OTUs" = nb.ITS.sig.all$., "Percentage of most abundant OTUs" = nb.ITS.sig.all.perc$nb.ITS.sig.all.perc,
                            row.names = row.names(nb.ITS.sig.all.perc))
rownames(ITS.affected.otus)[1] <- "Any Treatment"

write.table(x=ITS.affected.otus, file="out/ITS.affected.otus.csv", sep = ";", row.names = TRUE, col.names = NA)

# what if we make the condition that they have to be present in 60% of the replicates of any treatment in that fraction?
prev.any_F1 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F1 >= 60 | tmp_otu_prev_ttt$`1x_F1` >= 60 |
                        tmp_otu_prev_ttt$`10x_F1` >= 60 | tmp_otu_prev_ttt$`100x_F1` >= 60)
length(unique(rownames(prev.any_F1))) #517
otu.prev.any.F1 <- unique(rownames(prev.any_F1))
length(intersect(unique(sig.otu.F1$OTU), unique(rownames(prev.any_F1)))) #329
length(unique(sig.otu.F1$OTU)) #329 so is it all of them now?
not_overlapping_F1_any <- setdiff(unique(rownames(prev.any_F1)), unique(sig.otu.F1$OTU))

prev.any_F2 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F2 >= 60 | tmp_otu_prev_ttt$`1x_F2` >= 60 |
                        tmp_otu_prev_ttt$`10x_F2` >= 60 | tmp_otu_prev_ttt$`100x_F2` >= 60)
length(unique(rownames(prev.any_F2))) #359

prev.any_F3 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F3 >= 60 | tmp_otu_prev_ttt$`1x_F3` >= 60 |
                        tmp_otu_prev_ttt$`10x_F3` >= 60 | tmp_otu_prev_ttt$`100x_F3` >= 60)
length(unique(rownames(prev.any_F3))) #494

prev.any_F4 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F4 >= 60 | tmp_otu_prev_ttt$`1x_F4` >= 60 |
                        tmp_otu_prev_ttt$`10x_F4` >= 60 | tmp_otu_prev_ttt$`100x_F4` >= 60)
length(unique(rownames(prev.any_F4))) #367

prev.any_U <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_Unfiltered >= 60 | tmp_otu_prev_ttt$`1x_Unfiltered` >= 60 |
                       tmp_otu_prev_ttt$`10x_Unfiltered` >= 60 | tmp_otu_prev_ttt$`100x_Unfiltered` >= 60)
length(unique(rownames(prev.any_U))) #530


all.pres <- left_join("OTU" = all.df$OTU, "present U" = ifelse(otu.prev.any.F1 %in% all.df$OTU, 1, 0))
all.pres <- data.frame(matrix(ncol = 5, nrow = 473))
colnames(all.pres) <- c("U", "F1", "F2", "F3", "F4")
rownames(all.pres) <- all.df$OTU
all.pres[,1] <- apply(all.pres[, 1], 1, function(x) ifelse(rownames(x) %in% rownames(prev.any_U), 1, 0))
all.pres[,1] 
all.pres$U <- ifelse(rownames(all.pres) %in% rownames(prev.any_U), 1, 0)
all.pres$F1 <- ifelse(rownames(all.pres) %in% rownames(prev.any_F1), 1, 0)

prev.control_F2 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F2 >= 60) #from load_data script
length(unique(rownames(prev.control_F2))) #270
otu.pres.F2 <- intersect(unique(rownames(prev.control_F2)), unique(sum.present_F2$OTU))
length(otu.pres.F2) #244

prev.control_F3 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F3 >= 60) #from load_data script
length(unique(rownames(prev.control_F3))) #412
otu.pres.F3 <- intersect(unique(rownames(prev.control_F3)), unique(sum.present_F3$OTU))
length(otu.pres.F3) #378

prev.control_F4 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F4 >= 60) #from load_data script
length(unique(rownames(prev.control_F4))) #319
otu.pres.F4 <- intersect(unique(rownames(prev.control_F4)), unique(sum.present_F4$OTU))
length(otu.pres.F4) #294


save.image("/Users/carameyer/Desktop/PhD/exp1repo/ITS_env.RData")
