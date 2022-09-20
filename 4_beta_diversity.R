library(phyloseq)
library(ggplot2)
library(tibble)
library(dplyr)
library(vegan)
library(export)
library(devtools)
library(pairwiseAdonis)
library(usedist)

# rarefy phyloseq object to 11000 sequences
rare_ps = rarefy_even_depth(ps_ITS,sample.size = 2000, rngseed = 11)

# remove singletons and doubletons
rare_3_ps <- prune_taxa(taxa_sums(rare_ps) > 2, rare_ps)

# Bray Curtis ----
# select samples coordinates for ordination plot
bc_data_2D <- tibble("id" = as.numeric(rownames(bc_ordin[["vectors"]])),
                     "x" = as.numeric(bc_ordin[["vectors"]][,1]),
                     "y" = as.numeric(bc_ordin[["vectors"]][,2]))
# add treatment, color, order, shape, etc
bc_data_2D <- left_join(bc_data_2D,rare_3_ps@sam_data[,c("id","treatment","dose","fraction")], by = "id")

## calculate distance matrix and ordination ----
bc_dist = distance(rare_3_ps,"bray")
# calculate ordination
bc_ordin = ordinate(rare_3_ps, "PCoA", distance = bc_dist)

bc_disp <- with(bc_data_2D, betadisper(bc_dist, treatment))
bc_disp
anova(bc_disp)
TukeyHSD(bc_disp)
plot(bc_disp)
boxplot(bc_disp)


## plot PCoA ----

# data wrangling for plot


# select axis titles
bc_x_lab = round(bc_ordin[["values"]][["Relative_eig"]][1]*100,2)
bc_x_lab <- paste0("Axis 1: ",bc_x_lab," %")
bc_y_lab = round(bc_ordin[["values"]][["Relative_eig"]][2]*100,2)
bc_y_lab <- paste0("Axis 2: ",bc_y_lab," %")


bc_p <- ggplot(bc_data_2D, aes(x=x,y=y)) +
  geom_point(aes(shape = fraction, fill = dose, stroke = 1, size = 4)) +
  scale_fill_manual(breaks = c("Control", "1x", "10x", "100x"),
                    values=c("#009E73", "#F0E442", "#E69F00", "#D55E00")) +
  discrete_scale("shape", "shape", palette=function(n) { # to get shapes that work with "fill"
    stopifnot("more than 5 shapes not supported" = n <= 5)
    20 + seq_len(n)
  }) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5))) +
  theme_bw() +
  ggtitle("PCoA Bray-Curtis 16S experiment 1") + 
  ylab(bc_y_lab) +
  xlab(bc_x_lab) 
#stat_ellipse(aes(x=x, y=y, group=treatment),type = "norm")

bc_p

graph2ppt(bc_p, "out/PCoA_bray.pptx")

## PERMANOVA ----
adonis2(bc_dist ~ dose, data=bc_data_2D)

adonis2(bc_dist ~ fraction, data=bc_data_2D)

# bivariate (terms added sequentially)
perm_bc <- adonis2(bc_dist ~ dose * fraction, data=bc_data_2D)
perm_bc

perm_bc_df <- as.data.frame(perm_bc)
write.table(perm_bc_df, file="out/bray.curtis.permanova.csv", sep = ";", col.names = NA)

bc_adonis_pair <- pairwise.adonis(as.matrix(bc_dist), factors=as.factor(rare_3_ps@sam_data$treatment),p.adjust.m='BH')

bc_adonis_pair_unf = pairwise.adonis(as.matrix(bc_dist), factors=as.factor(rare_3_ps@sam_data$treatment),p.adjust.m='BH', reduce = 'Unfiltered')

## Plot distances and pairwise PERMANOVA for each fraction ----
fractions <- rare_3_ps@sam_data$fraction
U <- subset_samples(rare_3_ps,rare_3_ps@sam_data$fraction  == "Unfiltered")
dose_U <- U@sam_data$dose 
F1 <- subset_samples(rare_3_ps,rare_3_ps@sam_data$fraction  == "F1")
dose_F1 <- F1@sam_data$dose 
F2 <- subset_samples(rare_3_ps,rare_3_ps@sam_data$fraction  == "F2")
dose_F2 <- F2@sam_data$dose 
F3 <- subset_samples(rare_3_ps,rare_3_ps@sam_data$fraction  == "F3")
dose_F3 <- F3@sam_data$dose 
F4 <- subset_samples(rare_3_ps,rare_3_ps@sam_data$fraction  == "F4")
dose_F4 <- F4@sam_data$dose 

### Unfiltered ----
bc_dist_U <- dist_subset(bc_dist, fractions == "Unfiltered")
bc_calc_U <- dist_groups(as.dist(bc_dist_U), dose_U)
bc_calc_U_tab <- as_tibble(bc_calc_U)
bc_adonis_pair_U = pairwise.adonis(as.matrix(bc_dist_U), factors=as.factor(dose_U),p.adjust.m='bonferroni')
write.table(bc_adonis_pair_U, file="out/bray.curtis.unfiltered.csv", sep = ";", col.names = NA)

bc_U_box <- ggplot(bc_calc_U_tab, aes(x = Label, y = Distance)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.1,
              size = 1) + theme_bw()  +
  geom_boxplot(colour = "black", alpha=0.1, outlier.colour = NA) +
  xlab("")+ylab("")+
  ggtitle("Bray-Curtis Dissimilarities for Unfiltered")+
  theme(plot.title = element_text(size = rel(1.2)),
        legend.position="none")+
  theme(axis.text.x =  element_text(angle=45, hjust = 1),
        strip.text.x= element_text(size = 13),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size=12)) +
  scale_x_discrete(limits = rev)
bc_U_box

graph2ppt(bc_U_box, "out/bray.curtis.unfiltered.pptx")

### F1 ----
bc_dist_F1 <- dist_subset(bc_dist, fractions == "F1")
bc_calc_F1 <- dist_groups(as.dist(bc_dist_F1), dose_F1)
bc_calc_F1_tab <- as_tibble(bc_calc_F1)
bc_adonis_pair_F1 = pairwise.adonis(as.matrix(bc_dist_F1), factors=as.factor(dose_F1),p.adjust.m='bonferroni')
write.table(bc_adonis_pair_F1, file="out/bray.curtis.F1.csv", sep = ";", col.names = NA)

bc_F1_box <- ggplot(bc_calc_F1_tab, aes(x = Label, y = Distance)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.1,
              size = 1) + theme_bw()  +
  geom_boxplot(colour = "black", alpha=0.1, outlier.colour = NA) +
  xlab("")+ylab("")+
  ggtitle("Bray-Curtis Dissimilarities for F1")+
  theme(plot.title = element_text(size = rel(1.2)),
        legend.position="none")+
  theme(axis.text.x =  element_text(angle=45, hjust = 1),
        strip.text.x= element_text(size = 13),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size=12)) +
  scale_x_discrete(limits = rev)
bc_F1_box

graph2ppt(bc_F1_box, "out/bray.curtis.F1.pptx")

### F2 ----
bc_dist_F2 <- dist_subset(bc_dist, fractions == "F2")
bc_calc_F2 <- dist_groups(as.dist(bc_dist_F2), dose_F2)
bc_calc_F2_tab <- as_tibble(bc_calc_F2)
bc_adonis_pair_F2 = pairwise.adonis(as.matrix(bc_dist_F2), factors=as.factor(dose_F2),p.adjust.m='bonferroni')
write.table(bc_adonis_pair_F2, file="out/bray.curtis.F2.csv", sep = ";", col.names = NA)

# F2 is unbalanced because one replicate of treatment C_F2 is missing (low sequence count)
# because PERMANOVA can be sensitive to heterogeneity of variance in an unbalanced design, use betadisper to check homogeneity of variance in this fraction
bc_data_F2 <- subset(bc_data_2D, fraction =="F2")
bc_disp_F2 <- with(bc_data_F2, betadisper(bc_dist_F2, dose))
bc_disp_F2
anova(bc_disp_F2)
TukeyHSD(bc_disp_F2)
plot(bc_disp_F2)
boxplot(bc_disp_F2)

bc_F2_box <- ggplot(bc_calc_F2_tab, aes(x = Label, y = Distance)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.1,
              size = 1) + theme_bw()  +
  geom_boxplot(colour = "black", alpha=0.1, outlier.colour = NA) +
  xlab("")+ylab("")+
  ggtitle("Bray-Curtis Dissimilarities for F2")+
  theme(plot.title = element_text(size = rel(1.2)),
        legend.position="none")+
  theme(axis.text.x =  element_text(angle=45, hjust = 1),
        strip.text.x= element_text(size = 13),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size=12)) +
  scale_x_discrete(limits = rev)
bc_F2_box

graph2ppt(bc_F2_box, "out/bray.curtis.F2.pptx")

### F3 ----
bc_dist_F3 <- dist_subset(bc_dist, fractions == "F3")
bc_calc_F3 <- dist_groups(as.dist(bc_dist_F3), dose_F3)
bc_calc_F3_tab <- as_tibble(bc_calc_F3)
bc_adonis_pair_F3 = pairwise.adonis(as.matrix(bc_dist_F3), factors=as.factor(dose_F3),p.adjust.m='bonferroni')
write.table(bc_adonis_pair_F3, file="out/bray.curtis.F3.csv", sep = ";", col.names = NA)

bc_F3_box <- ggplot(bc_calc_F3_tab, aes(x = Label, y = Distance)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.1,
              size = 1) + theme_bw()  +
  geom_boxplot(colour = "black", alpha=0.1, outlier.colour = NA) +
  xlab("")+ylab("")+
  ggtitle("Bray-Curtis Dissimilarities for F3")+
  theme(plot.title = element_text(size = rel(1.2)),
        legend.position="none")+
  theme(axis.text.x =  element_text(angle=45, hjust = 1),
        strip.text.x= element_text(size = 13),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size=12)) +
  scale_x_discrete(limits = rev)
bc_F3_box

graph2ppt(bc_F3_box, "out/bray.curtis.F3.pptx")

### F4 ----
bc_dist_F4 <- dist_subset(bc_dist, fractions == "F4")
bc_calc_F4 <- dist_groups(as.dist(bc_dist_F4), dose_F4)
bc_calc_F4_tab <- as_tibble(bc_calc_F4)
bc_adonis_pair_F4 = pairwise.adonis(as.matrix(bc_dist_F4), factors=as.factor(dose_F4),p.adjust.m='bonferroni')
write.table(bc_adonis_pair_F4, file="out/bray.curtis.F4.csv", sep = ";", col.names = NA)

# F4 is unbalanced because one replicate of treatment xxF4 is missing (low sequence count)
# because PERMANOVA can be sensitive to heterogeneity of variance in an unbalanced design, use betadisper to check homogeneity of variance in this fraction
bc_data_F4 <- subset(bc_data_2D, fraction =="F4")
bc_disp_F4 <- with(bc_data_F4, betadisper(bc_dist_F4, dose))
bc_disp_F4
anova(bc_disp_F4)
TukeyHSD(bc_disp_F4)
plot(bc_disp_F4)
boxplot(bc_disp_F4)

bc_F4_box <- ggplot(bc_calc_F4_tab, aes(x = Label, y = Distance)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.1,
              size = 1) + theme_bw()  +
  geom_boxplot(colour = "black", alpha=0.1, outlier.colour = NA) +
  xlab("")+ylab("")+
  ggtitle("Bray-Curtis Dissimilarities for F4")+
  theme(plot.title = element_text(size = rel(1.2)),
        legend.position="none")+
  theme(axis.text.x =  element_text(angle=45, hjust = 1),
        strip.text.x= element_text(size = 13),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size=12)) +
  scale_x_discrete(limits = rev)
bc_F4_box

graph2ppt(bc_F4_box, "out/bray.curtis.F4.pptx")
