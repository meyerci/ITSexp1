library(dplyr)
library(ggplot2)
library(export)

tmp_ps = ps_ITS
# merge_samples will coerce sample data into NAs if they are not numeric
tmp_ps@sam_data$dose <- factor(tmp_ps@sam_data$dose, levels = c("Control", "1x", "10x", "100x"))
tmp_ps@sam_data$dose <- as.numeric(tmp_ps@sam_data$dose)
tmp_ps@sam_data$fraction <- factor(tmp_ps@sam_data$fraction, levels = c("Unfiltered", "F1", "F2", "F3", "F4"))
tmp_ps@sam_data$fraction <- as.numeric(tmp_ps@sam_data$fraction)

## aggregate by treatment
tmp_ps <- merge_samples(tmp_ps,"treatment")
## aggregate by taxa_rank
tmp_ps <- tax_glom(tmp_ps,taxrank = "Phylum")
## transform to relative abundance 
tmp_ps <- transform_sample_counts(tmp_ps,function(x) x / sum(x))

# extract relative abundance and taxonomy
tmp_data <- psmelt(tmp_ps)
tmp_data <- tmp_data[,c("Sample","Abundance","Phylum", "dose", "fraction")]
colnames(tmp_data) <- c("treatment","Abundance","Phylum", "dose", "fraction")

# rename phyla with < 5% abundance
tmp_data$Phylum[tmp_data$Abundance < 0.005] <- "< 0.5% abundance"

# Count phyla to set color palette
Count = length(unique(tmp_data$Phylum)); Count

col_vector<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
              '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', 
              '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', 
              '#000075', '#808080', '#ffffff', '#000000')

# Make Phyla names look nicer
tmp_data$Phylum <- gsub("p__", "", tmp_data$Phylum)

# Make factor
tmp_data$Phylum <- factor(tmp_data$Phylum)

# change sample data back to character factor and save data frame
ag_abund_phyla <- tmp_data 
ag_abund_phyla$dose <- recode(ag_abund_phyla$dose, '1' = "Control", '2' = "1x", '3' = "10x", '4' = "100x")
ag_abund_phyla$fraction <- recode(ag_abund_phyla$fraction, '1'= "Unfiltered", '2'="F1", '3'="F2", '4'="F3", '5'="F4")
ag_abund_phyla$dose <- factor(ag_abund_phyla$dose, levels = c("Sterile", "Control", "1x", "10x", "100x"))
ag_abund_phyla$fraction <- factor(ag_abund_phyla$fraction, levels = c("Sterile", "Unfiltered", "F1", "F2", "F3", "F4"))

phyla_fraction_plot <- ggplot(ag_abund_phyla, aes(x=dose, y=Abundance, fill=Phylum)) +
  facet_grid(~fraction) + 
  geom_bar(aes(), stat="identity", position="stack") + 
  labs (y= "Relative Abundance", x="Pesticide Dose") + 
  theme(axis.text.x =  element_text(angle=45, hjust = 1),
        strip.text.x= element_text(size = 13),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size=12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = col_vector)

phyla_fraction_plot

graph2ppt(phyla_fraction_plot ,file="out/Phylum_by_fraction.pptx")

phyla_dose_plot <- ggplot(ag_abund_phyla, aes(x=fraction, y=Abundance, fill=Phylum)) +
  facet_grid(~dose) + 
  geom_bar(aes(), stat="identity", position="stack") + 
  labs (y= "Relative Abundance", x="Fraction") + 
  theme(axis.text.x =  element_text(angle=45, hjust = 1),
        strip.text.x= element_text(size = 13),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size=12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = col_vector)

phyla_dose_plot

graph2ppt(phyla_dose_plot ,file="out/Phylum_by_dose.pptx")
