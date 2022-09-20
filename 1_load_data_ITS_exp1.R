setwd("/Users/carameyer/Desktop/PhD/exp1repo/ITS")
rm(list=ls())
load("/Users/carameyer/Desktop/PhD/exp1repo/ITS_env.RData")

library(phyloseq)
library(readr)
library(dplyr)
library(phytools)

# Load raw biom file and clean data

# import biom file
tmp_ps = import_biom("Data/otu_ITS_exp1.biom")
# remove extra taxa rank
tmp_ps@tax_table <- tmp_ps@tax_table[,1:7]
# set taxa rank names
colnames(tmp_ps@tax_table) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
# import mapping file
tmp_mapping_file <- read_csv("Data/mapping_file_ITS_exp1.csv")
tmp_mapping_file<- tmp_mapping_file[1:200,]
tmp_design = sample_data(tmp_mapping_file)
sample_names(tmp_design) <- tmp_design$id

#merge onto one phyloseq object
ps_ITS_raw = merge_phyloseq(tmp_ps,tmp_design)

# final ps_ITS_raw: 199 samples & 5937 OTUs

ps_ITS = ps_ITS_raw
# remove aberrant samples
ps_ITS <- prune_taxa(taxa_sums(ps_ITS) > 0, ps_ITS) # 4096
# remove samples with seq count < 2000
ps_ITS <- prune_samples(sample_sums(ps_ITS) > 2000, ps_ITS ) # 169 samples 4096 OTUs
ps_ITS <- prune_taxa(taxa_sums(ps_ITS) > 0, ps_ITS) # 4096 OTUs


## Filter most abundant OTUs
ITS_ps = ps_ITS
# calculate OTU relative abundance
ITS_df_otu <- as.data.frame(otu_table(ITS_ps))
ITS_df_otu_freq <- apply(ITS_df_otu, 2, FUN=function(x) x/sum(x)*100)
# keep only OTUs with a relative abundance of at least 0.1%
ITS <- apply(ITS_df_otu_freq, 1, FUN=function(x) sum(x>(0.1)))
# select OTUs above frequency threshold
ITS_otus_F1.1 <- rownames(ITS_df_otu[-which(ITS==0),])
# subset selected OTUs
ITS_ps_filter1.1 <- prune_taxa(taxa_names(ITS_ps) %in% ITS_otus_F1.1, ITS_ps) # 2462 OTUs

# keep only OTUs with a relative abundance of at least 0.5%
ITS <- apply(ITS_df_otu_freq, 1, FUN=function(x) sum(x>(0.5)))
# select OTUs above frequency threshold
ITS_otus_F1.5 <- rownames(ITS_df_otu[-which(ITS==0),])
# subset selected OTUs
ITS_ps_filter1.5 <- prune_taxa(taxa_names(ITS_ps) %in% ITS_otus_F1.5, ITS_ps) # 1012 OTUs

# Now filter out based on prevelance in treatments
# /!\ takes from 30min to 3h /!\

ITS_ps = ITS_ps_filter1.1 
# calculate OTUs prevalence in treatment (ttt)
ITS_df <- psmelt(ITS_ps)
ITS_otu_prev_ttt <- data.frame(matrix(ncol=length(unique(ITS_df$treatment)),
                                      nrow=length(unique(ITS_df$OTU)), 
                                      dimnames=list(unique(ITS_df$OTU),
                                                    unique(ITS_df$treatment))), check.names = FALSE)
for (i in unique(ITS_df$OTU)) {
  for (j in unique(ITS_df$treatment)) {
    ITS_otu_prev_ttt[i,j] <- sum(ITS_df$Abundance[ITS_df$OTU == i & ITS_df$treatment == j] > 0,
                                 na.rm = T) / length(ITS_df$Sample[ITS_df$OTU == i & ITS_df$treatment == j]) *100
  }
  
} 
rm(i,j)
# calculate maximum OTUs prevalence in treatment
ITS <- apply(ITS_otu_prev_ttt,1, FUN=function(x) max(x))
# select OTUs above a minimum prevalence in treatment threshold set to 60% 
ITS_otus_F2.1 <- rownames(ITS_otu_prev_ttt[which(ITS >= 60),])
# subset selected OTUs
ITS_ps <- prune_taxa(taxa_names(ITS_ps) %in% ITS_otus_F2.1, ITS_ps)
ps_ITS_most_abund.1 = ITS_ps # 295 OTUs in 169 samples

# look at OTU tables sorted by total count
# after relative abundane filter
otu_ITS_filtered <- as.data.frame(otu_table(ITS_ps_filter1.1))
otu_ITS_filtered$total_count <- rowSums(otu_ITS_filtered)
sorted_ITS_f1.1 <- otu_ITS_filtered[order(otu_ITS_filtered$total_count),]
# after prevalence filter
otu_ITS_f2.1 <- as.data.frame(otu_table(ps_ITS_most_abund.1))
otu_ITS_f2.1$total_count <- rowSums(otu_ITS_f2.1)
sorted_ITS_f2.1 <- otu_ITS_f2.1[order(otu_ITS_f2.1$total_count),]

save.image("/Users/carameyer/Desktop/PhD/exp1repo/ITS_env.RData")

# How many replicates are left in each treatment with the 169 samples?
nrep_by_ttt <- data.frame(sample_data(ITS_ps)) %>%
  count(treatment)
write.table(nrep_by_ttt, file="out/nrep_by_treatment_ITS.csv", sep = ";", col.names = TRUE, row.names = FALSE)

