library(dplyr)
library(xlsx)
library(readr)

tmp_mapping_file <- read_csv("Data/mapping_file_ITS_exp1.csv")[1:200,]

tmp_1 <- read.delim("Data/alpha_rarefaction_2000_0.txt", 
                    header = TRUE, sep = "\t", dec = ".", check.names = FALSE) %>%
  rename("id" = "") %>%
  merge(tmp_mapping_file, ., by = "id")

tmp_2 <- read.delim("Data/alpha_rarefaction_2000_1.txt", 
                    header = TRUE, sep = "\t", dec = ".", check.names = FALSE) %>%
  rename("id" = "") %>%
  merge(tmp_mapping_file, ., by = "id")

tmp_3 <- read.delim("Data/alpha_rarefaction_2000_2.txt", 
                    header = TRUE, sep = "\t", dec = ".", check.names = FALSE) %>%
  rename("id" = "") %>%
  merge(tmp_mapping_file, ., by = "id")

tmp_4 <- read.delim("Data/alpha_rarefaction_2000_3.txt", 
                    header = TRUE, sep = "\t", dec = ".", check.names = FALSE) %>%
  rename("id" = "") %>%
  merge(tmp_mapping_file, ., by = "id")

tmp_5 <- read.delim("Data/alpha_rarefaction_2000_4.txt", 
                    header = TRUE, sep = "\t", dec = ".", check.names = FALSE) %>%
  rename("id" = "") %>%
  merge(tmp_mapping_file, ., by = "id")

tmp_6 <- read.delim("Data/alpha_rarefaction_2000_5.txt", 
                    header = TRUE, sep = "\t", dec = ".", check.names = FALSE) %>%
  rename("id" = "") %>%
  merge(tmp_mapping_file, ., by = "id")

tmp_7 <- read.delim("Data/alpha_rarefaction_2000_6.txt", 
                    header = TRUE, sep = "\t", dec = ".", check.names = FALSE) %>%
  rename("id" = "") %>%
  merge(tmp_mapping_file, ., by = "id")

tmp_8 <- read.delim("Data/alpha_rarefaction_2000_7.txt", 
                    header = TRUE, sep = "\t", dec = ".", check.names = FALSE) %>%
  rename("id" = "") %>%
  merge(tmp_mapping_file, ., by = "id")

tmp_9 <- read.delim("Data/alpha_rarefaction_2000_8.txt", 
                    header = TRUE, sep = "\t", dec = ".", check.names = FALSE) %>%
  rename("id" = "") %>%
  merge(tmp_mapping_file, ., by = "id")

tmp_10 <- read.delim("Data/alpha_rarefaction_2000_9.txt", 
                     header = TRUE, sep = "\t", dec = ".", check.names = FALSE) %>%
  rename("id" = "") %>%
  merge(tmp_mapping_file, ., by = "id")

tmp_all <- bind_rows("1" = tmp_1, "2" = tmp_2, "3" = tmp_3, "4" = tmp_4,
                     "5" = tmp_5, "6" = tmp_6, "7" = tmp_7, "8" = tmp_8,
                     "9" = tmp_9, "10" = tmp_10, .id = "iteration") %>% 
  select(2,1,3:12) 

tmp_all$iteration <- as.numeric(tmp_all$iteration) 

tmp_all <- arrange(tmp_all, id, iteration)

tmp_mean <- tmp_all %>% 
  group_by(id) %>%  
  summarise_if(is.numeric, mean) %>% #get mean of each index for each sample
  merge(tmp_mapping_file, ., by = "id") %>% # add mapping data again (it was dropped in last argument)
  select(-starts_with("iter"))

alpha_div_ITS_exp1 <- tmp_mean

write.csv(alpha_div_ITS_exp1, file = "out/alpha_div_ITS_exp1.csv")

#write.xlsx(alpha_div_16S_exp1, file = "Data/alpha_div_16S_exp1.xlsx",
#   sheetName = "alpha_div_16S_exp1")

rm(list = names(.GlobalEnv)[grep("tmp",names(.GlobalEnv))])

save.image("/Users/carameyer/Desktop/PhD/exp1repo/ITS_env.RData")
