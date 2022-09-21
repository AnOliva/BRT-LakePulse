# Aim: Perform rarefaction on lakes present both sampled and present in upscaled dataset
# Performed in R V4.1.2

# working directory
setwd("E:/LakePulse/Chapter_2/scripts_github/1_microbes_preprocessing")

# librairies
library(vegan)
library(tidyr)
library(dplyr)

# data loading
taxonomy_16S <- read.delim("E:/LakePulse/Data_LakePulse/16S-based_metagenomics/bacteria-20210205T182937Z-001/bacteria/lp2017-2019watercolumn16s_melt.tsv")
upscaled.dataset <- read.csv("E:/LakePulse/Data_LakePulse/upscaled_dataset/upscaled.dataset.subset.txt", sep=";")

## normalization
# wide format with asv code and lake id
taxonomy_16S_wide =  spread(taxonomy_16S[, c("sample_id", "nseqs","asv_code")], asv_code, nseqs)
# transform NA into zeros
rownames(taxonomy_16S_wide) = taxonomy_16S_wide$sample_id
taxonomy_16S_wide = taxonomy_16S_wide[,-1]
taxonomy_16S_wide[is.na(taxonomy_16S_wide)] <- 0 # time consuming
# keep only lakes available in upscaled dataset
data.merged = merge(taxonomy_16S_wide, upscaled.dataset[,c("lakepulse_id"), drop = F], by.x = 0, by.y = "lakepulse_id")
rownames(data.merged) = data.merged$Row.names
data.merged = data.merged[ , ! names(data.merged) %in% c("data.merged", "Row.names")]

## Rarefaction
# rarefy using the minimum 
(raremax <- min(rowSums(data.merged)))
# set a seed
set.seed(123)
# rarefy function
Srare <- rrarefy(data.merged, raremax)

## Add all the infos associated with asv code
# transform in dataframe
Srare = as.data.frame(Srare)
# write rownames into a column id
Srare$sample_id = rownames(Srare)
# long format 
Srare_wide = gather(Srare,-length(Srare), key = "asv_code", value = "nseqs_rar")
# add full infos associated with each asv code
data.full = merge(taxonomy_16S, Srare_wide, by = c("sample_id", "asv_code"))

## Calculate relative abundance unrarefied and rarefied on the selected lakes
# rarefied relative abundance
relseqs_rar = data.full[,c("sample_id", "nseqs_rar")] %>%
  group_by(sample_id) %>%
  mutate_each(funs(./sum(.))) %>%
  cbind(clade = data.full$clade)
# rename columns
colnames(relseqs_rar) = c("sample_id", "relseqs_rar", "clade")
# unrarefied relative abundance
relseqs = data.full[,c("sample_id", "nseqs")] %>%
  group_by(sample_id) %>%
  mutate_each(funs(./sum(.))) %>%
  cbind(clade = data.full$clade)
# rename columns
colnames(relseqs) = c("sample_id", "relseqs", "clade")
# add the relative abiundances in final table
data.full = cbind(data.full,  relseqs[,"relseqs"], relseqs_rar[,"relseqs_rar"])

## reorder columns
data.full = data.full %>% select("sample_id","asv_code","sequence","kingdom","phylum","class","order","lineage","clade","tribe",
                     "nseqs","nseqs_rar","relseqs","relseqs_rar")

## export with the number of observations
write.table(data.full, file = paste0("results/rarefied_16S_", nrow(data.merged), ".txt"), sep = ";")


