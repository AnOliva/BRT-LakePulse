# Aim : Compare predictions with observations
# and build frequency table per mapping category
# Performed in R V4.1.2
# method reference: http://uc-r.github.io/gbm_regression
# Script author: Anais Oliva anais.oliva@usherbrooke.ca

# set working directory
setwd("E:/LakePulse/Chapter_2/3_BRT/results")

# libraries
library(tidyr)
library(dplyr)
library(ggplot2)

# observations
obs = read.csv("E:/LakePulse/Chapter_2/2_BRT_preprocessing/results/QC_dataset_16S_nseqs_rar.txt", sep=";")
# set idLatLong to column
obs$idLatLong = rownames(obs)
# gather obs
obs = obs[,c("idLatLong", "Cluster.3"), drop = F]
# rename columns
colnames(obs) = c("idLatLong","observed")

# predictions QC
preds.QC <- read.csv("E:/LakePulse/Chapter_2/3_BRT/results/QC_poisson/predictions.csv")
# select the cluster that was upscaled
preds.QC = preds.QC[,c("V2","medianr_Cluster.3")]
# rename columns
colnames(preds.QC) = c("idLatLong","predicted")

# predictions upscaled
preds.upscaled <- read.csv("E:/LakePulse/Chapter_2/3_BRT/results/upscaled_poisson/predictions.csv")
# select median
preds.upscaled = preds.upscaled[,c("V2","medianr")]
# rename columns
colnames(preds.upscaled) = c("idLatLong","predicted")

# classes of models
preds.QC$mod.type = "Poisson-QC"
preds.upscaled$mod.type = "Poisson-upscaled"

# combine predicted data
data.combined = plyr::rbind.fill(preds.QC, preds.upscaled)

# combine with obs
data.combined = merge(data.combined, 
                  obs,
                  by = c("idLatLong"))

# calculate residuals
data.combined$resid = data.combined$observed - data.combined$predicted
data.combined$resid.std = (data.combined$resid) / sd(data.combined$resid)

# rounding
data.combined$predicted = round(data.combined$predicted,0)


# coloring by category of standardized residuals
data.combined$resid.std.cat = as.factor(ifelse(data.combined$resid.std >= 1.96, "> 1.96",
                                           ifelse(data.combined$resid.std < 1.96 & data.combined$resid.std > -1.96, "-1.96 to 1.96",
                                                  ifelse(data.combined$resid.std < -1.96, "< -1.96", NA))))

# transform to factor and set levels
data.combined$mod.type = as.factor(data.combined$mod.type)
data.combined$mod.type <- factor(data.combined$mod.type, levels = c("Poisson-QC", "Poisson-upscaled"))
data.combined$resid.std.cat <- factor(data.combined$resid.std.cat, levels = c("> 1.96", "-1.96 to 1.96", "< -1.96"))

# logged
data.combined$observed.log = log10(data.combined$observed+1)
data.combined$predicted.log =  log10(data.combined$predicted+1)
# residuals
data.combined$resid.log = data.combined$observed.log - data.combined$predicted.log
data.combined$resid.std.log = (data.combined$resid.log) / sd(data.combined$resid.log)

data.combined$resid.std.cat.log = as.factor(ifelse(data.combined$resid.std.log >= 1.96, "> 1.96",
                                               ifelse(data.combined$resid.std.log < 1.96 & data.combined$resid.std.log > -1.96, "-1.96 to 1.96",
                                                      ifelse(data.combined$resid.std.log < -1.96, "< -1.96", NA))))
data.combined$resid.std.cat.log <- factor(data.combined$resid.std.cat.log, levels = c("> 1.96", "-1.96 to 1.96", "< -1.96"))

# plotting
tiff(filename = paste0("observed_vs_predicted_R2.tif"), width = 120, height = 70, units = "mm", res = 200)
ggplot(data.combined, aes(predicted, observed, col=as.factor(resid.std.cat))) +
  geom_point() +
  facet_wrap(. ~  mod.type, scale = "free", nrow = 3, ncol = 3)+ 
  geom_abline() +
  ggpubr::stat_cor(aes(label = ..rr.label..), color = "black", label.x = 30, label.y = 10)+
  theme_bw() +
  guides(color = guide_legend(title = "Standardized residuals"), shape = guide_legend(title = "Binary classification")) +
  theme(legend.position="bottom")
dev.off()
graphics.off()


tiff(filename = paste0("observed_vs_predicted_R2_log.tif"), width = 120, height = 70, units = "mm", res = 200)
ggplot(data.combined, aes(predicted.log, observed.log, col=as.factor(resid.std.cat.log))) +
  geom_point() +
  facet_wrap(. ~  mod.type, scale = "free", nrow = 3, ncol = 3)+ 
  geom_abline() +
  ggpubr::stat_cor(aes(label = ..rr.label..), color = "black",  label.x = 1.1, label.y = 0.7)+
  theme_bw() +
  guides(color = guide_legend(title = "Standardized residuals"), shape = guide_legend(title = "Binary classification")) +
  theme(legend.position="bottom")
dev.off()
graphics.off()

## extract subsets
write.table(data.combined[data.combined$mod.type == "Poisson-QC" ,], "data.combined.QC.csv", sep = ",", row.names =  F)
write.table(data.combined[data.combined$mod.type == "Poisson-upscaled" ,], "data.combined.upscaled.sub.csv", sep = ",", row.names =  F)

## extract upscaled
preds.upscaled <- read.csv("E:/LakePulse/Chapter_2/3_BRT/results/upscaled_poisson/predictions.csv")
preds.upscaled =  preds.upscaled[,c("V2", "medianr","rangeCI95")]
# renames columns
colnames(preds.upscaled) = c("idLatLong", "predictions","rangeCI95")
# export
write.table(preds.upscaled , "data.upscaled.csv", sep = ",", row.names =  F)



## table of frequency per ecozones
upscaled.dataset <- read.csv("E:/LakePulse/Data_LakePulse/upscaled_dataset/upscaled.dataset.txt", sep=";")

# merging
upscaled.data = merge(preds.upscaled, upscaled.dataset[,c("idLatLong", "ecozone")],by.x = "idLatLong", by.y = "idLatLong")
QC.data = merge(data.combined[data.combined$mod.type == "Poisson-QC" ,], upscaled.dataset[,c("idLatLong", "ecozone")],by.x = "idLatLong", by.y = "idLatLong")
upscaled.data.sub = merge(data.combined[data.combined$mod.type == "Poisson-upscaled" ,], upscaled.dataset[,c("idLatLong", "ecozone")],by.x = "idLatLong", by.y = "idLatLong")

# factor
upscaled.data$ecozone = as.factor(upscaled.data$ecozone)
QC.data$ecozone = as.factor(QC.data$ecozone)
upscaled.data.sub$ecozone = as.factor(upscaled.data.sub$ecozone)

# frequency table: count

freq.obs = QC.data %>%
  group_by(ecozone) %>%
  summarise("0-2" = sum(observed <= 2),
            "3-7" = sum(observed > 2 & observed <= 7),
            ">7" = sum(observed > 7 ))

freq.obs$total = rowSums(freq.obs[,-1])


freq.QC = QC.data %>%
  group_by(ecozone) %>%
  summarise("0-2" = sum(predicted <= 2),
            "3-7" = sum(predicted > 2 & predicted <= 7),
            ">7" = sum(predicted > 7 ))

freq.QC$total = rowSums(freq.QC[,-1])


freq.upscaled.sub = upscaled.data.sub %>%
  group_by(ecozone) %>%
  summarise("0-2" = sum(predicted <= 2),
            "3-7" = sum(predicted > 2 & predicted <= 7),
            ">7" = sum(predicted > 7 ))

freq.upscaled.sub$total = rowSums(freq.upscaled.sub[,-1])


freq.upscaled = upscaled.data %>%
  group_by(ecozone) %>%
  summarise("0-2" = sum(predictions <= 2),
            "3-7" = sum(predictions > 2 & predictions <= 7),
            ">7" = sum(predictions > 7 ))
  
freq.upscaled$total = rowSums(freq.upscaled[,-1])
  
 
# frequency table: perc
freq.obs[,c(2:4)] = round(freq.obs[,c(2:4)]/freq.obs$total*100,1)
freq.QC[,c(2:4)] = round(freq.QC[,c(2:4)]/freq.QC$total*100,1)
freq.upscaled.sub[,c(2:4)] = round(freq.upscaled.sub[,c(2:4)]/freq.upscaled.sub$total*100,1)
freq.upscaled[,c(2:4)] = round(freq.upscaled[,c(2:4)]/freq.upscaled$total*100,1)


#export
write.table(freq.obs, "freq.perc.obs.txt", row.names = F, sep = ";")
write.table(freq.QC, "freq.perc.preds.QC.txt", row.names = F, sep = ";")
write.table(freq.upscaled.sub, "freq.perc.preds.ups.txt", row.names = F, sep = ";")
write.table(freq.upscaled, "freq.perc.ups.txt", row.names = F, sep = ";")





