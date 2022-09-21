# Aim : check collinearity between independent variables
# Performed in R V4.1.2

# setting working directory
setwd("E:/LakePulse/Chapter_2/scripts_github/2_BRT_preprocessing/results")

## librairies
library(corrplot)
library(car)

## loading data
# combined independent variables 
variables_combined = read.csv("E:/LakePulse/Data_LakePulse/upscaled_dataset/upscaled.dataset.subset.txt", sep=";")
# species clusters
specie = read.csv("E:/LakePulse/Chapter_2/scripts_github/1_microbes_preprocessing/results/clusters.var.txt", sep=";")


### Pearson Correlatiob tests

## normalization
# keep only the lakes used for microbe preprocessing
data.merged = merge(variables_combined, specie, by.x = "lakepulse_id", by.y = 0)
# set lakeid to rownames
rownames(data.merged) = data.merged$idLatLong
# delete columns not useful
data.merged[, c("lakepulse_id","idLatLong",
                
                # geography-morphology
                "latitude","longitude","Alake","depth_avg", #"slope100m",
                
                # categories
                "province","team","AtHigh","AtMar","BorPla","BorShi",
                "MixPla", "MonCor", "PaciMar", "Prairies","SemiAPl","NorthEcoz",
                "WestMont","year", # "ecozone",
                
                # potential NPS
                "area_km2_forestry","area_km2_agriculture","area_km2_nat_landscapes",
                "area_km2_pasture","area_km2_grassland","area_km2_anthro", "area_km2_urban",
                "fraction_agriculture","fraction_nat_landscapes", "area_km2_agripas",
                "fraction_forestry","fraction_pasture","fraction_grassland", # "fraction_urban",
                "fraction_anthro","agripasNat","urbanNat","total_km2_area", # ,"fraction_agripas"
                
                "manure_wshed","animal_unit", "animal_unit_density","ecumene_agr","ecumene_pop",
                "npop", "pop_wshed",
                
                # climate-meteo
                "evaporation","precipitation","daysWithIce","fireNbr", # "Tair",
                
                # bio-optic
                "secchi_adjusted","chla", "aCDOM280","aCDOM355",#"aCDOM400", "TSS",
                
                # water quality
                "TP", "Wtemp", "salinity", "SRPNo3","TN"#,
                
                # clusters
                #"Cluster.1","Cluster.2","Cluster.3","Cluster.4"
                
)] = list(NULL)


# Pearson correlations
# create subdirectory if it does not exist
subDir <- "multicollinearity_upscaled"
ifelse(!dir.exists(file.path(subDir)), dir.create(file.path(subDir)), FALSE)

tiff(filename = paste0("multicollinearity_upscaled/multicollinearity_pearson.tif"), width = 95, height = 95, units = "mm", res = 300)
par(mar=c(0.5, 0.5, 0.5, 0.5))
corrplot(cor(data.merged[ , ! names(data.merged) %in% c("ecozone","Cluster.1","Cluster.2","Cluster.3","Cluster.4")],
             method="pearson", use="pairwise.complete.obs"),
         type="lower",
         method="color",
         outline = T,
         tl.cex = 0.8, tl.col = "black",
         addCoef.col = "black", number.digits = 2, number.cex = 0.6, 
         cl.cex = 0.8, 
         #addrect = 9, 
         #rect.lwd = 20, 
         col = colorRampPalette(c("darkcyan", "white","brown3"))(100),
         order = "hclust")
dev.off()
graphics.off()


### Checking GVIF

## data préparation and normlization
# set seed
set.seed(123)
# randomize dataset rows
random_index = sample(1:nrow(data.merged), nrow(data.merged))
data_randomized = data.merged[random_index, ] 
# environmental variables
indpt_var = data_randomized[ , ! names(data.merged) %in% c("Cluster.1","Cluster.2","Cluster.3","Cluster.4")]
# species variables
dpt_var = data_randomized[ , c("Cluster.1","Cluster.2","Cluster.3","Cluster.4")] 
# remove dpt var with less than 5% occurrence 
nbroccurrence = nrow(dpt_var)*0.05
# remove dpt var with more than 95% occurrence 
#nbroccurrence_end = nrow(dpt_var)*0.95
dpt_var = dpt_var[,colSums(dpt_var>0) > nbroccurrence ]

# create empty lists to store variables
dlist = list() 
for (i in 1:ncol(dpt_var)) {
  dlist[i] = assign(paste(colnames(dpt_var[i])),  as.vector(dpt_var[i]) )
}
# write corresponding names from df to list element
names(dlist) = colnames(dpt_var) 
# remove dataframes of genus
rm(list = colnames(dpt_var))
#list_omcdiag = list()
#list_imcdiag = list()
list_vif = list()
#list_model = list()
#list_pcor = list()

## checking GVIF

# POISSON
# for loop 
for(i in 1:ncol(dpt_var)) {
  
  indpt_var$dptvar = dlist[[i]] # extract the vector from the grid to the environmental table
  
  # glm model
  modeltest = glm(dptvar ~ ., data= indpt_var, family = poisson(link = "log"))
  # glm models can store
  #list_model[[i]] = modeltest
  # functions to test GVIF
  #list_omcdiag[[i]] = omcdiag(modeltest, na.rm = TRUE)
  #list_imcdiag[[i]] = imcdiag(modeltest, na.rm = TRUE, all = TRUE)
  list_vif[[i]] = vif(modeltest)
  #list_pcor[[i]] = pcor(na.omit(indpt_var), method = "spearman")
}
# write names of dpt variables
#names(list_omcdiag) = names(dlist)
names(list_vif) = names(dlist)
#names(list_model) = names(dlist)
#names(list_pcor) = names(dlist)
#names(list_imcdiag) = names(dlist)

# check the GVIF
P.GVIF = round(unlist(lapply(list_vif, function(x) x[,3][which.max(abs(x[,3]))])),2)

# BERNOULLI
# for loop 
for(i in 1:ncol(dpt_var)) {
  
  indpt_var$dptvar = dlist[[i]] # extract the vector from the grid to the environmental table
  #transform to P/A
  indpt_var$dptvar <- ifelse(indpt_var$dptvar > 0, 1,0)
  
  # glm model
  modeltest = glm(dptvar ~ ., data= indpt_var, family = binomial(link = "logit"))
  # glm models can store
  #list_model[[i]] = modeltest
  # functions to test GVIF
  #list_omcdiag[[i]] = omcdiag(modeltest, na.rm = TRUE)
  #list_imcdiag[[i]] = imcdiag(modeltest, na.rm = TRUE, all = TRUE)
  list_vif[[i]] = vif(modeltest)
  #list_pcor[[i]] = pcor(na.omit(indpt_var), method = "spearman")
}
# write names of dpt variables
#names(list_omcdiag) = names(dlist)
names(list_vif) = names(dlist)
#names(list_model) = names(dlist)
#names(list_pcor) = names(dlist)
#names(list_imcdiag) = names(dlist)

# check the GVIF
B.GVIF = round(unlist(lapply(list_vif, function(x) x[,3][which.max(abs(x[,3]))])),2)

# TRUNCATED POISSON
# for loop 
for(i in 1:ncol(dpt_var)) {
  
  indpt_var$dptvar = dlist[[i]] # extract the vector from the grid to the environmental table
  
  indpt_var.2 = indpt_var[indpt_var$dptvar > 0 ,]
  
  # glm model
  modeltest = glm(dptvar ~ ., data= indpt_var.2, family = poisson(link = "log"))
  # glm models can store
  #list_model[[i]] = modeltest
  # functions to test GVIF
  #list_omcdiag[[i]] = omcdiag(modeltest, na.rm = TRUE)
  #list_imcdiag[[i]] = imcdiag(modeltest, na.rm = TRUE, all = TRUE)
  list_vif[[i]] = vif(modeltest)
  #list_pcor[[i]] = pcor(na.omit(indpt_var), method = "spearman")
}
# write names of dpt variables
#names(list_omcdiag) = names(dlist)
names(list_vif) = names(dlist)
#names(list_model) = names(dlist)
#names(list_pcor) = names(dlist)
#names(list_imcdiag) = names(dlist)

# check the GVIF
TP.GVIF = round(unlist(lapply(list_vif, function(x) x[,3][which.max(abs(x[,3]))])),2)

# print VGIF
print(rbind(P.GVIF, B.GVIF, TP.GVIF))


## export table for further analysis
write.table(data.merged,
            "upscaled_dataset_16S_nseqs_rar.txt", sep=";")

