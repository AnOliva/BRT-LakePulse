# Aim : normalize data before BRT
# Performed in R V4.1.2

# create folder if missing
mainDir <- "E:/LakePulse/ARGs/BRT/3_BRT/results"
subDir <- "QC_bernoulli"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# set working directory
setwd(paste0(mainDir, "/",subDir))

# loading dataset
full.dataset <- read.csv("E:/LakePulse/ARGs/BRT/2_BRT_preprocessing/results/QC_dataset.txt", sep=";")

# names of independent variables (removing clusters)
indpdt.x = colnames(full.dataset[, ! names(full.dataset) %in% c( "Cluster.1","Cluster.2", "Cluster.3", "Cluster.4")])
# transform in P/A
full.dataset[ , ! names(full.dataset) %in% c(indpdt.x, "lakepulse_id")] = ifelse(full.dataset[ , ! names(full.dataset) %in% c(indpdt.x, "lakepulse_id")] > 0 , 1, 0)
# remove clusters occurent above 95% of the dataset
nbrlakes = nrow(full.dataset)*0.05
nbrlakes_end = nrow(full.dataset)*0.95
colSums(full.dataset[,c( "Cluster.1","Cluster.2", "Cluster.3", "Cluster.4")]>0) 
full.dataset = full.dataset[, ! names(full.dataset) %in% c( "Cluster.1","Cluster.2")]
rm(nbrlakes, nbrlakes_end)
# reproducibility
options(warn=1)
RNGkind(sample.kind = "Rounding")
set.seed(123)
# randomise the rows
random.index <- sample(1:nrow(full.dataset), nrow(full.dataset))
full.dataset <- full.dataset[random.index, ]

# ecozones separated
ecozones = full.dataset[,c("ecozone"), drop = F]
full.dataset = full.dataset[, ! names(full.dataset) %in% c( "ecozone")]

# number of cores
nbrcore = 5