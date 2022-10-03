# Aim : normalize data before BRT
# Performed in R V4.1.2

# create folder if missing
mainDir <- "/3_BRT/results"
subDir <- "aCDOM"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# set working directory
setwd(paste0(mainDir, "/",subDir))

# loading dataset
full.dataset <- read.csv("/2_BRT_preprocessing/results/QC_dataset_16S_nseqs_rar.txt", sep=";")

# remove clusters
full.dataset = full.dataset[, ! names(full.dataset) %in% c( "Cluster.1","Cluster.2", "Cluster.3", "Cluster.4")]
# names of independent variables
indpdt.x = colnames(full.dataset[, ! names(full.dataset) %in% c('aCDOM400', 'TSS')])
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