# Aim : normalize data before BRT
# Performed in R V4.1.2

# create folder if missing
mainDir <- "E:/LakePulse/Chapter_2/scripts_github/3_BRT/results"
subDir <- "upscaled_poisson"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# set working directory
setwd(paste0(mainDir, "/",subDir))

# loading dataset
full.dataset <- read.csv("E:/LakePulse/Chapter_2/scripts_github/2_BRT_preprocessing/results/upscaled_dataset_16S_nseqs_rar.txt", sep=";")

# reorder by rownames (ATTENTION: important to keep same results as in old dataset)
full.dataset = full.dataset[order(rownames(full.dataset)), ] 
# keep the cluster required for analysis
full.dataset = full.dataset[,colnames(full.dataset[ , ! names(full.dataset) %in% c("Cluster.1","Cluster.2","Cluster.4")] )] # out: "Cluster.3"

# names of independent variables (removing clusters)
indpdt.x = colnames(full.dataset[, ! names(full.dataset) %in% c( "Cluster.1","Cluster.2", "Cluster.3", "Cluster.4")])
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
nbrcore = 6