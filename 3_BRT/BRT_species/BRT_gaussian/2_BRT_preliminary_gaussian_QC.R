# Aim : BRT (5 consecutive models) to select best tuning parameters
# Performed in R V4.1.2

# libraries
library(bigmemory)
library(foreach)
library(dplyr)
# load normalization script
source("E:/LakePulse/ARGs/BRT/3_BRT/BRT_gaussian/1_normalization_gaussian_QC.R")

## parameters of hypergrid to be tuned
depdt.name <- colnames(full.dataset[ , ! names(full.dataset) %in% indpdt.x, drop = F])
lr = .03
tc = c(1:3)
bf = c(.7, .8, .9, 1)
ss = c(10,20,30,40)

## Create a hypergrid of parameters to be tested in BRT model
hyper.grid <- expand.grid(
  dptvar = depdt.name, 
  shrinkage = lr,
  tree.complexity = tc,
  bag.fraction = bf,
  stepsize = ss,
  optimal.trees = 0,               
  min.RMSE = 0,
  #Time = NA ,
  dev.tot = 0,
  dev.res = 0,
  dev.exp = 0
)
# transform columns in character format
hyper.grid$dptvar <- as.character(hyper.grid$dptvar)
# total number of combinations
(nbrcombi <- nrow(hyper.grid)) 

## create big matrices in shared memory between cores
mat <- as.matrix(hyper.grid[,c(2:5)])
rownames(mat) <- hyper.grid[,1]
hgm <- as.big.matrix(mat, backingfile = "hgm")

mat <- as.matrix(full.dataset)
rownames(mat) <- rownames(full.dataset)
full.dataset <- as.big.matrix(mat, backingfile = "full.dataset")

# Create class which holds multiple results for each loop iteration.
# Each loop iteration populates two properties: $result1 and $result2.
# For a great tutorial on S3 classes, see: 
# http://www.cyclismo.org/tutorial/R/s3Classes.html#creating-an-s3-class
multiResultClass <- function(min.RMSE=NULL, optimal.trees=NULL,dev.tot=NULL, dev.res=NULL, nbrobs=NULL)
{
  me <- list(
    min.RMSE = min.RMSE,
    optimal.trees = optimal.trees,
    dev.tot = dev.tot,
    dev.res = dev.res,
    nbrobs = nbrobs
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

# remove unuseful objects
rm(random.index, bf, ss, tc, lr, depdt.name, mat)
# reset memory before BRT
gc(reset = T)



## Loop that tests each dependent variables with each possible combination set in the hyper.grid
# it uses the nbr of cores specified in element "nbrcore"

# register the number of cores to be used simultaneously
cl <- parallel::makeCluster(nbrcore)
doParallel::registerDoParallel(cl)

# timer ON for the model
ptm <- proc.time()

# 1: beginning for i loop
BRT.model1.pre <- foreach(i = 1:nbrcombi) %:% # nbr of tests in hypergrid
  foreach(j = 1:5, # 5 models per hypergrid test
          #.multicombine = T,
          .packages = c("dismo",
                        "bigmemory"),
          .errorhandling = 'pass',
          .verbose = F) %dopar% {
            
            # attach the big matrices
            hgm <- attach.big.matrix("hgm.desc")
            full.dataset <- attach.big.matrix("full.dataset.desc")
            
            # add ecozones
            full.dataset = as.data.frame(cbind(full.dataset[,], ecozones))
            full.dataset$ecozone = as.factor(full.dataset$ecozone)
            
            # create the class object that will store the results
            result <- multiResultClass()
            
            # family type
            BRTfamily <- "gaussian"
            # initial number of trees
            inittrees <- 1
            # number of folds
            nbrfolds <- 10
            # maximum trees
            maxitrees <- 3000
            
            # reproducibility
            options(warn=1)
            RNGkind(sample.kind = "Rounding")
            set.seed(j) # different seeds each time
            
            # train GBM model
            BRT.preliminary <- try(gbm.step( 
              data = as.data.frame(full.dataset[!(is.na(full.dataset[,rownames(hgm[,0])[i]])), ]), # select rows that are not NA for each dependent variable
              gbm.x = indpdt.x,
              gbm.y = rownames(hgm[,0])[i],
              family = BRTfamily,
              tree.complexity = hgm[,"tree.complexity"][i],
              learning.rate = hgm[,"shrinkage"][i],
              bag.fraction = hgm[,"bag.fraction"][i],
              n.trees = inittrees,
              step.size = hgm[,"stepsize"][i],
              n.folds = nbrfolds,
              n.cores = nbrcore,
              max.trees = maxitrees,
              prev.stratify = T,
              silent = T,
              plot.main = F,
              keep.fold.models = F,                 # keep the fold models from cross valiation
              keep.fold.vector = F,                 # allows the vector defining fold membership to be kept
              keep.fold.fit = T         # allows the predicted values for observations from CV to be kept  
            ), TRUE)
            
            # store the results
            # minimal RMSE found
            result$min.RMSE <- sqrt(min(as.numeric(unlist(BRT.preliminary$cv.values), na.rm = TRUE)))
            # the optimal nbr of trees
            result$optimal.trees <- BRT.preliminary$n.trees
            # total initial deviance
            result$dev.tot <- BRT.preliminary$self.statistics$mean.null
            # redisual deviance
            result$dev.res <- BRT.preliminary$cv.statistics$deviance.mean
            # nbr of observations used
            result$nbrobs <- nrow(as.data.frame(full.dataset[!(is.na(full.dataset[,rownames(hgm[,0])[i]])), ]))
            # return the result
            return(result)
            
            # delete the model
            rm(BRT.preliminary)
            
            # 1: end for i loop 
          }

BRTmodel1.time <- proc.time() - ptm
proc.time() - ptm # ti

# stop the simultaneous use of cores
parallel::stopCluster(cl)
doParallel::stopImplicitCluster()
closeAllConnections()

# reset memory
gc(reset = T)


## retrieve the 5 consecutive model for each parameter combination
# beginning for k loop: for each parameter combination.
for (k in 1:length(BRT.model1.pre)) {
  # beginning for l loop: for each consecutive model, check the number of parameters.
  for (l in 1:length(BRT.model1.pre[[k]])) {
    # beginning if loop: if prameters are empty because the model did not compute, do not process.
    if(length(BRT.model1.pre[[k]][[l]]) < 4) { 
      BRT.model1.pre[[k]][l] <- list(NULL)
    } else {
      BRT.model1.pre[[k]][[l]] <- BRT.model1.pre[[k]][[l]]
      # ending if loop
    }
    # ending for l loop
  }
  # remove NULL elements: a.k.a. parameter combination with less than 5 consecutive models
  BRT.model1.pre[[k]] <- BRT.model1.pre[[k]][!sapply(BRT.model1.pre[[k]],is.null)]
  
  ## average the 5 consecutive model parameters for each combination and store in hypergrid
  if (length(BRT.model1.pre[[k]]) > 3) {
    
    # optimal number of trees
    hyper.grid$optimal.trees[k] <- arrange(as.data.frame(rlist::list.rbind(BRT.model1.pre[[k]])) %>%
                                                   summarise(optimal.trees = mean(as.numeric(optimal.trees), na.rm = TRUE)))
    # minimal RMSE
    hyper.grid$min.RMSE[k] <- arrange(as.data.frame(rlist::list.rbind(BRT.model1.pre[[k]])) %>%
                                              summarise(min.RMSE = mean(as.numeric(min.RMSE), na.rm = TRUE)))
    # total deviance
    hyper.grid$dev.tot[k] <- arrange(as.data.frame(rlist::list.rbind(BRT.model1.pre[[k]])) %>%
                                             summarise(dev.tot = mean(as.numeric(dev.tot), na.rm = TRUE)))
    # residual deviance
    hyper.grid$dev.res[k] <- arrange(as.data.frame(rlist::list.rbind(BRT.model1.pre[[k]])) %>%
                                             summarise(dev.res = mean(as.numeric(dev.res), na.rm = T)))
    # explained deviance
    hyper.grid$dev.exp[k] <- as.numeric(hyper.grid$dev.tot[k]) - as.numeric(hyper.grid$dev.res[k])
    
    # nbr of observations
    hyper.grid$nbrobs[k] <- arrange(as.data.frame(rlist::list.rbind(BRT.model1.pre[[k]])) %>%
                                        summarise(nbrobs = paste0(unique(nbrobs)[1])))
    
    
  } else {
    hyper.grid$optimal.trees[k] <- NA
    hyper.grid$min.RMSE[k] <- NA
    hyper.grid$dev.tot[k] <- NA
    hyper.grid$dev.res[k] <- NA
    hyper.grid$dev.exp[k] <- NA
    hyper.grid$nbrobs[k] <- NA
  }
  
}

## transform columns into numeric
hyper.grid$optimal.trees <- as.numeric(hyper.grid$optimal.trees)
hyper.grid$min.RMSE <- as.numeric(hyper.grid$min.RMSE)
hyper.grid$dev.tot <- as.numeric(hyper.grid$dev.tot)
hyper.grid$dev.res <- as.numeric(hyper.grid$dev.res)
hyper.grid$nbrobs <- as.numeric(hyper.grid$nbrobs)

## compute percentage of explained deviance
hyper.grid$p.dev.exp <- round((as.numeric(hyper.grid$dev.exp) / as.numeric(hyper.grid$dev.tot) * 100), 1)

## remove NA rows in hypergrid if there is any
hyper.grid.preliminary <- hyper.grid[!is.na(hyper.grid$min.RMSE), ]

## select the smallest min.RMSE by species and store them in a new table
preliminary.selection <- hyper.grid.preliminary[ hyper.grid.preliminary$min.RMSE == ave(hyper.grid.preliminary$min.RMSE, hyper.grid.preliminary$dptvar, FUN=min), ] 

## round values
preliminary.selection$optimal.trees <- round(preliminary.selection$optimal.trees, 0)
preliminary.selection$min.RMSE <- round(preliminary.selection$min.RMSE, 2)
preliminary.selection$dev.tot <- round(preliminary.selection$dev.tot, 6)
preliminary.selection$dev.res <- round(preliminary.selection$dev.res, 6)
preliminary.selection$dev.exp <- round(preliminary.selection$dev.exp, 3)

## export hypergrid and tuned parameter table
write.table(hyper.grid.preliminary, "hyper.grid.preliminary.txt", sep=";", row.names = FALSE) # print the hypergrid and remove the vectors in column 1
write.table(preliminary.selection, "preliminary.selection.txt", sep=";", row.names = FALSE) # save a copy



