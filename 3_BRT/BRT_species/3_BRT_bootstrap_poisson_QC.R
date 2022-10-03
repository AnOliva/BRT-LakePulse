# Aim : Fit BRT with additional bagging method
# Performed in R V4.1.2
# method reference: http://uc-r.github.io/gbm_regression
# Script author: Anais Oliva anais.oliva@usherbrooke.ca

## librairies
library(bigmemory)
library(foreach)
library(dplyr)
library(ggplot2)
## create partial_plot folder if it does not exist
plotDir <- "partial_plot" # plot directory
# check if plot folder if present or create it
if (file.exists(plotDir)){
  "Folder plotDir already exists"
} else {
  dir.create(file.path(getwd(), plotDir))
}

## loading data
# tuned parameters
preliminary.selection <- read.csv("preliminary.selection.txt", sep=";")
# non-sampled lakes
variables_combined <- read.csv("/2_BRT_preprocessing/results/variables_combined.txt", sep=";")
# keep observations in sampled ecozones
variables_combined = variables_combined[variables_combined$ecozone == "Prairies" | 
                                          variables_combined$ecozone == "WestMont" |
                                          variables_combined$ecozone == "Boreal Plains" |
                                          variables_combined$ecozone == "Boreal Shield" |
                                          variables_combined$ecozone == "Atlantic Maritime" |
                                          variables_combined$ecozone == "Atlantic Highlands" |
                                          variables_combined$ecozone == "Mixedwood Plains" ,]
# remove lakes used in full.dataset
variables_combined = variables_combined[!(variables_combined$idLatLong %in% rownames(full.dataset)) ,]
# remove lakes with NA for idLatLong
variables_combined = variables_combined[!variables_combined$idLatLong == "" ,]
# rename rownames to idLatLong
rownames(variables_combined) <- variables_combined$idLatLong
# keep same columns as in full.dataset
variables_combined <- variables_combined[ ,  names(variables_combined) %in% indpdt.x]
# ecozones separated
ecozones.var = variables_combined[,c("ecozone"), drop = F]
variables_combined = variables_combined[, ! names(variables_combined) %in% c( "ecozone")]

## create big matrices in shared memory between cores
mat <- as.matrix(preliminary.selection[,-1]) 
row.names(mat) <- preliminary.selection[,1]
preliminary.selection <- as.big.matrix(mat, backingfile = "preliminary.selection")

mat <- as.matrix(full.dataset)
full.dataset <- as.big.matrix(mat, backingfile = "full.dataset")

mat <- as.matrix(variables_combined)
variables_combined <- as.big.matrix(mat, backingfile = "variables_combined")

# Create class which holds multiple results for each loop iteration.
# Each loop iteration populates two properties: $result1 and $result2.
# For a great tutorial on S3 classes, see: 
# http://www.cyclismo.org/tutorial/R/s3Classes.html#creating-an-s3-class
multiResultClass <- function(RI.final=NULL,
                             prediction.final=NULL,
                             # interaction.final=NULL, 
                             #performance.final=NULL,
                             prediction.new=NULL,
                             pplotlist.final = NULL,
                             varlist.final = NULL,
                             min.RMSE=NULL,
                             optimal.trees=NULL

)
{
  me <- list(
    RI.final = RI.final,
    prediction.final = prediction.final,
    #interaction.final = interaction.final,
    #performance.final = performance.final,
    prediction.new = prediction.new,
    pplotlist.final = pplotlist.final,
    varlist.final = varlist.final,
    min.RMSE = min.RMSE,
    optimal.trees = optimal.trees
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

## remove objects
rm(mat, random.index)

## reset memory before BRT
gc(reset = T)


## Loop that model BRT models using best tuned parameters with bagging
# register the number of cores to be used simultaneously
cl <- parallel::makeCluster(nbrcore)
doParallel::registerDoParallel(cl)

# timer ON for the model
ptm <- proc.time() 

# 3: beginning for l loop
BRT.model.2 <- foreach(i = 1:nrow(preliminary.selection)) %:% # nbr of dpt variables
  foreach(j = 1:1000, # 1000 models per dpt variable
          #.multicombine = T,
          .packages = c("gbm",
                        "dismo",
                        "bigmemory"
          ),
          .errorhandling = 'pass',
          .verbose = F) %dopar% {
            
            # import big.matrix
            full.dataset <- attach.big.matrix("full.dataset.desc")
            preliminary.selection <- attach.big.matrix("preliminary.selection.desc")
            variables_combined <- attach.big.matrix("variables_combined.desc")
            
            # add ecozones
            full.dataset = as.data.frame(cbind(full.dataset[,], ecozones))
            full.dataset$ecozone = as.factor(full.dataset$ecozone)
            
            variables_combined = as.data.frame(cbind(variables_combined[,], ecozones.var))
            variables_combined$ecozone = as.factor(variables_combined$ecozone)
            
            # create the class object that will store the results
            result <- multiResultClass()
            
            # family type
            BRTfamily <- "poisson"
            # initial number of trees
            inittrees <- 1
            # number of folds
            nbrfolds <- 10
            # maximum trees
            maxitrees <- 3000
            
            # sample dataset
            sample.dataset = as.data.frame(full.dataset[!(is.na(full.dataset[,rownames(preliminary.selection[,0])[i]])),])
            
            # reproducibility
            options(warn=1)
            RNGkind(sample.kind = "Rounding")
            set.seed(j) # different seeds each time
            sample.dataset = as.data.frame(sample.dataset[sample(1:nrow(sample.dataset),nrow(sample.dataset), replace = TRUE),]) 
            
            # reproducibility
            options(warn=1)
            RNGkind(sample.kind = "Rounding")
            set.seed(j) # different seeds each time
            
            # train GBM model
            BRT.final.boot <- try(gbm.step(
              # keep non-NA observations and random sampling
              data =  sample.dataset,                                                                                                                   
              gbm.x = indpdt.x,
              gbm.y = rownames(preliminary.selection[,0])[i],
              family = BRTfamily,
              tree.complexity = preliminary.selection[,"tree.complexity"][i],
              learning.rate = preliminary.selection[,"shrinkage"][i],
              bag.fraction = preliminary.selection[,"bag.fraction"][i],
              n.trees = inittrees,
              step.size = preliminary.selection[,"stepsize"][i],
              n.folds = 10,
              n.cores = nbrcore,
              max.trees = maxitrees,
              prev.stratify = T,
              silent = T,
              plot.main = F,
              keep.fold.models = F,                 # keep the fold models from cross validation
              keep.fold.vector = F,                 # allows the vector defining fold membership to be kept
              keep.fold.fit = T
            ), TRUE)
            
            
            # store results
            # minimal RMSE median over 1000 bootstraps
            result$min.RMSE <- sqrt(min(as.numeric(unlist(BRT.final.boot$cv.values), na.rm = TRUE)))
            # optimal number of trees median over 1000 bootstrap
            result$optimal.trees <- BRT.final.boot$n.trees
            # predictions
            result$prediction.final <- as.data.frame(cbind(exp(BRT.final.boot[["fold.fit"]]),rownames(BRT.final.boot[["gbm.call"]][["dataframe"]]),
                                                           V3=j))
            # predictions
            result$prediction.new <- as.data.frame(cbind(predict.gbm(BRT.final.boot, 
                                                                     as.data.frame(variables_combined[,indpdt.x, drop = F]), 
                                                                     n.trees = BRT.final.boot$gbm.call$best.trees,
                                                                     rownames(variables_combined),
                                                                     type = "response"),
                                                         rownames(variables_combined)))
            
            # relative influence
            result$RI.final <- summary(BRT.final.boot, las = 1, plotit = FALSE)
            
            # partial plot predictions
            pplotlist = list()
            varlist = list()
            
            # 11: begin for p loop
            for(p in 1: length(indpdt.x)) {
              pplotlist[[p]] <- plot.gbm(BRT.final.boot, i.var = p, n.trees = BRT.final.boot[["n.trees"]], return.grid =  TRUE) # gbm.step
              varlist[[p]] <- c(colnames(pplotlist[[p]][1]), "preds")
              # 11: end for p loop
            }
            
            result$pplotlist.final = pplotlist
            result$varlist.final = varlist
            
            
            
            # delete the model
            rm(BRT.final.boot)
            
            # return stored item
            return(result)
            
            # 3: end for l loop
          }

# stop timer
BRTmodel2.time <- proc.time() - ptm
proc.time() - ptm

# stop the cluster
parallel::stopCluster(cl)
doParallel::stopImplicitCluster()
closeAllConnections()

# reset memory
gc(reset = T)


## average the parameters
# create lists for predictions
RI.final <- list()
prediction.final <- list()
prediction.new <- list()
plot.prediction.final <- list()
plot.prediction.med.final = list()
final.selection = as.data.frame(preliminary.selection[])
# beginning for k loop: for each indpt variables
for (k in 1:length(BRT.model.2)) {
  # beginning for l loop: for each consecutive model, check the number of parameters.
  for (l in 1:length(BRT.model.2[[k]])) {
    # beginning if loop: if prameters are empty because the model did not compute, do not process.
    if(length(BRT.model.2[[k]][[l]]) < 4) { 
      BRT.model.2[[k]][l] <- list(NULL)
    } else {
      BRT.model.2[[k]][[l]] <- BRT.model.2[[k]][[l]]
      # ending if loop
    }
    # ending for l loop
  }
  # remove NULL elements: a.k.a. parameter combination with less than 5 consecutive models
  BRT.model.2[[k]] <- BRT.model.2[[k]][!sapply(BRT.model.2[[k]],is.null)]
  ## average the 5 consecutive model parameters for each combination and store in hypergrid
  # optimal number of trees
  final.selection$optimal.trees.boot[k] <- as.numeric(arrange(as.data.frame(rlist::list.rbind(BRT.model.2[[k]])) %>%
                                                                summarise(optimal.trees = median(as.numeric(optimal.trees), na.rm = TRUE))))
  # minimal RMSE
  final.selection$min.RMSE.boot[k] <- round(as.numeric(arrange(as.data.frame(rlist::list.rbind(BRT.model.2[[k]])) %>%
                                                           summarise(min.RMSE = median(as.numeric(min.RMSE), na.rm = TRUE)))),2)
  # prediction old
  # extract all the predictions in a df
  prediction.final[[k]] = do.call("rbind", unlist(lapply(BRT.model.2[[k]], `[`, 2), recursive = FALSE))
  # transform into numeric
  prediction.final[[k]]$V1 <- as.numeric(prediction.final[[k]]$V1)
  # remove number behind points
  prediction.final[[k]][,"V2"] <- try(gsub("\\..*","",unlist(prediction.final[[k]][,"V2"], use.names = FALSE)), TRUE) 
  # transform into factor
  prediction.final[[k]][,"V2"] <- as.factor(prediction.final[[k]][,"V2"])
  # average predictions by bootstrap
  prediction.final[[k]] <- aggregate(prediction.final[[k]]$V1, list(prediction.final[[k]]$V2, prediction.final[[k]]$V3), FUN=mean) 
  
  # average of predictions by row while avoiding the NA
  prediction.final[[k]] <- prediction.final[[k]]  %>%     
    group_by(Group.1) %>%
    summarize(medianr = median(x, na.rm = TRUE), meansr = mean(x, na.rm = TRUE),
              lowquant = quantile(x, probs=c(0.025), na.rm = TRUE), highquant = quantile(x, probs=c(0.975), na.rm = TRUE))
  prediction.final[[k]]$rangeCI95 = prediction.final[[k]]$highquant - prediction.final[[k]]$lowquant
  # rename
  colnames(prediction.final[[k]]) = c("V2","medianr","meansr","lowquant",  "highquant", "rangeCI95")
  # prediction new
  # extract all the predictions in a df
  prediction.new[[k]] = do.call("rbind", unlist(lapply(BRT.model.2[[k]], `[`, 3), recursive = FALSE))
  # transform into numeric
  prediction.new[[k]]$V1 <- as.numeric(prediction.new[[k]]$V1)
  # remove number behind points
  prediction.new[[k]][,"V2"] <- try(gsub("\\..*","",unlist(prediction.new[[k]][,"V2"], use.names = FALSE)), TRUE) 
  # transform into factor
  prediction.new[[k]][,"V2"] <- as.factor(prediction.new[[k]][,"V2"])
  # average of predictions by row while avoiding the NA
  prediction.new[[k]] <- prediction.new[[k]]  %>%     
    group_by(V2) %>%
    summarize(medianr = median(V1, na.rm = TRUE), meansr = mean(V1, na.rm = TRUE),
              lowquant = quantile(V1, probs=c(0.025), na.rm = TRUE), highquant = quantile(V1, probs=c(0.975), na.rm = TRUE))
  prediction.new[[k]]$rangeCI95 = prediction.new[[k]]$highquant - prediction.new[[k]]$lowquant
  # relative influence
  RI.final[[k]] = do.call("rbind", unlist(lapply(BRT.model.2[[k]], `[`, 1), recursive = FALSE))
  # transform into factor
  RI.final[[k]][["var"]] = as.factor(RI.final[[k]][["var"]])
  # relative influence median for each dependent variable stored in list
  RI.final[[k]] <- RI.final[[k]] %>%
    group_by(var) %>%
    summarise_each(funs(median(., na.rm = TRUE)))
  # ordering decreasing by relative influence
  RI.final[[k]] <- RI.final[[k]] [order(-RI.final[[k]]$rel.inf), ]
  
  # deviance explained by the prediction median
  data = cbind(na.omit(full.dataset[][order(row.names(full.dataset[])), ][,rownames(final.selection)[k]]),
               prediction.final[[k]][order(prediction.final[[k]]$V2), ])
  
  final.selection$dev.res.boot[k] <- dismo::calc.deviance(obs=na.omit(full.dataset[][order(row.names(full.dataset[])), ][,rownames(final.selection)[k]]),
                                                           pred= prediction.final[[k]][order(prediction.final[[k]]$V2), ]$medianr, family = "poisson", calc.mean=TRUE)
  # percentage of deviance explained by the prediction median
  final.selection$p.dev.exp.boot[k] <- (1-final.selection$dev.res.boot[k] /final.selection$dev.tot[k])*100
  
  # plot
  # unlist into list
  plot.prediction = list()
  plot.prediction.med = list()
  list.plot = unlist(lapply(BRT.model.2[[k]], `[`, 4), recursive = FALSE)
  # for each indpt variable, take the median
  for (w in 1:length(indpdt.x)) {
    
    plot.prediction[[w]] = do.call("rbind", unlist(lapply(list.plot, `[`, w), recursive = FALSE))
    
    plot.prediction.med[[w]]  <- plot.prediction[[w]]  %>%     
      group_by(.[[indpdt.x[w]]]) %>%
      summarize(medianr = median(y, na.rm = TRUE), 
                lowquant = quantile(y, probs=c(0.025), na.rm = TRUE), highquant = quantile(y, probs=c(0.975), na.rm = TRUE))
    
    colnames(plot.prediction.med[[w]]) <- c("x.val", "y.val","lowquant","highquant")
    plot.prediction.med[[w]]$var = indpdt.x[w]
    plot.prediction.med[[w]]$x.val = as.numeric(plot.prediction.med[[w]]$x.val)
    
    colnames(plot.prediction[[w]]) <- c("x.val.full", "y.val.full")
    plot.prediction[[w]]$var = indpdt.x[w]
    plot.prediction[[w]]$x.val.full = as.numeric(plot.prediction[[w]]$x.val.full)
    
  }
  plot.prediction.med.final[[k]] = plot.prediction.med
  plot.prediction.final[[k]] = plot.prediction
  
}
## round values
final.selection$min.RMSE.boot <- round(final.selection$min.RMSE.boot, 2)
final.selection$dev.res.boot <- round(final.selection$dev.res.boot, 2)
final.selection$p.dev.exp.boot <- round(final.selection$p.dev.exp.boot, 1)
# export
write.table(final.selection, "final.selection.txt", sep=";", row.names = TRUE)

# rename the lists with dpt names
names(RI.final) <- rownames(final.selection) 
names(prediction.final) <- rownames(final.selection) 
names(prediction.new) <- rownames(final.selection) 
names(plot.prediction.final) <- rownames(plot.prediction.final) 


## relative influence table
# extract all species into dataframe
species.rel.imp <- do.call(rbind.data.frame, RI.final) 
# remove specific strings
species.rel.imp$species <- gsub("\\.[0-9]*$","",rownames(species.rel.imp))   #gsub("\\..*","",rownames(species.rel.imp))
# transverse table
species.rel.imp.sp <- tidyr::spread(species.rel.imp, var, rel.inf)
# round numeric values
species.rel.imp.sp[,-1] = round(species.rel.imp.sp[,-1],2)
# export table
write.table(species.rel.imp.sp, "species.relative.influence.txt", sep=";", row.names = F)




#--- Predictions ---#
# adding suffix to column names 
for (j in 1:length(prediction.final)) {
  
  colnames(prediction.final[[j]]) <- c( "V2", paste(colnames(prediction.final[[j]][, -which(names(prediction.final[[j]]) == "V2")]),rownames(final.selection)[j] ,sep="_"))
  colnames(prediction.new[[j]]) <- c( "V2", paste(colnames(prediction.new[[j]][, -which(names(prediction.new[[j]]) == "V2")]),rownames(final.selection)[j] ,sep="_"))
  
}
predictions.prob <- prediction.final %>% purrr::reduce(full_join, by="V2") # extract the prediction in a dataframe
# add idlake
rownames(predictions.prob) <- predictions.prob$V2
# round numeric column
predictions.prob <- data.frame(lapply(predictions.prob,    # Using Base R functions
                                      function(x) if(is.numeric(x)) round(x, 2) else x))
#--- Predictions in probabilities for new dataset ---#
predictions.merged <- prediction.new %>% purrr::reduce(full_join, by='V2') # extract the prediction in a dataframe
# add idlake
rownames(predictions.merged) <- predictions.merged$V2
# round numeric column
predictions.merged <- data.frame(lapply(predictions.merged,    # Using Base R functions
                                        function(x) if(is.numeric(x)) round(x, 2) else x))
# rowbinding
predictions.merged <- rbind(predictions.prob, predictions.merged)
# export predictions and categories
write.csv(predictions.merged, "predictions.csv", row.names = F)

## plots
for (i in 1:nrow(preliminary.selection)){
  
  # gather
  df.full = do.call(rbind, plot.prediction.final[[i]])
  df.sum = do.call(rbind, plot.prediction.med.final[[i]])
  
  # merge to get the relative influence
  df.sum = merge(df.sum, RI.final[[i]],  by = "var")
  df.full = merge(df.full, RI.final[[i]],  by = "var")
  
  # create new column combining name of var and rel inf
  df.sum$var.rel = paste(df.sum$var, "-", round(df.sum$rel.inf, 1))
  df.full$var.rel = paste(df.full$var, "-", round(df.full$rel.inf, 1))
  RI.final[[i]]$var.rel = paste(RI.final[[i]]$var, "-", round(RI.final[[i]]$rel.inf, 1))
  
  # resort the factor for decreasing rel inf
  df.full <- arrange(transform(df.full,
                               var.rel=factor(var.rel,levels=RI.final[[i]]$var.rel)),var.rel)
  
  df.sum = arrange(transform(df.sum,
                             var.rel=factor(var.rel,levels=RI.final[[i]]$var.rel)),var.rel)
  
  # plot the partial plots
  tiff(filename = paste0("partial_plot/", rownames(preliminary.selection[,0])[i], "_partial_plot.tif"), width = 100, height = 50, units = "mm", res = 300)
  
  pplot5 <- ggplot(df.sum, aes(x = x.val, y = y.val)) +
    geom_point(size = 0.3, alpha = 0.1, col = "grey", data= df.full, aes(x = x.val.full, y = y.val.full)) +
    geom_point(data = subset(df.sum, var == "ecozone"), size = 0.3, col = "black") +
    geom_point(data = subset(df.sum, var == "ecozone"), aes(y = lowquant), col = "cornflowerblue", alpha = 0.5, size = 0.3) +
    geom_point(data = subset(df.sum, var == "ecozone"), aes(y = highquant), col = "cornflowerblue", alpha = 0.5, size = 0.3) +
    
    geom_smooth(data = subset(df.sum, var != "ecozone"),method = "loess", se = F, size = 0.5, col = "black") +
    geom_smooth(data = subset(df.sum, var != "ecozone"),method = "loess", se = F, aes(y = lowquant), col = "cornflowerblue", alpha = 0.5, size = 0.4) +
    geom_smooth(data = subset(df.sum, var != "ecozone"),method = "loess", se = F, aes(y = highquant), col = "cornflowerblue", alpha = 0.5, size = 0.4) +
    facet_wrap(~var.rel, scale = "free_x", ncol= 4)+
    #coord_cartesian(ylim = c(0, max(unlist(pplotlist))))+
    labs(x = "Independent vars",
         y = "Predictions")+
    theme_bw() +
    #ggtitle("") + #rownames(preliminary.selection[,0])[l]
    theme(panel.grid.major = element_line(color = 'black', linetype = 'dotted', size = 0.1),
          plot.title=element_text(family='', face='bold', size=6, hjust=0.5, vjust=0.5),
          axis.text=element_text(size=4),
          axis.title=element_text(size=4),
          strip.text = element_text(size = 4))
  
  print(pplot5)
  dev.off()
  graphics.off()
  
  
  
}































#### export interactions
#interaction.final.df <- do.call(rbind.data.frame, interaction.final) 
# remove specific strings
#interaction.final.df$species <- gsub("\\..*","",rownames(interaction.final.df))
# concatenate text from two cells in one
#interaction.final.df$var1.2.names <- paste(interaction.final.df$var1.names,interaction.final.df$var2.names,sep="-")
#interaction.final.df <- interaction.final.df[ , -which(names(interaction.final.df) %in% c("var1.names","var2.names"))]
# transverse table
#interaction.final.df.sp <- t(spread(interaction.final.df, var1.2.names, int.size))
#colnames(interaction.final.df.sp) <- interaction.final.df.sp[1,]
#interaction.final.df.sp <- interaction.final.df.sp[-1,]
# export result
#write.csv(interaction.final.df, "interactions_predictive_variables.csv", row.names = TRUE)


# rename the list of performances

#performance.final.df <- as.data.frame(performance.final) 
#rownames(performance.final.df) <- performance.final.df[,1]
#performance.final.df <- performance.final.df[,!grepl(".V1",colnames(performance.final.df))]
#performance.final.df %>% select(-contains(".V1"))
#names(performance.final.df) <- preliminary.selection$speciesnames 
# transform in numeric
#performance.final.df2 <- as.data.frame(lapply(performance.final.df, function(x) as.numeric(as.character(x))))
#rownames(performance.final.df2) <- rownames(performance.final.df)
# round
#performance.final.df2 <- round(performance.final.df2,2)
#export result
#write.csv(performance.final.df2, "model_bootstrap_performance.csv", row.names = TRUE)




