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
variables_combined <- read.csv("/upscaled.dataset.txt", sep=";")
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
BRT.model.2 <- foreach(j = 1:1000, # 1000 models per dpt variable
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
            
            # name of dpt var
            dpt.y = colnames(full.dataset[ , ! names(full.dataset) %in% indpdt.x, drop = F])
            
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
            sample.dataset = as.data.frame(full.dataset[!(is.na(full.dataset[ , dpt.y, drop = F])),])
            
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
              gbm.y = dpt.y,
              family = BRTfamily,
              tree.complexity = preliminary.selection[,"tree.complexity"],
              learning.rate = preliminary.selection[,"shrinkage"],
              bag.fraction = preliminary.selection[,"bag.fraction"],
              n.trees = inittrees,
              step.size = preliminary.selection[,"stepsize"],
              n.folds = nbrfolds,
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
final.selection = read.csv("preliminary.selection.txt", sep=";")
# beginning for l loop: for each consecutive model, check the number of parameters.
for (l in 1:length(BRT.model.2)) {
  # beginning if loop: if prameters are empty because the model did not compute, do not process.
  if(length(BRT.model.2[[l]]) < 4) { 
    BRT.model.2[l] <- list(NULL)
  } else {
    BRT.model.2[[l]] <- BRT.model.2[[l]]
    # ending if loop
  }
  # ending for l loop
}
# remove NULL elements: a.k.a. parameter combination with less than 5 consecutive models
BRT.model.2 <- BRT.model.2[!sapply(BRT.model.2,is.null)]
## average the 5 consecutive model parameters for each combination and store in hypergrid
# optimal number of trees
final.selection$optimal.trees.boot <- as.numeric(arrange(as.data.frame(rlist::list.rbind(BRT.model.2)) %>%
                                                           summarise(optimal.trees = median(as.numeric(optimal.trees), na.rm = TRUE))))
# minimal RMSE
final.selection$min.RMSE.boot <- as.numeric(arrange(as.data.frame(rlist::list.rbind(BRT.model.2)) %>%
                                                      summarise(min.RMSE = median(as.numeric(min.RMSE), na.rm = TRUE))))
# relative influence
RI.final = do.call("rbind", unlist(lapply(BRT.model.2, `[`, 1), recursive = FALSE))
# relative influence median for each dependent variable stored in list
RI.final <- RI.final %>%
  group_by(var) %>%
  summarise_each(funs(median(., na.rm = TRUE)))
# ordering decreasing by relative influence
RI.final <- RI.final [order(-RI.final$rel.inf), ]

# plot
# unlist into list
plot.prediction = list()
plot.prediction.med = list()
list.plot = unlist(lapply(BRT.model.2, `[`, 4), recursive = FALSE)
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
plot.prediction.med.final = plot.prediction.med
plot.prediction.final = plot.prediction
# prediction old
# extract all the predictions in a df
prediction.final = do.call("rbind", unlist(lapply(BRT.model.2, `[`, 2), recursive = FALSE))
# prediction new
# extract all the predictions in a df
prediction.new = do.call("rbind", unlist(lapply(BRT.model.2, `[`, 3), recursive = FALSE))

# remove df
rm(BRT.model.2, ecozones.var, ecozones, multiResultClass, cl)
# reset memory
gc(reset = T)

# prediction old
# transform into numeric
prediction.final$V1 <- as.numeric(prediction.final$V1)
# remove number behind points
prediction.final[,"V2"] <- try(gsub("\\..*","",unlist(prediction.final[,"V2"], use.names = FALSE)), TRUE) 
# transform into factor
prediction.final[,"V2"] <- as.factor(prediction.final[,"V2"])
# average predictions by bootstrap
prediction.final <- aggregate(prediction.final$V1, list(prediction.final$V2, prediction.final$V3), FUN=mean) 
# average of predictions by row while avoiding the NA
prediction.final <- data.table::setDT(prediction.final)[,list(medianr=median(x, na.rm = TRUE), 
                                                              meansr=mean(x, na.rm = TRUE), 
                                                              lowquant=quantile(x, probs=c(0.025), na.rm = TRUE),
                                                              highquant = quantile(x, probs=c(0.975), na.rm = TRUE)), by=Group.1]

# calculate CI range
prediction.final$rangeCI95 = prediction.final$highquant - prediction.final$lowquant
# rename
colnames(prediction.final) = c("V2","medianr","meansr","lowquant",  "highquant", "rangeCI95")


# transform into numeric
prediction.new$V1 <- as.numeric(prediction.new$V1)
# remove number behind points
prediction.new[,"V2"] <- try(gsub("\\..*","",unlist(prediction.new[,"V2"], use.names = FALSE)), TRUE) 
# average of predictions by row while avoiding the NA
prediction.new <- data.table::setDT(prediction.new)[,list(medianr=median(V1, na.rm = TRUE), 
                                                          meansr=mean(V1, na.rm = TRUE), 
                                                          lowquant=quantile(V1, probs=c(0.025), na.rm = TRUE),
                                                          highquant = quantile(V1, probs=c(0.975), na.rm = TRUE)), by=V2]
# reset memory
gc(reset = T)
# calculate CI range
prediction.new$rangeCI95 = prediction.new$highquant - prediction.new$lowquant
# deviance explained by the prediction median
data = cbind(na.omit(full.dataset[][order(row.names(full.dataset[])), ][,final.selection$dptvar]),
             prediction.final[order(prediction.final$V2), ])

final.selection$dev.res.boot <- dismo::calc.deviance(obs=na.omit(full.dataset[][order(row.names(full.dataset[])), ][,final.selection$dptvar]),
                                                     pred= prediction.final[order(prediction.final$V2), ]$medianr, family = "poisson", calc.mean=TRUE)
# percentage of deviance explained by the prediction median
final.selection$p.dev.exp.boot <- (1-final.selection$dev.res.boot /final.selection$dev.tot)*100

## round values
final.selection$min.RMSE.boot <- round(final.selection$min.RMSE.boot, 2)
final.selection$dev.res.boot <- round(final.selection$dev.res.boot, 2)
final.selection$p.dev.exp.boot <- round(final.selection$p.dev.exp.boot, 1)
# export
write.table(final.selection, "final.selection.txt", sep=";", row.names = F)

## relative influence table
# transverse table
species.rel.imp.sp <- tidyr::spread(RI.final, var, rel.inf)
# data.frame transform
species.rel.imp.sp = as.data.frame(species.rel.imp.sp)
# add rownames
species.rel.imp.sp$dptvar = final.selection$dptvar
# export table
write.table(species.rel.imp.sp, "species.relative.influence.txt", sep=";", row.names = F)

#--- Predictions ---#
# rowbinding
predictions.merged <- rbind(prediction.final, prediction.new)
# round numeric column
predictions.merged <- data.frame(lapply(predictions.merged,    # Using Base R functions
                                        function(x) if(is.numeric(x)) round(x, 0) else x))
# export predictions and categories
write.csv(predictions.merged, "predictions.csv", row.names = F)

## plots
# gather
df.full = do.call(rbind, plot.prediction.final)
df.sum = do.call(rbind, plot.prediction.med.final)

# merge to get the relative influence
df.sum = merge(df.sum, RI.final,  by = "var")
df.full = merge(df.full, RI.final,  by = "var")

# create new column combining name of var and rel inf
df.sum$var.rel = paste(df.sum$var, "-", round(df.sum$rel.inf, 1))
df.full$var.rel = paste(df.full$var, "-", round(df.full$rel.inf, 1))
RI.final$var.rel = paste(RI.final$var, "-", round(RI.final$rel.inf, 1))

# resort the factor for decreasing rel inf
df.full <- arrange(transform(df.full,
                             var.rel=factor(var.rel,levels=RI.final$var.rel)),var.rel)

df.sum = arrange(transform(df.sum,
                           var.rel=factor(var.rel,levels=RI.final$var.rel)),var.rel)

# plot the partial plots
tiff(filename = paste0("partial_plot/", final.selection$dptvar, "_partial_plot.tif"), width = 100, height = 50, units = "mm", res = 300)

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




# reset memory
gc(reset = T)




























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




