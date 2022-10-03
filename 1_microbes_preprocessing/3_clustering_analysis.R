# Aim : clustering analysis of potential pathogens
# Performed in R V4.1.2

# set working directory
setwd("/1_microbes_preprocessing/results")

# librairies
library(dplyr)
library(rcompanion)
library(GoodmanKruskal)
library("dendextend")
library(ade4)
library(ggplot2)
library("ggplotify")


## data loading
specie <- read.csv("epathogen_16S_nseqs_rar.txt", sep=";")


## normalisation of data
# rownames set to clades
rownames(specie) = specie$clade
# remove epathogen columns which are not useful
specie <- specie[,-c(1:6)]
# transpose table to get species as column name
specie <- t(specie)
# tale lakeid (rownames) and remove X and space + replace . by -
cnamesfix <- rownames(specie)
cnamesfix <- gsub("X", "", cnamesfix)
cnamesfix <- gsub("[.]", "-", cnamesfix)
rownames(specie) <- cnamesfix
# transform to matrix
specie = as.data.frame(specie)
# remove column with less than two occurrence
nbrlakes = nrow(specie)*0.05
nbrlakes_end = nrow(specie)*0.95
specie = specie[,colSums(specie>0) > nbrlakes & colSums(specie>0) < nbrlakes_end]
# transposed matrix
specie.t = as.data.frame(t(specie))


##  CRAMER V
list.cramer = list()
# for loop
for (i in 1:ncol(specie)) {
  
  # remove the column to compare
  newdata = specie#[,-i]
  
  # select the column to compare
  varcomp = specie[,i, drop = F]
  
  # create empty vector
  vec.cramer = vector()
  
  for (j in 1:ncol(newdata)) {
    
    cont.table = base::table(cbind(varcomp, newdata[,j, drop = F]))
    
    vec.cramer[j] = try(cramerV(cont.table,
                                bias.correct =  T), T)
    
    #names(vec.cramer)[j] = paste0(colnames(varcomp), "__", colnames(newdata[,j, drop = F]))
    names(vec.cramer)[j] = paste0(colnames(newdata[,j, drop = F]))
    
    
  } # j ending
  
  # store vector in list
  list.cramer[[i]] = vec.cramer
  
  # rename model list
  names(list.cramer)[i] = colnames(varcomp)
  
  
} # i ending
data.cramer = data.frame(t(sapply(list.cramer,c)))
# remove Infinity values
data.cramer[data.cramer == Inf] <- 0
# remove text column
data.cramer[grepl('Error', data.cramer)] = 0
# remove low values
#data.cramer[data.cramer < 0.30] <- NA
# remove values equal to 1 after checking that it corresponds to specie-specie
data.cramer[data.cramer == 1] <- 0
# keep colnames
rowdata.cramer = rownames(data.cramer)
# transform numeric
data.cramer = apply(data.cramer,2,as.numeric)
# set colnames
rownames(data.cramer) = rowdata.cramer
# corrplot
tiff(paste0("corrplot_cramer.tiff"), width = 300, height = 200, units = "mm", res = 200)
corrplot::corrplot(data.cramer,
                   is.corr = F,
                   method = 'square',
                   order = 'hclust',
                   type = 'lower',
                   addrect = 1,
                   diag = FALSE,
                   col = colorRampPalette(c("darkseagreen", "brown3"))(100),
                   cl.cex = 0.6,
                   tl.cex = 0.6, tl.col = "black")
dev.off()
graphics.off()


### Goodman Kruskal
# function to calculate
GKmatrix1 <- GKtauDataframe(specie, dgts = 3)
# transform max value (2) in 0
GKmatrix1[GKmatrix1 == 2] <- 0
# corrplot
tiff(paste0("corrplot_kruskal.tiff"), width = 300, height = 200, units = "mm", res = 200)
corrplot::corrplot(GKmatrix1,
                   is.corr = F,
                   method = 'square',
                   order = 'hclust',
                   type = 'lower',
                   addrect = 1,
                   diag = FALSE,
                   col = colorRampPalette(c("darkseagreen", "brown3"))(100),
                   cl.cex = 0.6,
                   tl.cex = 0.6, tl.col = "black")


dev.off()
graphics.off()


### Dissimilarity analysis
## PA matrix
specie.t.PA <- ifelse(specie.t > 0, 1, 0)
# Hellinger + euclidean
#specie.t.chi = decostand(specie.t, "hel")
#specie.t.D16 = dist(specie.t.chi)

## Sorensen
specie.t.S8 = dist.binary(specie.t.PA, method = 5)
# choose hclust method
hclust.S8 = hclust(specie.t.S8, method = "ward.D2")
# dendrogram
dend = as.dendrogram(hclust.S8)
# plot
par(mar=c(2,1,1,8), cex.axis = 0.7, cex.sub = 0.7)
dend %>% color_branches(k=4, col = c("darkorange3","brown3","dodgerblue","forestgreen"
                                     )) %>% plot(horiz=TRUE)
# add horiz rect (tests with different number of clusters)
#dend %>% rect.dendrogram(k=10,horiz=TRUE)
# add column of clusters
specie.t[ , paste0("Cluster")] <- as.vector(cutree(dend, k=4))
# plot simplified and horizontal dendrogram
tiff("cluster_histo_simp.tif", width = 260, height = 100, units = "mm", res = 300)

dend = as.dendrogram(hclust.S8)
dend2 = dend %>% set("labels_col", value = c("darkorange3","brown3","dodgerblue","forestgreen"), k=4)
par(mar=c(6,2,2,2), cex.axis = 0.7, cex.sub = 0.7)

dend %>% set("labels_col", value = c("white")) %>%
  color_branches(k=4, col = c("darkorange3","brown3","dodgerblue","forestgreen"))  %>% #,"darkgoldenrod"
  plot(horiz=F, xlab = "")

# Add axis labels rotated by 45 degrees
text(seq_along(labels(dend)), par("usr")[3] - 0, labels = labels(dend), srt = 45, adj = 1, xpd = TRUE,
     col =labels_col(dend2))

dev.off()


## sum abundance of all species within each cluster
# create a list
list.specie.t = list()
# for loop
for (i in 1:ncol(specie.t %>% dplyr::select(starts_with("Cluster")))) {
  
  # select dataset cluster associated to species
  data.clus = specie.t %>% dplyr::select(starts_with("Cluster"))
  
  # select only the cluster
  data.clus.var =  data.clus[,i]
  
  # sum all counts part of the same clusters
  data.clus.ag = aggregate(specie.t, by = list(data.clus.var), FUN = sum, na.rm = TRUE)
  
  # remove cluster column
  data.clus.ag = data.clus.ag %>% dplyr::select(-starts_with("Cluster"))
  
  # change colnames of the cluster
  colnames(data.clus.ag)[1] = colnames(data.clus)[i]
  
  # change rownames to each group cluster
  data.clus.ag[,1] = paste0(colnames(data.clus.ag)[1], ".",data.clus.ag[,1])
  rownames(data.clus.ag) = data.clus.ag[,1]
  data.clus.ag = data.clus.ag[,-1]
  
  
  # list clusters
  list.specie.t[[i]] = as.data.frame(t(data.clus.ag))
  
  # tale lakeid (rownames) and remove X and space + replace . by -
  cnamesfix <- rownames(list.specie.t[[i]])
  cnamesfix <- gsub("X", "", cnamesfix)
  cnamesfix <- gsub("[.]", "-", cnamesfix)
  rownames(list.specie.t[[i]]) <- cnamesfix
  
  
}
# table with abundances summed per cluster
cluster.table = as.data.frame(list.specie.t)


## write all species names contained in each cluster in a dataframe
# selecte cluster column with species names
data.names = dplyr::select(specie.t,contains("Cluster"))
# create a list
list.namespecies = list()
# for loop
for ( i in 1:ncol(data.names)) {
  
  # spliut the species between the clusters
  namessp = split(rownames(data.names),data.names[,i])
  
  # max number of species found per cluster
  mx <- max(lengths(namessp))
  # write clustered species in row
  namessp2 = t(data.frame(lapply(namessp, `length<-`, mx)))
  # join the species rowed by collapsing with separator
  list.namespecies[[i]] = apply(namessp2, 1, paste, collapse=" - ")
  # remove NA
  list.namespecies[[i]] = gsub("- NA", "", list.namespecies[[i]])
  # remove space between species
  list.namespecies[[i]] = gsub(" ", "", list.namespecies[[i]])
  # write a dataframe for the cluster
  list.namespecies[[i]] = as.data.frame(list.namespecies[[i]])
  # replace rownames to cluster names
  rownames(list.namespecies[[i]]) = paste0(colnames(data.names)[i], ".", gsub("X", "", rownames(list.namespecies[[i]])))
  # write cluster name in new column
  list.namespecies[[i]]$met.cluster.gr = rownames(list.namespecies[[i]])
  # change the colnames
  colnames(list.namespecies[[i]]) = c("Clades", "met.cluster.gr")
  
}
# dataframe of names
cluster.names = as.data.frame(do.call(plyr::rbind.fill, list.namespecies))
# final export
write.table(cluster.names, "clusters.names.txt", sep = ";")


## outliers set to 0 - .99 quantiles
# function to set at specific intervals
fun <- function(x){
  quantiles <- quantile( x, c(0, .99 ) )
 x[ x < quantiles[1] ] <- quantiles[1]
  x[ x > quantiles[2] ] <- quantiles[2]
  x
}

# plot histograms with and without outliers set
tiff(paste0("histogram_outliers.tiff"), width = 200, height = 100, units = "mm", res = 300)
par(mfrow=c(2,4), mar = c(2, 2, 2, 2))
hist(cluster.table$Cluster.1, main = "Non-transformed 1.1", xlab = "Number of relative abundance counts",
     breaks = sqrt(nrow(cluster.table)))
hist(cluster.table$Cluster.2, main = "Non-transformed 2.1", xlab = "Number of relative abundance counts",
     breaks = sqrt(nrow(cluster.table)))
hist(cluster.table$Cluster.3, main = "Non-transformed 3.1", xlab = "Number of relative abundance counts",
     breaks = sqrt(nrow(cluster.table)))
hist(cluster.table$Cluster.4, main = "Non-transformed 3.3", xlab = "Number of relative abundance counts",
     breaks = sqrt(nrow(cluster.table)))

hist(fun(cluster.table$Cluster.1), main = "Outliers set to 0.99", xlab = "Number of relative abundance counts",
     breaks = sqrt(nrow(cluster.table)))
hist(fun(cluster.table$Cluster.2), main = "Outliers set to 0.99", xlab = "Number of relative abundance counts",
     breaks = sqrt(nrow(cluster.table)))
hist(fun(cluster.table$Cluster.3), main = "Outliers set to 0.99", xlab = "Number of relative abundance counts",
     breaks = sqrt(nrow(cluster.table)))
hist(fun(cluster.table$Cluster.4), main = "Outliers set to 0.99", xlab = "Number of relative abundance counts",
     breaks = sqrt(nrow(cluster.table)))
dev.off()
graphics.off()
#  outlier corrections
cluster.table$Cluster.1 = round(fun(cluster.table$Cluster.1), 0)
cluster.table$Cluster.2 = round(round(fun(cluster.table$Cluster.2), 0), 0)
cluster.table$Cluster.3 = round(round(fun(cluster.table$Cluster.3), 0), 0)
cluster.table$Cluster.4 = round(round(fun(cluster.table$Cluster.4), 0), 0)
# export
write.table(cluster.table, "clusters.var.txt", sep = ";")


## frequency histogram per cluster 
# keep names of cluster.table to be removed
rem.names = names(cluster.table[ , ! names(cluster.table) %in% c( "Row.names", "idLatLong")])
# long format
cluster.table.g = tidyr::gather(cluster.table, all_of(rem.names), key = "cluster.gr", value = "rel.ab")
# transofmr to numeric
cluster.table.g$rel.ab = as.numeric(cluster.table.g$rel.ab)
# plotting histograms
p3 = ggplot(cluster.table.g, aes(x = rel.ab)) + #, fill = quantile
  geom_histogram(bins = 70, fill = "cornsilk4", color = "white") + #aes(fill = dts)binwidth = 0.001 binwidth = 0.1
  facet_wrap(cluster.gr ~ ., scale = "free", nrow = 5) + 
  xlab("Rarefied counts") +
  ylab("Number of lakes") +
theme_minimal() 
 # theme(strip.background=element_blank())

# function to change the color of the boxplots in facet_wrap
dummy <- ggplot(data = cluster.table.g, aes(x = rel.ab))+ facet_wrap(cluster.gr ~ ., scale = "free", nrow = 5) + 
  geom_rect(aes(fill=cluster.gr), xmin=-Inf, xmax=Inf,ymin=-Inf, ymax=Inf) +
  scale_fill_manual(values=c("brown3",
                             "forestgreen",
                             "darkorange3",
                             "dodgerblue")) +
  theme_minimal()

g1 <- ggplotGrob(p3)
g2 <- ggplotGrob(dummy)

gtable_select <- function (x, ...) 
{
  matches <- c(...)
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  x
}

panels <- grepl(pattern="panel", g2$layout$name)
strips <- grepl(pattern="strip-t", g2$layout$name)

g2$layout$t[panels] <- g2$layout$t[panels] - 1
g2$layout$b[panels] <- g2$layout$b[panels] - 1

new_strips <- gtable_select(g2, panels | strips)
grid::grid.newpage()
grid::grid.draw(new_strips)

gtable_stack <- function(g1, g2){
  g1$grobs <- c(g1$grobs, g2$grobs)
  g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
  g1$layout <- rbind(g1$layout, g2$layout)
  g1
}

new_plot <- gtable_stack(g1, new_strips)
grid::grid.newpage()
grid::grid.draw(new_plot)

p4 = as.ggplot(new_plot)

# function multipanel plots
lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 

# multipanel plots
tiff("cluster_histo_with0.tif", width = 200, height = 250, units = "mm", res = 300)

par(mar=c(2,1,1,8), cex.axis = 0.7, cex.sub = 0.7) # parameters for dendrogram
m <- rbind(c(1, 2), c(1, 2)) # layout for dendrogram
layout(m)

# layout for ggplot + dend
lay_out(list(dend %>% color_branches(k=4, col = c(
  "darkorange3",
  "brown3",
  "dodgerblue",
  "forestgreen"
)) %>% plot(horiz=TRUE), 1:4, 1:2),
        
        list(p4, 1:4, 3:4))

dev.off()





