# Aim: Extract potential pathogens from rarefied dataset based on partial match algorithm
# between taxonomical names and the ePathogen (and CFIA) Canadian database names.
# Performed in R V4.1.2


# working directory
setwd("E:/LakePulse/Chapter_2/scripts_github/1_microbes_preprocessing")

# librairies
library(dplyr)
library(tidyr)

# data loading
data.taxon <- read.csv("E:/LakePulse/Chapter_2/scripts_github/1_microbes_preprocessing/results/rarefied_16S_413.txt", sep=";") # rarefied data
epathogen <- read.csv("E:/LakePulse/Data_pathogen/epathogen-result.csv") # ePathogen
#CFIA <- read_excel("E:/LakePulse/Data_pathogen/CFIA/List_of_pests_regulated_by_Canada.xlsx") # CFIA


## normalization
# keep only species from RG2 classification and above (for human or animal) 
epathogen <- subset(epathogen, Human.classification == "RG2"
                            |Human.classification == "RG3"
                            |Human.classification == "RG4"
                            |Animal.classification == "RG2"
                            |Animal.classification == "RG3"
                            |Animal.classification == "RG4")
# keep only bacterial pathogens to be matched with 16S results
epathogen <- subset(epathogen, Agent.type == "Bacteria")
# normalize the names of ePathogen
epathogen$Name <- gsub("[0-9]+", "",epathogen$Name)
epathogen$Name <- gsub(" spp.", "", epathogen$Name)
epathogen$Name <- gsub("_", " ", epathogen$Name) 
epathogen$Name <- gsub("-", " ", epathogen$Name)
epathogen$Name <- gsub(" Excluding.*","",epathogen$Name)

## if CFIA is included
# remove disease informations from CFIA
#CFIA = CFIA[,c("Pest type",  "Genus")]
# rename columns
#colnames(CFIA) = c("Agent.type","Name")
# merge CFIA and epathogen databases
#epathogen = full_join(epathogen, CFIA[,c("Name"), drop = F], by = "Name")

# remove empty row without "name"
epathogen = epathogen[!apply(is.na(epathogen), 1, all), ]  
# remove duplicates to gain processing time
data.taxon_sub <- data.taxon[!duplicated(data.taxon[ , c("clade")]), ]
epathogen <- epathogen[!duplicated(epathogen[ , c("Name")]), ]


## Partial match functions from:
# https://www.r-bloggers.com/2012/09/merging-data-sets-based-on-partially-matched-data-elements/

signature=function(x){
  sig=paste(unlist(strsplit(tolower(x)," ")),collapse='')
  return(sig)
}

partialMatch=function(x,y,levDist=0.1){
  xx=data.frame(sig=sapply(x, signature),row.names=NULL)
  yy=data.frame(sig=sapply(y, signature),row.names=NULL)
  xx$raw=x
  yy$raw=y
  xx=subset(xx,subset=(sig!=''))
  xy=merge(xx,yy,by='sig',all=T)
  matched=subset(xy,subset=(!(is.na(raw.x)) & !(is.na(raw.y))))
  matched$pass="Duplicate"
  todo=subset(xy,subset=(is.na(raw.y)),select=c(sig,raw.x))
  colnames(todo)=c('sig','raw')
  todo$partials= as.character(sapply(todo$sig, agrep, yy$sig,max.distance = levDist,value=T))
  todo=merge(todo,yy,by.x='partials',by.y='sig')
  partial.matched=subset(todo,subset=(!(is.na(raw.x)) & !(is.na(raw.y))),select=c("sig","raw.x","raw.y"))
  partial.matched$pass="Partial"
  matched=rbind(matched,partial.matched)
  un.matched=subset(todo,subset=(is.na(raw.x)),select=c("sig","raw.x","raw.y"))
  if (nrow(un.matched)>0){
    un.matched$pass="Unmatched"
    matched=rbind(matched,un.matched)
  }
  matched=subset(matched,select=c("raw.x","raw.y","pass"))
  
  return(matched)
}



## partial matches
matches = partialMatch(epathogen$Name, na.omit(data.taxon_sub$clade))

## ATTENTION: Manual check of the partial match.
matches <- rbind(matches[!(matches$pass=="Partial" ),],
                 matches[matches$raw.x == "Burkholderia",],
                 matches[matches$raw.x == "Escherichia",],
                 matches[matches$raw.x == "Hafnia",],
                 matches[matches$raw.x == "Methylobacterium",])

# keep the potential pathogens from the taxonomical dataset based on matches
epathogen.subset <- data.taxon[data.taxon$clade %in% matches$raw.y,]

# nbr of different ASV in total
(nbr.asv = sum(length(unique(data.taxon$asv_code))))
# nbr of different potentially pathogenic ASV 
(nbr.asv.pp = sum(length(unique(epathogen.subset$asv_code))))

## melt by row (= sample_id or lake) at the clade taxonomical level (ignoring sequence, asv code and tribe)
epathogen.subset <- epathogen.subset[ , !(names(epathogen.subset) %in% c("sequence", "tribe", "asv_code"))] %>%
  group_by(sample_id, kingdom, phylum, class,
           order, lineage, clade) %>%
  summarise(across(everything(), sum))

## replace the name of the taxonomy by epathogen
epathogen.subset[epathogen.subset$clade == "Methylobacterium-Methylorubrum",]$clade <- "Methylobacterium"
epathogen.subset[epathogen.subset$clade == "Escherichia-Shigella",]$clade <- "Escherichia"
epathogen.subset[epathogen.subset$clade == "Burkholderia-Caballeronia-Paraburkholderia",]$clade <- "Burkholderia"
epathogen.subset[epathogen.subset$clade == "Hafnia-Obesumbacterium",]$clade <- "Hafnia"

## create wide format data for each type of data: nseqs, nseqs_rqr, relseqs, relseqs_rar
epathogen.subset_nseqs <- spread(epathogen.subset[, !(names(epathogen.subset) %in% c("relseqs",
                                                                                     "relseqs_rar",
                                                                                     "nseqs_rar"))], 
                                 key ="sample_id", 
                                 value = "nseqs",
                                 fill = 0 # replace the NA by 0
)

epathogen.subset_relseqs <- spread(epathogen.subset[, !(names(epathogen.subset) %in% c("nseqs",
                                                                              "relseqs_rar",
                                                                              "nseqs_rar"))], 
                                   key ="sample_id", 
                                   value = "relseqs",
                                   fill = 0 # replace the NA by 0
)

epathogen.subset_nseqs_rar <- spread(epathogen.subset[, !(names(epathogen.subset) %in% c("relseqs_rar", 
                                                                                "relseqs", 
                                                                                "nseqs"))], 
                                     key ="sample_id", 
                                     value = "nseqs_rar",
                                     fill = 0 # replace the NA by 0
)

epathogen.subset_relseqs_rar <- spread(epathogen.subset[, !(names(epathogen.subset) %in% c("nseqs_rar", 
                                                                                  "relseqs", 
                                                                                  "nseqs"))], 
                                       key ="sample_id", 
                                       value = "relseqs_rar",
                                       fill = 0 # replace the NA by 0

)


## write output files
write.table(epathogen.subset_nseqs, file = "results/epathogen_16S_nseqs.txt", sep = ";")
write.table(epathogen.subset_relseqs, file = "results/epathogen_16S_relseqs.txt", sep = ";")

write.table(epathogen.subset_nseqs_rar, file = "results/epathogen_16S_nseqs_rar.txt", sep = ";")
write.table(epathogen.subset_relseqs_rar, file = "results/epathogen_16S_relseqs_rar.txt", sep = ";")


