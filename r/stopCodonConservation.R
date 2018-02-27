setwd("~/Documents/david/scripts_git/out/stops")

## Long lived 
longLivedStops <- read.delim("./allLongLivedWithStops.tsv", header = FALSE, stringsAsFactors = FALSE)
longLivedStops_filtered <- longLivedStops[longLivedStops$V2 %in% c("TAA", "TGA", "TAG"),]
sum(longLivedStops_filtered$V2 == "TAA")/length(longLivedStops_filtered$V2)
sum(longLivedStops_filtered$V2 == "TGA")/length(longLivedStops_filtered$V2)
sum(longLivedStops_filtered$V2 == "TAG")/length(longLivedStops_filtered$V2)

## Ortholog file and column names
coln <- c("humanID", "orthologID", "humanStop", "orthologStop")
one2oneOrthologs_human_chicken  <- read.csv("./human_chicken_one2one.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
one2oneOrthologs_human_mouse    <- read.csv("./human_mouse_one2one.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
one2oneOrthologs_human_dolphin  <- read.csv("./human_dolphin_one2one.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
one2oneOrthologs_human_gar      <- read.csv("./human_gar_one2one.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

colnames(one2oneOrthologs_human_chicken)  <- coln
colnames(one2oneOrthologs_human_mouse)    <- coln
colnames(one2oneOrthologs_human_dolphin)  <- coln
colnames(one2oneOrthologs_human_gar)      <- coln

## stop + dN/dS files
#primates
csa_dnds_stop <- read.csv("./primates/csa.csv", header = FALSE, stringsAsFactors = FALSE)
ggo_dnds_stop <- read.csv("./primates/ggo.csv", header = FALSE, stringsAsFactors = FALSE)
mmu_dnds_stop <- read.csv("./primates/mmu.csv", header = FALSE, stringsAsFactors = FALSE)
pan_dnds_stop <- read.csv("./primates/pan.csv", header = FALSE, stringsAsFactors = FALSE)
soe_dnds_stop <- read.csv("./primates/soe.csv", header = FALSE, stringsAsFactors = FALSE)
tsy_dnds_stop <- read.csv("./primates/tsy.csv", header = FALSE, stringsAsFactors = FALSE)
primate_dnds_files <- list(csa_dnds_stop, ggo_dnds_stop, mmu_dnds_stop, pan_dnds_stop, soe_dnds_stop, tsy_dnds_stop)

#rodents
rodents_dir <- dir("./rodents/", full.names = TRUE)
for(filename in rodents_dir){
  rodent <- read.csv(filename, header = FALSE, stringsAsFactors = FALSE)
  filterdnds(rodent)
}


## Filter/plot dnds files
filterdnds <- function(x){
  filteredDnds <- x[x$V5 > 0 & x$V6 > 0.01 & x$V5 < 2 & x$V6 < 2 & x$V4 < 10 & x$V4 < 1 & x$V7 %in% stops,]
  
  # allTGA = c(allTGA, filteredDnds$V4[filteredDnds$V7 == "TGA"])
  # allTAA = c(allTAA, filteredDnds$V4[filteredDnds$V7 == "TAA"])
  # allTAG = c(allTAG, filteredDnds$V4[filteredDnds$V7 == "TAG"])
  boxplot(log(filteredDnds$V4) ~ filteredDnds$V7, data = filteredDnds, main = titles[which(files == file)])
}

for(primate in primate_dnds_files){
  filterdnds(as.data.frame(primate))
}


## Expression files and their names
expressionFiles <- dir("~/data/expression/", full.names = TRUE)
expressionNames <- c("adipose", "adrenalGland", "artery", "blood", "brainCortex", "brainHippocampus", "colon", "esophagus", "heart", "liver", "lung", "nerve", "ovary", "pancreas", "prostate", "skeletalMuscle", "smallIntestine", "spleen", "testis", "thyroid", "uterus", "vagina")

## Legal stopcodons
stopCodons <- as.factor(c("TAG", "TAA", "TGA"))
one2oneOrthologs_human_chicken <- one2oneOrthologs_human_chicken[one2oneOrthologs_human_chicken$humanStop %in% stopCodons & one2oneOrthologs_human_chicken$orthologStop %in% stopCodons,]
one2oneOrthologs_human_dolphin <- one2oneOrthologs_human_dolphin[one2oneOrthologs_human_dolphin$humanStop %in% stopCodons & one2oneOrthologs_human_dolphin$orthologStop %in% stopCodons,]
one2oneOrthologs_human_gar     <- one2oneOrthologs_human_gar[one2oneOrthologs_human_gar$humanStop %in% stopCodons & one2oneOrthologs_human_gar$orthologStop %in% stopCodons,]
one2oneOrthologs_human_mouse   <- one2oneOrthologs_human_mouse[one2oneOrthologs_human_mouse$humanStop %in% stopCodons & one2oneOrthologs_human_mouse$orthologStop %in% stopCodons,]


## Stop codon distribution
getCodonFreqOneSpecies <- function(x){
  TAG <- length(x[x$orthologStop == "TAG",1])/length(x$orthologStop)
  TGA <- length(x[x$orthologStop == "TGA",1])/length(x$orthologStop)
  TAA <- length(x[x$orthologStop == "TAA",1])/length(x$orthologStop)
  return(c(TAG, TAA, TGA))
}


#human TGA == ortholog TGA ratio
tgaConservation <- sum(orthologsFiltered$humanStop == "TGA" & orthologsFiltered$orthologStop == "TGA")/sum(orthologsFiltered$humanStop == "TGA")
#human TAA == ortholog TAA ratio
taaConservation <- sum(orthologsFiltered$humanStop == "TAA" & orthologsFiltered$orthologStop == "TAA")/sum(orthologsFiltered$humanStop == "TAA")
#human TAG == ortholog TAG ratio
tagConservation <- sum(orthologsFiltered$humanStop == "TAG" & orthologsFiltered$orthologStop == "TAG")/sum(orthologsFiltered$humanStop == "TAG")

## Function to calculate stopcodon frequency
getCodonFreq <- function(x){
  stopcodonFreq <- data.frame(TAG = numeric(0), TAA = numeric(0), TGA = numeric(0))
  for (stop in stopCodons){
    stopToTga <- round(sum(x$humanStop == stop & x$orthologStop == "TGA")/sum(x$humanStop == stop), 2)
    stopToTag <- round(sum(x$humanStop == stop & x$orthologStop == "TAG")/sum(x$humanStop == stop), 2)
    stopToTaa <- round(sum(x$humanStop == stop & x$orthologStop == "TAA")/sum(x$humanStop == stop), 2)
    stopcodonFreq[stop,] <- c(stopToTag, stopToTaa, stopToTga)
  }
  return(stopcodonFreq)
}

# Stopcodon frequency in all ortholog genes
stopFrequency <- getCodonFreq(orthologsFiltered)


calculateHighlyExpressedStops <- function(x){
  stopFrequency <- getCodonFreq(x)
  ## Create placeholder tables to be filled
  stopcodonTrans <- data.frame(
    TAG = numeric(0), 
    TAA = numeric(0), 
    TGA = numeric(0))
  stopcodonTrans["TAG",] <- c(0,0,0)
  stopcodonTrans["TAA",] <- c(0,0,0)
  stopcodonTrans["TGA",] <- c(0,0,0)
  
  stopcodonConservationHighlyExpressed <- data.frame(
    TAG = numeric(0),
    TAA = numeric(0),
    TGA = numeric(0))
  
  totals <- data.frame(
    orthologs = numeric(0),
    TAG = numeric(0),
    TAA = numeric(0),
    TGA = numeric(0),
    TAGConserved = numeric(0),
    TAAConserved = numeric(0),
    TGAConserved = numeric(0))
  
  
  
  for (file in expressionFiles){
    print(file)
    expressionTable <- read.csv(file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
    expressionTable$Gencode.Id <- apply(expressionTable, 1, function(x) strsplit(x[1], "\\.")[[1]][1])
    highlyExpressed <- x[x$humanID %in% expressionTable$Gencode.Id,]
    
    stopcodonTrans <- stopcodonTrans + getCodonFreq(highlyExpressed)
    print(getCodonFreq(highlyExpressed))
    totals[expressionNames[which(expressionFiles == file)],] <-  c(
      length(highlyExpressed$humanStop),
      sum(highlyExpressed$humanStop == "TAG"),
      sum(highlyExpressed$humanStop == "TAA"),
      sum(highlyExpressed$humanStop == "TGA"),
      sum(highlyExpressed$humanStop == "TAG" & highlyExpressed$orthologStop == "TAG"),
      sum(highlyExpressed$humanStop == "TAA" & highlyExpressed$orthologStop == "TAA"),
      sum(highlyExpressed$humanStop == "TGA" & highlyExpressed$orthologStop == "TGA")
    )
    stopcodonConservationHighlyExpressed[expressionNames[which(expressionFiles == file)],] <- c(
      round(sum(highlyExpressed$humanStop == "TAG" & highlyExpressed$orthologStop == "TAG")/sum(highlyExpressed$humanStop == "TAG"),2), 
      round(sum(highlyExpressed$humanStop == "TAA" & highlyExpressed$orthologStop == "TAA")/sum(highlyExpressed$humanStop == "TAA"),2), 
      round(sum(highlyExpressed$humanStop == "TGA" & highlyExpressed$orthologStop == "TGA")/sum(highlyExpressed$humanStop == "TGA"),2))
  }
  return(list(stopcodonConservationHighlyExpressed, round(stopcodonTrans/length(stopcodonConservationHighlyExpressed$TAG),2), totals, stopFrequency))
}
highlyExpressedStops <- round(stopcodonTrans/length(stopcodonConservationHighlyExpressed$TAG),2)


human_garList     <- calculateHighlyExpressedStops(one2oneOrthologs_human_gar)
human_chickenList <- calculateHighlyExpressedStops(one2oneOrthologs_human_chicken)
human_dolphinList <- calculateHighlyExpressedStops(one2oneOrthologs_human_dolphin)
human_mouseList   <- calculateHighlyExpressedStops(one2oneOrthologs_human_mouse)

as.data.frame(human_mouseList[2]) - as.data.frame(human_mouseList[4])
as.data.frame(human_dolphinList[2]) - as.data.frame(human_dolphinList[4])
as.data.frame(human_chickenList[2]) - as.data.frame(human_chickenList[4])
as.data.frame(human_garList[2]) - as.data.frame(human_garList[4])

human_gar_taa     <- one2oneOrthologs_human_gar[one2oneOrthologs_human_gar$humanStop == "TAA" & one2oneOrthologs_human_gar$orthologStop == "TAA", 1]
human_chicken_taa <- one2oneOrthologs_human_chicken[one2oneOrthologs_human_chicken$humanStop == "TAA" & one2oneOrthologs_human_chicken$orthologStop == "TAA", 1]
human_mouse_taa   <- one2oneOrthologs_human_mouse[one2oneOrthologs_human_mouse$humanStop == "TAA" & one2oneOrthologs_human_mouse$orthologStop == "TAA", 1]
human_dolphin_taa <- one2oneOrthologs_human_dolphin[one2oneOrthologs_human_dolphin$humanStop == "TAA" & one2oneOrthologs_human_dolphin$orthologStop == "TAA", 1]

geneUniverseTaa <- unique(c(human_chicken_taa, human_gar_taa, human_dolphin_taa, human_mouse_taa))
geneUniverseAll <- unique(c(one2oneOrthologs_human_chicken$humanID, one2oneOrthologs_human_dolphin$humanID, one2oneOrthologs_human_gar$humanID, one2oneOrthologs_human_mouse$humanID))

calculateIntersectForStop <- function(stop){
  intersection = Reduce(intersect, list(
  one2oneOrthologs_human_chicken[one2oneOrthologs_human_chicken$humanStop == stop & one2oneOrthologs_human_chicken$orthologStop == stop, 1],
  one2oneOrthologs_human_dolphin[one2oneOrthologs_human_dolphin$humanStop == stop & one2oneOrthologs_human_dolphin$orthologStop == stop, 1],
  one2oneOrthologs_human_gar[one2oneOrthologs_human_gar$humanStop == stop & one2oneOrthologs_human_gar$orthologStop == stop, 1],
  one2oneOrthologs_human_mouse[one2oneOrthologs_human_mouse$humanStop == stop & one2oneOrthologs_human_mouse$orthologStop == stop, 1]
  ))
  return(intersection)
}

nonConserved = Reduce(intersect, list(
  one2oneOrthologs_human_chicken[one2oneOrthologs_human_chicken$humanStop != one2oneOrthologs_human_chicken$orthologStop,1],
  one2oneOrthologs_human_dolphin[one2oneOrthologs_human_dolphin$humanStop != one2oneOrthologs_human_dolphin$orthologStop,1],
  one2oneOrthologs_human_gar[one2oneOrthologs_human_gar$humanStop != one2oneOrthologs_human_gar$orthologStop,1],
  one2oneOrthologs_human_mouse[one2oneOrthologs_human_mouse$humanStop != one2oneOrthologs_human_mouse$orthologStop,1]
))
# write.csv(nonConserved, file = "~/data/stopCodons/goTerms/nonConserved.csv", row.names = F)

conservedTga <- calculateIntersectForStop("TGA")
# write.csv(conservedTga, file = "~/data/stopCodons/goTerms/conservedTga.csv", row.names = F)
conservedTag <- calculateIntersectForStop("TAG")
# write.csv(conservedTag, file = "~/data/stopCodons/goTerms/conservedTag.csv", row.names = F)

# sum(geneUniverseAll %in% conserved_taa)

# write.csv(geneUniverseAll, file = "~/data/stopCodons/goTerms/allGenes.csv", row.names = F)

test1 <- c("aaa", "bbb", "ccc", "ddd")
test2 <- c("bbb", "ddd", "eee", "fff")

conserved_taa <- Reduce(intersect, list(human_gar_taa, human_chicken_taa, human_mouse_taa, human_dolphin_taa))
# write.csv(conserved_taa, file = "~/data/stopCodons/goTerms/conserved_taa.csv", row.names = F)
for (file in expressionFiles){
  print(file)
  expressionTable <- read.csv(file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  expressionTable$Gencode.Id <- apply(expressionTable, 1, function(x) strsplit(x[1], "\\.")[[1]][1])
  print(intersect(expressionTable$Gencode.Id, conserved_taa))
}
