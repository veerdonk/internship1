setwd("~/Documents/david/scripts_git/out/stops")




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

## Expression files and their names
expressionFiles <- dir("~/data/expression/", full.names = TRUE)
expressionNames <- c("adipose", "adrenalGland", "artery", "blood", "brainCortex", "brainHippocampus", "colon", "esophagus", "heart", "liver", "lung", "nerve", "ovary", "pancreas", "prostate", "skeletalMuscle", "smallIntestine", "spleen", "testis", "thyroid", "uterus", "vagina")

## Legal stopcodons
stopCodons <- as.factor(c("TAG", "TAA", "TGA"))
orthologsFiltered <- one2oneOrthologs[one2oneOrthologs$humanStop %in% stopCodons & one2oneOrthologs$orthologStop %in% stopCodons,]

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
# highlyExpressedStops <- round(stopcodonTrans/length(stopcodonConservationHighlyExpressed$TAG),2)


human_garList     <- calculateHighlyExpressedStops(one2oneOrthologs_human_gar)
human_chickenList <- calculateHighlyExpressedStops(one2oneOrthologs_human_chicken)
human_dolphinList <- calculateHighlyExpressedStops(one2oneOrthologs_human_dolphin)
human_mouseList   <- calculateHighlyExpressedStops(one2oneOrthologs_human_mouse)




as.data.frame(human_mouseList[4]) - as.data.frame(human_mouseList[2])
as.data.frame(human_dolphinList[2]) - as.data.frame(human_dolphinList[4])
as.data.frame(human_chickenList[2]) - as.data.frame(human_chickenList[4])
as.data.frame(human_garList[2]) - as.data.frame(human_garList[4])
