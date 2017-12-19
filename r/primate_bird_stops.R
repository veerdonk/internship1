# install.packages("ggnetwork")
# install.packages("ggplot2")
library("ggplot2")
library("ggnetwork")

setwd("~/Documents/david/scripts_git/out/stops")
#open file
# alltypes <- read.csv("./human_chicken.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
alltypes <- read.csv("./human_chicken_one2one.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
secondCodon <- read.csv("./human_chicken_second_codon.tsv", sep="\t", header = FALSE, stringsAsFactors = FALSE)
colnames(alltypes) <- c("humanID", "chickenID", "humanStop", "chickenStop")

#allowed stop codons
stopCodons <- as.factor(c("TAG", "TAA", "TGA"))

#filter on codons in stopCodons
alltypesGoodStops <- alltypes[alltypes$humanStop %in% stopCodons & alltypes$chickenStop %in% stopCodons,]

#conserved stop codons between human & chicken (all)
sum(alltypesGoodStops$humanStop == alltypesGoodStops$chickenStop)

#human TGA == chicken TGA
tgaConservation <- sum(alltypesGoodStops$humanStop == "TGA" & alltypesGoodStops$chickenStop == "TGA")/sum(alltypesGoodStops$humanStop == "TGA")
#human TAA == chicken TAA
taaConservation <- sum(alltypesGoodStops$humanStop == "TAA" & alltypesGoodStops$chickenStop == "TAA")/sum(alltypesGoodStops$humanStop == "TAA")
#human TAG == chicken TAG
tagConservation <- sum(alltypesGoodStops$humanStop == "TAG" & alltypesGoodStops$chickenStop == "TAG")/sum(alltypesGoodStops$humanStop == "TAG")

stopcodonTrans <- data.frame(TAG = numeric(0), TAA = numeric(0), TGA = numeric(0))

codonFreq <- getCodonFreq(alltypesGoodStops)

for (stop in stopCodons){
  stopToTga <- round(sum(alltypesGoodStops$humanStop == stop & alltypesGoodStops$chickenStop == "TGA")/sum(alltypesGoodStops$humanStop == stop), 5)
  stopToTag <- round(sum(alltypesGoodStops$humanStop == stop & alltypesGoodStops$chickenStop == "TAG")/sum(alltypesGoodStops$humanStop == stop), 5)
  stopToTaa <- round(sum(alltypesGoodStops$humanStop == stop & alltypesGoodStops$chickenStop == "TAA")/sum(alltypesGoodStops$humanStop == stop), 5)
  print(stop)
  print(stopToTaa)
  stopcodonTrans[stop,] <- c(stopToTag, stopToTaa, stopToTga)
}

getCodonFreq <- function(x){
  stopcodonTest <- data.frame(TAG = numeric(0), TAA = numeric(0), TGA = numeric(0))
  for (stop in stopCodons){
    stopToTga <- round(sum(x$humanStop == stop & x$chickenStop == "TGA")/sum(x$humanStop == stop), 2)
    stopToTag <- round(sum(x$humanStop == stop & x$chickenStop == "TAG")/sum(x$humanStop == stop), 2)
    stopToTaa <- round(sum(x$humanStop == stop & x$chickenStop == "TAA")/sum(x$humanStop == stop), 2)
    stopcodonTest[stop,] <- c(stopToTag, stopToTaa, stopToTga)
  }
  return(stopcodonTest)
}
getCodonFreq(alltypesGoodStops)

as.array(as.matrix(stopcodonTrans))


#aa = R
aaCodons <- as.factor(c("CGT",
                        "CGC",
                        "CGA",
                        "CGG"))
#filter on amino acid
secondCodonFiltered <- secondCodon[secondCodon$V3 %in% aaCodons & secondCodon$V4 %in% aaCodons,]

sum(secondCodonFiltered$V3 == secondCodonFiltered$V4)

cod1 <- "CGT"
cod2 <- "CGC"
cod3 <- "CGA"
cod4 <- "CGG"


#codon 1 > codon 1 preservation
sum(secondCodonFiltered$V3 == cod1 & secondCodonFiltered$V4 == cod1)/sum(secondCodonFiltered$V3 == cod1)*100
#codon 2 > codon 2 preservation
sum(secondCodonFiltered$V3 == cod2 & secondCodonFiltered$V4 == cod2)/sum(secondCodonFiltered$V3 == cod2)*100
#codon 3 > codon 3 preservation
sum(secondCodonFiltered$V3 == cod3 & secondCodonFiltered$V4 == cod3)/sum(secondCodonFiltered$V3 == cod3)*100
#codon 4 > codon 4 preservation
sum(secondCodonFiltered$V3 == cod4 & secondCodonFiltered$V4 == cod4)/sum(secondCodonFiltered$V3 == cod4)*100

#expression data - top 100 expressed in whole blood
expr_wholeBlood <- read.csv("~/data/expression/top100expressedWholeBlood.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
expr_brainCortex <- read.csv("~/data/expression/top100expressed_brain.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)


#fix IDs
expr_wholeBlood$Gencode.Id <- apply(expr_wholeBlood, 1, function(x) strsplit(x[1], "\\.")[[1]][1])
expr_brainCortex$Gencode.Id <- apply(expr_brainCortex, 1, function(x) strsplit(x[1], "\\.")[[1]][1])
eQTL_blood$Gencode.Id <- apply(eQTL_blood, 1, function(x) strsplit(x[1], "\\.")[[1]][1]) 

#filter stopcodon table on genes in expression table
highlyExpressedBlood <- alltypesGoodStops[alltypesGoodStops$humanID %in% expr_wholeBlood$Gencode.Id,]
highlyExpressedBrain <- alltypesGoodStops[alltypesGoodStops$humanID %in% expr_brainCortex$Gencode.Id,]
highlyExpressedBloodQTL <- alltypesGoodStops[alltypesGoodStops$humanID %in% eQTL_blood$Gencode.Id,]

#frequencies in all genes
freq<- c(sum(alltypesGoodStops$humanStop == "TAA")/length(alltypesGoodStops$humanID)*100,
sum(alltypesGoodStops$humanStop == "TGA")/length(alltypesGoodStops$humanID)*100,
sum(alltypesGoodStops$humanStop == "TAG")/length(alltypesGoodStops$humanID)*100)

for (stop in stopCodons){
  print(stop)
  print(sum(highlyExpressedBrain$humanStop == stop & highlyExpressedBrain$chickenStop == stop))
  print(sum(highlyExpressedBlood$humanStop == stop & highlyExpressedBlood$chickenStop == stop))
}
## load expression tables ##
expressionFiles <- dir("~/data/expression/", full.names = TRUE)
expressionNames <- c("adipose", "adrenalGland", "artery", "blood", "brainCortex", "brainHippocampus", "colon", "esophagus", "heart", "liver", "lung", "nerve", "ovary", "pancreas", "prostate", "skeletalMuscle", "smallIntestine", "spleen", "testis", "thyroid", "uterus", "vagina")
expressionTables <- c()

stopcodonConservationHighlyExpressed <- data.frame(TAG = numeric(0), TAA = numeric(0), TGA = numeric(0))
totals <- data.frame(orthologs = numeric(0), TAG = numeric(0), TAA = numeric(0), TGA = numeric(0), TAGConserved = numeric(0), TAAConserved = numeric(0), TGAConserved = numeric(0))


stopcodonTrans <- data.frame(TAG = numeric(0), TAA = numeric(0), TGA = numeric(0))
stopcodonTrans["TAG",] <- c(0,0,0)
stopcodonTrans["TAA",] <- c(0,0,0)
stopcodonTrans["TGA",] <- c(0,0,0)



for (file in expressionFiles){
  print(file)
  table <- read.csv(file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  table$Gencode.Id <- apply(table, 1, function(x) strsplit(x[1], "\\.")[[1]][1])
  highlyExpressed <- alltypesGoodStops[alltypesGoodStops$humanID %in% table$Gencode.Id,]
  
  
  
  stopcodonTrans <- stopcodonTrans + getCodonFreq(highlyExpressed)
  totals[expressionNames[which(expressionFiles == file)],] <-  c(
    length(highlyExpressed$humanStop),
    sum(highlyExpressed$humanStop == "TAG"),
    sum(highlyExpressed$humanStop == "TAA"),
    sum(highlyExpressed$humanStop == "TGA"),
    sum(highlyExpressed$humanStop == "TAG" & highlyExpressed$chickenStop == "TAG"),
    sum(highlyExpressed$humanStop == "TAA" & highlyExpressed$chickenStop == "TAA"),
    sum(highlyExpressed$humanStop == "TGA" & highlyExpressed$chickenStop == "TGA")
  )
  stopcodonConservationHighlyExpressed[expressionNames[which(expressionFiles == file)],] <- c(
    round(sum(highlyExpressed$humanStop == "TAG" & highlyExpressed$chickenStop == "TAG")/sum(highlyExpressed$humanStop == "TAG"),2), 
    round(sum(highlyExpressed$humanStop == "TAA" & highlyExpressed$chickenStop == "TAA")/sum(highlyExpressed$humanStop == "TAA"),2), 
    round(sum(highlyExpressed$humanStop == "TGA" & highlyExpressed$chickenStop == "TGA")/sum(highlyExpressed$humanStop == "TGA"),2))
}

conserv <- stopcodonConservationHighlyExpressed

round(stopcodonTrans/22, 2)

median(stopcodonConservationHighlyExpressed$TAG)
median(stopcodonConservationHighlyExpressed$TAA)
median(stopcodonConservationHighlyExpressed$TGA)

stopcodonConservationHighlyExpressed[stopcodonConservationHighlyExpressed$TGA > tgaConservation,]
stopcodonConservationHighlyExpressed[stopcodonConservationHighlyExpressed$TAA < taaConservation,]
stopcodonConservationHighlyExpressed[stopcodonConservationHighlyExpressed$TAG > tagConservation,]
# pdf("stopcodonConservation.pdf")

boxplot(totals[2:4], main = "top 100 genes of 22 tissues (counts)", ylab = "number of genes", xlab = "stopcodon")

boxplot(stopcodonConservationHighlyExpressed, main = "stop conservation in highly expressed", ylab = "conservation", xlab = "stop codon")
x0s <- 1:3 - 0.4
x1s <- 1:3 + 0.4
y0s <- c(tagConservation, taaConservation, tgaConservation)
segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red")

# dev.off()