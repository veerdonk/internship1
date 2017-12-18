# install.packages("ggnetwork")
# install.packages("ggplot2")
library("ggplot2")
library("ggnetwork")

setwd("~/Documents/david/scripts_git/out/stops")
#open file
alltypes <- read.csv("./human_chicken.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
secondCodon <- read.csv("./human_chicken_second_codon.tsv", sep="\t", header = FALSE, stringsAsFactors = FALSE)
colnames(alltypes) <- c("humanID", "chickenID", "humanStop", "chickenStop")

#allowed stop codons
stopCodons <- as.factor(c("TAG", "TAA", "TGA"))

#filter on codons in stopCodons
alltypesGoodStops <- alltypes[alltypes$humanStop %in% stopCodons & alltypes$chickenStop %in% stopCodons,]

#conserved stop codons between human & chicken (all)
sum(alltypesGoodStops$humanStop == alltypesGoodStops$chickenStop)

#human TGA == chicken TGA
sum(alltypesGoodStops$humanStop == "TGA" & alltypesGoodStops$chickenStop == "TGA")/sum(alltypesGoodStops$humanStop == "TGA")*100
#human TAA == chicken TAA
sum(alltypesGoodStops$humanStop == "TAA" & alltypesGoodStops$chickenStop == "TAA")/sum(alltypesGoodStops$humanStop == "TAA")*100
#human TAG == chicken TAG
sum(alltypesGoodStops$humanStop == "TAG" & alltypesGoodStops$chickenStop == "TAG")/sum(alltypesGoodStops$humanStop == "TAG")*100

stopcodonTrans <- data.frame(TAG = numeric(0), TAA = numeric(0), TGA = numeric(0))

for (stop in stopCodons){
  stopToTga <- sum(alltypesGoodStops$humanStop == stop & alltypesGoodStops$chickenStop == "TGA")
  stopToTag <- sum(alltypesGoodStops$humanStop == stop & alltypesGoodStops$chickenStop == "TAG")
  stopToTaa <- sum(alltypesGoodStops$humanStop == stop & alltypesGoodStops$chickenStop == "TAA")
  stopcodonTrans[stop,] <- c(stopToTag, stopToTaa, stopToTga)
}

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

#frequencies in higly expressed genes
bloodFreq <- c(sum(highlyExpressedBlood$humanStop == "TAA")/length(highlyExpressedBlood$humanID)*100,
sum(highlyExpressedBlood$humanStop == "TGA")/length(highlyExpressedBlood$humanID)*100,
sum(highlyExpressedBlood$humanStop == "TAG")/length(highlyExpressedBlood$humanID)*100)

sum(highlyExpressedBlood$chickenStop == "TAA")/length(highlyExpressedBlood$chickenID)*100
sum(highlyExpressedBlood$chickenStop == "TGA")/length(highlyExpressedBlood$chickenID)*100
sum(highlyExpressedBlood$chickenStop == "TAG")/length(highlyExpressedBlood$chickenID)*100

# conserved stops in highly expressed genes
#brain
taabrain <- sum(highlyExpressedBrain$humanStop == "TAA" & highlyExpressedBrain$chickenStop == "TAA")/sum(highlyExpressedBrain$humanStop == "TAA")*100
tagbrain <- sum(highlyExpressedBrain$humanStop == "TAG" & highlyExpressedBrain$chickenStop == "TAG")/sum(highlyExpressedBrain$humanStop == "TAG")*100
tgabrain <- sum(highlyExpressedBrain$humanStop == "TGA" & highlyExpressedBrain$chickenStop == "TGA")/sum(highlyExpressedBrain$humanStop == "TGA")*100

#blood
taablood <- sum(highlyExpressedBlood$humanStop == "TAA" & highlyExpressedBlood$chickenStop == "TAA")/sum(highlyExpressedBlood$humanStop == "TAA")*100
tagblood <- sum(highlyExpressedBlood$humanStop == "TAG" & highlyExpressedBlood$chickenStop == "TAG")/sum(highlyExpressedBlood$humanStop == "TAG")*100
tgablood <- sum(highlyExpressedBlood$humanStop == "TGA" & highlyExpressedBlood$chickenStop == "TGA")/sum(highlyExpressedBlood$humanStop == "TGA")*100

#blood eQTL
sum(highlyExpressedBloodQTL$humanStop == "TAA" & highlyExpressedBloodQTL$chickenStop == "TAA")/sum(highlyExpressedBloodQTL$humanStop == "TAA")*100

for (stop in stopCodons){
  print(stop)
  print(sum(highlyExpressedBrain$humanStop == stop & highlyExpressedBrain$chickenStop == stop))
  print(sum(highlyExpressedBlood$humanStop == stop & highlyExpressedBlood$chickenStop == stop))
}
## load expression tables ##
expressionFiles <- dir("~/data/expression/", full.names = TRUE)
expressionNames <- c("adipose", "artery", "blood", "brain", "colon", "heart", "liver", "lung", "spleen", "testis", "thyroid", "uterus")
expressionTables <- c()

stopcodonConservationHighlyExpressed <- data.frame(TAG = numeric(0), TAA = numeric(0), TGA = numeric(0))

for (file in expressionFiles){
  print(file)
  table <- read.csv(file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  table$Gencode.Id <- apply(table, 1, function(x) strsplit(x[1], "\\.")[[1]][1])
  highlyExpressed <- alltypesGoodStops[alltypesGoodStops$humanID %in% table$Gencode.Id,]
  
  stopcodonConservationHighlyExpressed[expressionNames[which(expressionFiles == file)],] <- c(
    sum(highlyExpressed$humanStop == "TAG" & highlyExpressed$chickenStop == "TAG")/sum(highlyExpressed$humanStop == "TAG")*100, 
    sum(highlyExpressed$humanStop == "TAA" & highlyExpressed$chickenStop == "TAA")/sum(highlyExpressed$humanStop == "TAA")*100, 
    sum(highlyExpressed$humanStop == "TGA" & highlyExpressed$chickenStop == "TGA")/sum(highlyExpressed$humanStop == "TGA")*100)
}

expressionTables[1]
