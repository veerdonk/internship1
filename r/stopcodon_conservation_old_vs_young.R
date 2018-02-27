long_lived <- read.csv("~/data/orthologs/old/conservedStopsOld.csv")
short_lived <- read.csv("~/data/orthologs/young/conservedStopYoungAll.csv")
long_lived$longevity <- "long"
short_lived$longevity <- "short"
merged_non_t <- merge(long_lived, short_lived, all = TRUE)

merged_non_t[c(5,6), -1]

library(reshape2)
library(ggplot2)

number_stops_melt <- melt(merged_non_t[c(5,6), -1],id.vars="longevity")
ggplot(number_stops_melt, aes(x=variable, y=value, fill=factor(longevity))) +
  geom_bar(stat="identity", position = "dodge")+
  scale_fill_discrete("Longevity")+
  xlab("stop codon")+ylab("number of conserved stops")

conservedFraction <- melt(merged_non_t[c(1,2),-1], id.vars = "longevity")
ggplot(conservedFraction, aes(x=variable, y=value, fill=factor(longevity)))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_discrete("Longevity")+
  xlab("stop codon")+ylab("percentage of total conserved stops")

conservedOfAll <- melt(merged_non_t[c(3,4), -1], id.vars = "longevity")
ggplot(conservedOfAll, aes(x=variable, y=value, fill=factor(longevity)))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_discrete("Longevity")+
  xlab("stop codon")+ylab("percentage of all stops")

#############################################################################################################
## significant/longlived/shortlived
expressionFiles <- dir("~/data/expression/", full.names = TRUE)
expressionNames <- c("adipose", "adrenalGland", "artery", "blood", "brainCortex", "brainHippocampus", "colon", "esophagus", "heart", "liver", "lung", "nerve", "ovary", "pancreas", "prostate", "skeletalMuscle", "smallIntestine", "spleen", "testis", "thyroid", "uterus", "vagina")

significantStops <- read.delim("~/data/stopCodons/significantStops.tsv", header = FALSE)
colnames(significantStops) <- c("geneid", "stop")

stopsOld <- read.delim("~/data/orthologs/old/idToConservedStop.csv", header = FALSE)
colnames(stopsOld) <- c("geneid", "stop")

stopsYoung <- read.delim("~/data/orthologs/young/idToConservedStop.csv", header = FALSE)
colnames(stopsYoung) <- c("geneid", "stop")

idMapping <- read.delim("~/data/orthologs/primates/hsa_csa_mapping.txt")

stopcodonConservationHighlyExpressedOld <- data.frame(TAG = numeric(0), TAA = numeric(0), TGA = numeric(0))
stopcodonConservationHighlyExpressedYoung <- data.frame(TAG = numeric(0), TAA = numeric(0), TGA = numeric(0))





for (file in expressionFiles){
  print(file)
  expressionTable <- read.csv(file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  expressionTable$Gencode.Id <- apply(expressionTable, 1, function(x) strsplit(x[1], "\\.")[[1]][1])
  
  highlyExpressedOld <- stopsOld[stopsOld$geneid %in% expressionTable$Gencode.Id,]
  highlyExpressedYoung <- stopsYoung[idMapping[stopsYoung$geneid,1] %in% expressionTable$Gencode.Id,]
  
  stopcodonConservationHighlyExpressedYoung[expressionNames[which(expressionFiles ==file)],] <- c(
    sum(highlyExpressedYoung$stop == "TAG")/length(highlyExpressedYoung$geneid),
    sum(highlyExpressedYoung$stop == "TAA")/length(highlyExpressedYoung$geneid),
    sum(highlyExpressedYoung$stop == "TGA")/length(highlyExpressedYoung$geneid)
  )
  
  
  stopcodonConservationHighlyExpressedOld[expressionNames[which(expressionFiles == file)],] <- c(
    sum(highlyExpressedOld$stop == "TAG")/length(highlyExpressedOld$geneid),
    sum(highlyExpressedOld$stop == "TAA")/length(highlyExpressedOld$geneid),
    sum(highlyExpressedOld$stop == "TGA")/length(highlyExpressedOld$geneid)
  )
  
}
oldTaa <- sum(stopsOld$stop == "TAA")/length(stopsOld$geneid)
oldTga <- sum(stopsOld$stop == "TGA")/length(stopsOld$geneid)
oldTag <- sum(stopsOld$stop == "TAG")/length(stopsOld$geneid)

youngTaa <- sum(stopsYoung$stop == "TAA")/length(stopsOld$geneid)
youngTga <- sum(stopsYoung$stop == "TGA")/length(stopsYoung$geneid)
youngTag <- sum(stopsYoung$stop == "TAG")/length(stopsYoung$geneid)

boxplot(stopcodonConservationHighlyExpressedOld, main = "stop conservation in highly expressed (long-lived)", ylab = "fraction", xlab = "stop codon")
x0s <- 1:3 - 0.4
x1s <- 1:3 + 0.4
y0s <- c(oldTag, oldTaa, oldTga)
segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red")

boxplot(stopcodonConservationHighlyExpressedYoung, main = "stop conservation in highly expressed (short-lived)", ylab = "fraction", xlab = "stop codon")
y1s <- c(youngTag, youngTaa, youngTga)
segments(x0 = x0s, x1 = x1s, y0 = y1s, col = "red")

for (file in expressionFiles){
  print(file)
  expressionTable <- read.csv(file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  expressionTable$Gencode.Id <- apply(expressionTable, 1, function(x) strsplit(x[1], "\\.")[[1]][1])
  print(intersect(expressionTable$Gencode.Id, significantStops$geneid))
}

#######################################################################################################################
## stops for significant genes
short <- c("oro", "bacu", "tut", "ehe", "mgp", "cca", "pad", "fab", "aan", "tsy", "soe", "csa", "cap", "sto", "mau", "rno", "mus", "mpa", "ola", "xma", "amx", "pfo", "nfu", "min")

long <- c("phm", "pcr", "tal", "hle", "gga", "afo", "apla", "ggo", "pan", "mmu", "can", "fda", "cla", "ode", "dor", "pem", "icu", "omy", "oni", "myd", "pale", "efu" )

sigGenes <- read.csv("/home/dylan/Documents/david/scripts_git/internship1/r/model/results/B_All Data.csv")
colnames(sigGenes) <- c("gene", "pval", "adj_pval")

stopsOfSigGenes <- read.csv("~/Documents/david/scripts_git/out/significant/significantstops.csv", header = FALSE, stringsAsFactors = FALSE)
stopsOfSigGenes[1,1] <- "organism"
colnames(stopsOfSigGenes) <- stopsOfSigGenes[1,]
stopsOfSigGenes <- stopsOfSigGenes[-1,]
stops_t <- as.data.frame(t(stopsOfSigGenes))


legal_stops <- c("TGA", "TAG", "TAA")

TAA <- rowSums(stops_t == "TAA", na.rm = TRUE)
TAG <- rowSums(stops_t == "TAG", na.rm = TRUE)
TGA <- rowSums(stops_t == "TGA", na.rm = TRUE)
all <- TAA + TAG + TGA
boxplot(TAA[TAA/all > 0], TAG[TAG/all > 0], TGA[TGA/all > 0], ylim = c(0,35))

stopCounts <- data.frame(TAA, TGA, TAG)
stopCounts <- stopCounts[-1,]

##lifespan stops
oldStops <- t(stopsOfSigGenes)[ t(stopsOfSigGenes$organism) %in% long,]
TAAold <- rowSums(oldStops == "TAA", na.rm = TRUE)
TAGold <- rowSums(oldStops == "TAG", na.rm = TRUE)
TGAold <- rowSums(oldStops == "TGA", na.rm = TRUE)
allold <- TAAold + TGAold + TAGold
boxplot(TAAold[TAAold/allold > 0], TGAold[TGAold/allold > 0], TAGold[TGAold/allold > 0], ylim = c(0,30))

youngStops <- t(stopsOfSigGenes)[t(stopsOfSigGenes$organism) %in% short,]
TAAyoung <- rowSums(youngStops == "TAA", na.rm = TRUE)
TAGyoung <- rowSums(youngStops == "TAG", na.rm = TRUE)
TGAyoung <- rowSums(youngStops == "TGA", na.rm = TRUE)
allyoung <- TAAyoung + TGAyoung + TAGyoung
boxplot(TAAyoung[TAAyoung/allyoung > 0], TGAyoung[TGAyoung/allyoung > 0], TAGyoung[TGAyoung/allyoung > 0], ylim = c(0,30))

######################################################################################################################
## Stops for all genes
allGenesStops <- read.csv("~/Documents/david/scripts_git/out/stops/stopsOfAllGenes.csv", header = FALSE, stringsAsFactors = FALSE)
allGenesStops[1,1] <- "organism"
colnames(allGenesStops) <- allGenesStops[1,]
allGenesStops <- allGenesStops[-1,]

allGenesStops_t <- as.data.frame(t(allGenesStops))
colnames(allGenesStops_t) <- allGenesStops_t[1,]
allGenesStops_t <- allGenesStops_t[-1,]

TAA_all <- rowSums(allGenesStops_t == "TAA", na.rm = TRUE)
TAG_all <- rowSums(allGenesStops_t == "TAG", na.rm = TRUE)
TGA_all <- rowSums(allGenesStops_t == "TGA", na.rm = TRUE)
all_all <- TAA_all + TAG_all + TGA_all
boxplot(TAA_all[TAA_all/all_all > 0], TAG_all[TAG_all/all_all > 0], TGA_all[TGA_all/all_all > 0], ylim = c(0,35))

######################################################################################################################
## dN/dS
dndsOfSigGenes <- read.csv("~/Documents/david/scripts_git/out/significant/significantdnds.csv",row.names = 1, header = TRUE, stringsAsFactors = FALSE)
# dndsOfSigGenes[1,1] <- "organism"
# colnames(dndsOfSigGenes) <- dndsOfSigGenes[1,]
# dndsOfSigGenes <- dndsOfSigGenes[-1,]
dnds_t <- as.data.frame(t(dndsOfSigGenes))
colnames(dnds_t) <- dnds_t[1,]
dnds_t <- dnds_t[-1,]

avgDnds <- colSums(dndsOfSigGenes, na.rm = TRUE)/apply(dndsOfSigGenes, 2, function(x) length(which(!is.na(x))))
na.omit(avgDnds[avgDnds > 1.5])

hist(log(avgDnds), breaks = 50, freq = FALSE)
lines(density(log(na.omit(avgDnds))), col="red")

longlivedDnds <- dndsOfSigGenes[rownames(dndsOfSigGenes) %in% long,]
avgLonglivedDnds <- colSums(longlivedDnds, na.rm = TRUE)/apply(longlivedDnds, 2, function(x) length(which(!is.na(x))))
shortlivedDnds <- dndsOfSigGenes[rownames(dndsOfSigGenes) %in% short,]
avgShortlivedDnds <- colSums(shortlivedDnds, na.rm = TRUE)/apply(shortlivedDnds, 2, function(x) length(which(!is.na(x))))

boxplot(log(avgLonglivedDnds), log(avgShortlivedDnds))

longlivedPositive <- avgLonglivedDnds[avgLonglivedDnds > 1]
longlivedPositive <- longlivedPositive[!is.na(longlivedPositive)]
shortlivedPositive <- avgShortlivedDnds[avgShortlivedDnds > 1]
shortlivedPositive <- shortlivedPositive[!is.na(shortlivedPositive)]
allPositive <- avgDnds[avgDnds > 1]
allPositive <- allPositive[!is.na(allPositive)]

onlyLongPositive <- longlivedPositive[!names(longlivedPositive) %in% names(allPositive)]
onlyShortPositive <- shortlivedPositive[!names(shortlivedPositive) %in% names(allPositive)]
shortlivedDnds$ZNF6110

intersectLongShort <- longlivedPositive[names(longlivedPositive) %in% names(shortlivedPositive)]

goodShortGenes <- onlyShortPositive[lapply(names(onlyShortPositive), function(x) sum(!is.na(shortlivedDnds[x]))) > 12]
goodLongGenes <- onlyLongPositive[lapply(names(onlyLongPositive), function(x) sum(!is.na(longlivedDnds[x]))) > 11]
goodIntersectGenes <- allPositive[lapply(names(allPositive), function(x) sum(!is.na(dndsOfSigGenes[x]))) > 32]

goodIntersectGenes[names(goodIntersectGenes) %in% ageing_genes$symbol]

shortlivedDnds$RAMP3
dndsOfSigGenes$MINOS1.NBL1

