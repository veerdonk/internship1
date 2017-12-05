setwd("~/Documents/david/scripts_git/out")
labels <- c("gene", "ref", "ort", "dnds", "dn", "ds", "perc")

##########################################################################################
#######################################  PRIMATES  #######################################
##########################################################################################
primateFiles <- dir("./names/primates", full.names = TRUE)

for (primate in primateFiles){
  primateData <- read.csv(primate, sep = "\t", header = FALSE)
  primateDataDnsdsFixed <- primateData[primateData$V5>0 & primateData$V6 >0.001 & primateData$V5<2 & primateData$V6<2 & primateData$V4<10 & primateData$V4>0.001,]
  primateDataDnsdsFixed$V7 <- primateDataDnsdsFixed$V4/max(primateDataDnsdsFixed$V4)*100
  colnames(primateDataDnsdsFixed) <- labels
  
  if (exists("primatePerc") == TRUE){
    dftomerge <- data.frame(primateDataDnsdsFixed$gene, primateDataDnsdsFixed$perc)
    colnames(dftomerge) <- c("gene", "perc")
    primatePerc <- merge(primatePerc, dftomerge, by = "gene", all = TRUE)
   
  }else{
    primatePerc <- data.frame(primateDataDnsdsFixed$gene, primateDataDnsdsFixed$perc)
    colnames(primatePerc) <- c("gene", "perc")
  }
}
primatePerc$means <- rowMeans(primatePerc[,-1]) 
primatePerc$median <- apply(primatePerc[,2:(length(names(primatePerc))-1)], 1, median)


##########################################################################################
#######################################  RODENTS  ########################################
##########################################################################################
rodentFiles <- dir("./names/rodents/both", full.names = TRUE)

for (file in rodentFiles){
  data <- read.csv(file, sep="\t", header = FALSE)
  filteredDnds <- data[data$V5 > 0 & data$V6 > 0.01 & data$V5 < 2 & data$V6 < 2 & data$V4 < 10 & data$V4 > 0.001 ,]
  filteredDnds$V7 <- filteredDnds$V4/max(filteredDnds$V4)*100
  colnames(filteredDnds) <- labels
  
  if (exists("merged") == TRUE){
    dfToMerge <- data.frame(filteredDnds$gene, filteredDnds$perc)
    colnames(dfToMerge) <- c("gene", file)
    merged <- merge(merged, dfToMerge, by.x = "gene", by.y = "gene", all = TRUE)
  }
  else {
    merged <- data.frame(filteredDnds$gene, filteredDnds$perc)
    colnames(merged) <- c("gene", "perc")
  }
}
test <- merge(merged, dfToMerge, by = "gene", all = TRUE)
merged$means <- rowMeans(merged[,-1])
merged$median <- apply(merged[,2:(length(names(merged))-1)], 1, median)

## change this ##
rodentPercentiles <- unique(merged)
rm(merged)

