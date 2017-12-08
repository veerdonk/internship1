##########################################################################################
#######################################    ALL    ########################################
##########################################################################################
setwd("~/Documents/david/scripts_git/out")
labels <- c("gene", "ref", "ort", "dnds", "dn", "ds", "perc")
# 
# rodents <- dir("./names/rodents/both", full.names = TRUE)
# files <- rodents

# primates <- dir("./names/primates", full.names = TRUE)
# files <- primates

# bats <- dir("./names/bats_with_num", full.names = TRUE)
# files <- bats

# fish <- dir("./names/fish", full.names = TRUE)
# files <- fish

whales <- dir("./names/whales", full.names = TRUE)
files <- whales

# birds <- dir("./names/birds/test", full.names = TRUE)
# files <- birds

####### WITH STOP CODONS #######
# stopRodents <- dir("./stops/rodents/", full.names = TRUE) # are .csv
# stopPrim <- dir("./stops/primates/", full.names = TRUE)
# labels <- labels <- c("gene", "ref", "ort", "dnds", "dn", "ds", "stopcodon", "perc")
# files <- stopRodents

for (fileName in files){
  data <- read.csv(fileName, sep="\t", header = FALSE)
  filteredDnds <- data[data$V5 > 0 & data$V6 > 0.01 & data$V5 < 2 & data$V6 < 2 & data$V4 < 10 & data$V4 > 0.001 ,]
  filteredDnds$perc <- filteredDnds$V4/max(filteredDnds$V4)*100
  colnames(filteredDnds) <- labels
  
  if (exists("merged") == TRUE){
    dfToMerge <- data.frame(filteredDnds$gene, filteredDnds$perc)
    colnames(dfToMerge) <- c("gene", fileName)
    merged <- merge(merged, dfToMerge, by.x = "gene", by.y = "gene", all = FALSE)
  }
  else {
    merged <- data.frame(filteredDnds$gene, filteredDnds$perc)
    colnames(merged) <- c("gene", fileName)
  }
}
merged$means <- rowMeans(merged[,-1])
merged$median <- apply(merged[,2:(length(names(merged))-1)], 1, median)

# rodentPerc <- unique(merged)
# primatesPerc <- unique(merged)
# batsPerc <- unique(merged)
# fishPerc <- unique(merged)
# whalesPerc <- unique(merged)
# birdsPerc <- unique(merged)
rm(merged)


