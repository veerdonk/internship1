setwd("~/Documents/david/scripts_git/out/stops")
coln <- c("org1ID", "org2ID", "org1stop", "org2stop")
one2oneOrthologs_human_chicken  <- read.csv("./human_chicken_one2one.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
one2oneOrthologs_human_mouse    <- read.csv("./human_mouse_one2one.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(one2oneOrthologs_human_mouse)    <- coln

colnames(one2oneOrthologs_human_chicken)  <- coln

filteredStops <- filterIllegalStops(one2oneOrthologs_human_chicken)
filtrMouse <- filterIllegalStops(one2oneOrthologs_human_mouse)
getCodonFreq(filteredStops)
getCodonConservation(filteredStops)

test1 <- head(filteredStops)
test2 <- head(filtrMouse)

calculateIntersectForStop <- function(stop){
  intersection = Reduce(intersect, list(
    one2oneOrthologs_human_chicken[one2oneOrthologs_human_chicken$org1stop == stop & one2oneOrthologs_human_chicken$org2stop == stop, 1],
    one2oneOrthologs_human_mouse[one2oneOrthologs_human_mouse$org1stop == stop & one2oneOrthologs_human_mouse$org2stop == stop, 1]
  ))
  print(length(intersection))
}
calculateIntersectForStop("TAG")
# install.packages("rlist")
# library("rlist")
######################################################################################################
stopCodons <- as.factor(c("TAG", "TAA", "TGA"))

getCodonFreq <- function(x){
  stopcodonFreq <- data.frame(TAG = numeric(0), TAA = numeric(0), TGA = numeric(0))
  for (stop in stopCodons){
    stopToTga <- round(sum(x$org1stop == stop & x$org2stop == "TGA")/sum(x$org1stop == stop), 2)
    stopToTag <- round(sum(x$org1stop == stop & x$org2stop == "TAG")/sum(x$org1stop == stop), 2)
    stopToTaa <- round(sum(x$org1stop == stop & x$org2stop == "TAA")/sum(x$org1stop == stop), 2)
    stopcodonFreq[stop,] <- c(stopToTag, stopToTaa, stopToTga)
  }
  return(stopcodonFreq)
}

getCodonConservation <- function(inTable){
  conservedTga <- round(sum(inTable$org1stop == "TGA" & inTable$org2stop == "TGA")/length(inTable$org1ID), 2)
  conservedTag <- round(sum(inTable$org1stop == "TAG" & inTable$org2stop == "TAG")/length(inTable$org1ID), 2)
  conservedTaa <- round(sum(inTable$org1stop == "TAA" & inTable$org2stop == "TAA")/length(inTable$org1ID), 2)
  return(c(conservedTga, conservedTag, conservedTaa))
}

filterIllegalStops <- function(inTable){
  return(inTable[inTable$org1stop %in% stopCodons & inTable$org2stop %in% stopCodons,])
}

