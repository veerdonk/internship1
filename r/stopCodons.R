setwd("~/Documents/david/scripts_git/out/stops")

files <- dir("./primates/", full.names = TRUE)
titles <- c("Green monkey", "Gorilla", "Macaque", "Olive baboon", "Squirrel monkey", "Phillipine tarsier")
files <- dir("./rodents/", full.names = TRUE)
titles <- c("Brazillian guinea pig", "Long-tailed chinchilla", "Ord's kangaroo rat", "Damara mole-rat", "Golden hamster", "Mouse", "Degu", "NA deer mouse", "Rat", "Squirrel")

prim <- read.csv(files[10], header = FALSE, stringsAsFactors = FALSE)

stops <- as.factor(c("TAG", "TAA", "TGA"))


filteredDnds <- prim[prim$V5 > 0 & prim$V6 > 0.01 & prim$V5 < 2 & prim$V6 < 2 & prim$V4 < 10 & prim$V4 > 0.001 & prim$V7 %in% stops,]
filteredDnds$V7 <- as.factor(filteredDnds$V7)

boxplot(filteredDnds$V4 ~ filteredDnds$V7, data = filteredDnds, log = "y")

# pdf("stopCodonSelection.pdf")

for (file in files){
  prim <- read.csv(file, header = FALSE, stringsAsFactors = FALSE)
  
  stops <- as.factor(c("TAG", "TAA", "TGA"))

  
  filteredDnds <- prim[prim$V5 > 0 & prim$V6 > 0.01 & prim$V5 < 2 & prim$V6 < 2 & prim$V4 < 10 & prim$V4 > 0.001 & prim$V7 %in% stops,]
  filteredDnds$V7 <- as.factor(filteredDnds$V7)
  print(file)
  print(count(filteredDnds$V7))
  
  boxplot(log(filteredDnds$V4) ~ filteredDnds$V7, data = filteredDnds, main = titles[which(files == file)])
}
# dev.off()

stops <- as.factor(c("TAG", "TAA", "TGA"))
labels <- c("gene", "ref", "ort", "dnds", "dn", "ds", "stopcodon")
gorilla <- read.csv(files[3], header = FALSE, stringsAsFactors = FALSE)

mouse <- read.csv(files[6], header = FALSE, stringsAsFactors = FALSE)

colnames(gorilla) <- labels
colnames(mouse) <- labels

gorillaFiltered <- gorilla[gorilla$stopcodon %in% stops,c("gene", "stopcodon")]
mouseFiltered <- mouse[mouse$stopcodon %in% stops,c("gene", "stopcodon")]
# mouseFiltered$stopcodon <- as.factor(mouseFiltered$stopcodon)
# gorillaFiltered$stopcodon <- as.factor(gorillaFiltered$stopcodon)

gorillaMouse <- merge(gorillaFiltered, mouseFiltered, by = "gene")
gorillaMouse$stopcodon.x <- as.factor(gorillaMouse$stopcodon.x)
gorillaMouse$stopcodon.y <- as.factor(gorillaMouse$stopcodon.y)

same <- gorillaMouse[gorillaMouse$stopcodon.x == gorillaMouse$stopcodon.y,]
diff <- gorillaMouse[gorillaMouse$stopcodon.x != gorillaMouse$stopcodon.y,]
library("plyr")

count(same$stopcodon.x)
count(diff$stopcodon.x)
count(diff$stopcodon.y)

match(gorillaMouse$stopcodon.x, gorillaMouse$stopcodon.y)

rodentFiles <- dir("./rodents/", full.names = TRUE)
rodentNames <- titles <- c("Brazillian guinea pig", "Long-tailed chinchilla", "Ord's kangaroo rat", "Damara mole-rat", "Golden hamster", "Mouse", "Degu", "NA deer mouse", "Rat", "Squirrel")

labels <- c("gene", "ref", "ort", "dnds", "dn", "ds", "stopcodon")
stops <- as.factor(c("TAG", "TAA", "TGA"))
for (rodentFile in rodentFiles){
  rodentData <- read.csv(rodentFile , header = FALSE, stringsAsFactors = FALSE)
  colnames(rodentData) <- labels
  rodentDataFiltered <- rodentData[rodentData$stopcodon %in% stops, c("gene", "stopcodon")]
  colnames(rodentDataFiltered) <- c("gene", rodentNames[which(rodentFiles == rodentFile)])
  if (exists("tempMergeTable") == TRUE){
    tempMergeTable <- merge(tempMergeTable, rodentDataFiltered, by = "gene")
  }
  else{
    tempMergeTable <- rodentDataFiltered
  }
}
rodentStops <- tempMergeTable
rm(tempMergeTable)
