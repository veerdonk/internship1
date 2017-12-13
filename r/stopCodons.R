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
  dnds.aov <- aov(log(filteredDnds$V4) ~ filteredDnds$V7, data = filteredDnds)
  print(summary(dnds.aov))
  # plot(dnds.aov, 2)
  print(count(filteredDnds$V7))
  
  boxplot(log(filteredDnds$V4) ~ filteredDnds$V7, data = filteredDnds, main = titles[which(files == file)])
}
# dev.off()


test <- PlantGrowth
pairwise.t.test(filteredDnds$V4, filteredDnds$V7, p.adjust.method = "BH")

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
# rodentNames <- c("Green monkey", "Gorilla", "Macaque", "Olive baboon", "Squirrel monkey", "Phillipine tarsier")

labels <- c("gene", "ref", "ort", "dnds", "dn", "ds", "stopcodon")
stops <- as.factor(c("TAG", "TAA", "TGA"))
rodentFiles <- rodentFiles[1]
for (rodentFile in rodentFiles){
  rodentData <- read.csv(rodentFile , header = FALSE, stringsAsFactors = FALSE)
  colnames(rodentData) <- labels
  rodentDataFiltered <- rodentData[rodentData$stopcodon %in% stops, c("gene", "stopcodon")]# & rodentData$dn > 0 & rodentData$ds > 0.01 & rodentData$dn < 2 & rodentData$ds < 2
  # rodentDataFiltered[,2] <- as.factor(rodentDataFiltered[,2])
  rodentDataFiltered[,1] <- as.factor(rodentDataFiltered[,1])
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

primateFiles <- dir("./primates/", full.names = TRUE)
file1 <- read.csv(primateFiles[3], header = FALSE, stringsAsFactors = FALSE)
colnames(file1) <- labels
file1 <- file1[file1$stopcodon %in% stops, c("gene", "stopcodon")]
# file1$gene <- as.factor(file1$gene)
# file1$stopcodon <- as.factor(file1$stopcodon)


overlap <- rodentStops[rodentStops$gene %in% file1$gene,]
overlap2 <- file1[file1$gene %in% rodentStops$gene,]

calculateFrequencies <- function(x){
  # same <- sum(overlap[overlap$gene == x[1],][-1] == x[2])/length(overlap[overlap$gene == x[1],][-1])*100
  # otherStops <- which(stops != x[2])

  stop1 <- sum(overlap[overlap$gene == x[1],][-1] == as.character(stops[1]))/length(overlap[overlap$gene == x[1],][-1])*100
  stop2 <- sum(overlap[overlap$gene == x[1],][-1] == as.character(stops[2]))/length(overlap[overlap$gene == x[1],][-1])*100
  stop3 <- sum(overlap[overlap$gene == x[1],][-1] == as.character(stops[3]))/length(overlap[overlap$gene == x[1],][-1])*100
  c(x[1], x[2], stop1, stop2, stop3)
}

stopcodonFreq <- data.frame(t(apply(overlap2, 1, calculateFrequencies)), row.names = "gene")
colnames(stopcodonFreq) <- c("primateStop", as.character(stops[1]), as.character(stops[2]), as.character(stops[3]))

functionName <- function(x){
  as.numeric(x[which(as.character(stops) == x[1])+1])
}

stopcodonFreq[stopcodonFreq$primateStop == "TAA",]

sum(apply(stopcodonFreq[stopcodonFreq$primateStop == "TAA",], 1, functionName))/length(stopcodonFreq[stopcodonFreq$primateStop == "TAA",1])
sum(apply(stopcodonFreq[stopcodonFreq$primateStop == "TAG",], 1, functionName))/length(stopcodonFreq[stopcodonFreq$primateStop == "TAG",1])
sum(apply(stopcodonFreq[stopcodonFreq$primateStop == "TGA",], 1, functionName))/length(stopcodonFreq[stopcodonFreq$primateStop == "TGA",1])

hk_genes <- read.table("/home/dylan/data/genes/HK_genes.txt", stringsAsFactors = FALSE)
ribosomal_genes <- read.table("/home/dylan/data/genes/ribosomal_genes.txt", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
ribo_freq <- stopcodonFreq[row.names(stopcodonFreq) %in% ribosomal_genes$Approved.Symbol,]
hk_freq <- stopcodonFreq[row.names(stopcodonFreq) %in% hk_genes$V1,]

taa <- c()
tag <- c()
tga <- c()

hk_taa <- c()
hk_tag <- c()
hk_tga <- c()

hk_freqs <- c()
all_freqs <- c()
primateFiles <- dir("./primates/", full.names = TRUE)
for (primateFile in primateFiles){
  primate <- read.csv(primateFile, header = FALSE, stringsAsFactors = FALSE)
  colnames(primate) <- labels
  primate <- primate[primate$stopcodon %in% stops, c("gene", "stopcodon")]

  overlap <- rodentStops[rodentStops$gene %in% primate$gene,]
  overlap2 <- primate[primate$gene %in% rodentStops$gene,]
  
  stopcodonFreq <- data.frame(t(apply(overlap2, 1, calculateFrequencies)), row.names = "gene")
  colnames(stopcodonFreq) <- c("primateStop", as.character(stops[1]), as.character(stops[2]), as.character(stops[3]))
  
  hk_freq <- stopcodonFreq[row.names(stopcodonFreq) %in% hk_genes$V1,]
  
  taa <- c(taa, sum(apply(stopcodonFreq[stopcodonFreq$primateStop == "TAA",], 1, functionName))/length(stopcodonFreq[stopcodonFreq$primateStop == "TAA",1]))
  tag <- c(tag, sum(apply(stopcodonFreq[stopcodonFreq$primateStop == "TAG",], 1, functionName))/length(stopcodonFreq[stopcodonFreq$primateStop == "TAG",1]))
  tga <- c(tga, sum(apply(stopcodonFreq[stopcodonFreq$primateStop == "TGA",], 1, functionName))/length(stopcodonFreq[stopcodonFreq$primateStop == "TGA",1]))
  
  hk_taa <- c(hk_taa, sum(apply(hk_freq[hk_freq$primateStop == "TAA",], 1, functionName))/length(hk_freq[hk_freq$primateStop == "TAA",1]))
  hk_tag <- c(hk_tag, sum(apply(hk_freq[hk_freq$primateStop == "TAG",], 1, functionName))/length(hk_freq[hk_freq$primateStop == "TAG",1]))
  hk_tga <- c(hk_tga, sum(apply(hk_freq[hk_freq$primateStop == "TGA",], 1, functionName))/length(hk_freq[hk_freq$primateStop == "TGA",1]))
  
  one <- sum(hk_freq$primateStop == "TAA")/length(hk_freq$primateStop)*100
  two <- sum(hk_freq$primateStop == "TAG")/length(hk_freq$primateStop)*100
  three <- sum(hk_freq$primateStop == "TGA")/length(hk_freq$primateStop)*100
  
  hk_freqs <- c(hk_freqs, c(one, two, three))
  one <- sum(stopcodonFreq$primateStop == "TAA")/length(stopcodonFreq$primateStop)*100
  two <- sum(stopcodonFreq$primateStop == "TAG")/length(stopcodonFreq$primateStop)*100
  three <- sum(stopcodonFreq$primateStop == "TGA")/length(stopcodonFreq$primateStop)*100
  all_freqs <- c(hk_freqs, c(one, two, three))
}

setwd("~/Documents/david/scripts_git/out/stops")
stops <- as.factor(c("TAG", "TAA", "TGA"))
cavia <- read.csv("./rodents/rno.csv", header = FALSE, stringsAsFactors = FALSE)
monkey <- read.csv("./primates/ggo.csv", header = FALSE, stringsAsFactors = FALSE)
cavia <- cavia[cavia$V7 %in% stops,]
monkey <- monkey[monkey$V7 %in% stops,]


caviaStops <- cavia[,c(1,7)]
monkeyStops <- monkey[,c(1,7)]

colnames(caviaStops) <- c("gene", "stopCavia")
colnames(monkeyStops) <- c("gene", "stopMonkey")

#combine primate + rodent
stopsMerged <- merge(monkeyStops, caviaStops, by = "gene")

#number of stopcodons that are the same
tagSame <- length(stopsMerged[stopsMerged$stopMonkey == "TAG" & stopsMerged$stopCavia == "TAG",1])
tgaSame <- length(stopsMerged[stopsMerged$stopMonkey == "TGA" & stopsMerged$stopCavia == "TGA",1])
taaSame <- length(stopsMerged[stopsMerged$stopMonkey == "TAA" & stopsMerged$stopCavia == "TAA",1])

#percentage of taa that is same
taaSame/length(stopsMerged[stopsMerged$stopCavia== "TAA",1])*100
taaSame/length(stopsMerged[stopsMerged$stopMonkey== "TAA",1])*100

#of tga/tag
tgaSame/length(stopsMerged[stopsMerged$stopCavia == "TGA",1])*100
tgaSame/length(stopsMerged[stopsMerged$stopMonkey == "TGA",1])*100

tagSame/length(stopsMerged[stopsMerged$stopCavia == "TAG",1])*100
tagSame/length(stopsMerged[stopsMerged$stopMonkey == "TAG",1])*100

#number of TAA -> TGA/TAG for monkey
paste("monkey TAA -> TAG:", length(stopsMerged[stopsMerged$stopMonkey == "TAA" & stopsMerged$stopCavia == "TAG",1]), "=", length(stopsMerged[stopsMerged$stopMonkey == "TAA" & stopsMerged$stopCavia == "TAG",1])/length(stopsMerged[stopsMerged$stopMonkey == "TAA",1])*100, "%")
paste("monkey TAA -> TGA:", length(stopsMerged[stopsMerged$stopMonkey == "TAA" & stopsMerged$stopCavia == "TGA",1]), "=", length(stopsMerged[stopsMerged$stopMonkey == "TAA" & stopsMerged$stopCavia == "TGA",1])/length(stopsMerged[stopsMerged$stopMonkey == "TAA",1])*100, "%")

#number of TAA -> TGA/TAG for cavia
paste("monkey TAA -> TAG:", length(stopsMerged[stopsMerged$stopCavia== "TAA" & stopsMerged$stopMonkey == "TAG",1]), "=", length(stopsMerged[stopsMerged$stopCavia== "TAA" & stopsMerged$stopMonkey == "TAG",1])/length(stopsMerged[stopsMerged$stopCavia== "TAA",1])*100, "%")
paste("monkey TAA -> TGA:", length(stopsMerged[stopsMerged$stopCavia== "TAA" & stopsMerged$stopMonkey == "TGA",1]), "=", length(stopsMerged[stopsMerged$stopCavia== "TAA" & stopsMerged$stopMonkey == "TGA",1])/length(stopsMerged[stopsMerged$stopCavia== "TAA",1])*100, "%")

