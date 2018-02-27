 
#library(KEGGREST)
# install.packages("Matrix")
# install.packages("https://www.jstatsoft.org/index.php/jss/article/downloadSuppFile/v033i01/glmnet_1.1-4.tar.gz", repos=NULL, method="libcurl")
# install.packages("https://cran.r-project.org/src/contrib/doMC_1.3.5.tar.gz", repos=NULL, method="libcurl")

library(GenomicRanges)
library(glmnet)
library(doMC)
library(gplots)
library( matrixStats)

# library(ggplot2)
# library(viridis)
# require(cowplot)

# Load a bunch of colors
Tcolors  <- colorRampPalette(c("cyan","black", "gold"))(n = 99)
RBcolors <- colorRampPalette(c("red","white", "blue"))(n = 99)
myColors = c("black","green","lightblue","orange","red","blue")

# Load data
allData = read.table("/home/dylan/data/orthologs/allOrganismsOrthologyTable.csv", sep = ",", head = T, row.names = 2, stringsAsFactors  = T)
data_matrix = allData[,-1]

# Phenotype data
phenoData = read.table("data/AnimalGroups.csv",sep = ",", head = T)
phenoData$colors = myColors[phenoData$Group] 
rownames(phenoData) = phenoData$Organism
phenoData$AgeColor = Tcolors[round(phenoData$Age)+1]                #Add some color

#Ageing data
ageingData = read.table("data/ageing_genes.csv", sep = ",", head = T, stringsAsFactors  = T)


#Align names
phenoData = phenoData[colnames(data_matrix),]
colnames(data_matrix ) == rownames(phenoData)

# Order by age... just to see
ageOrder = order(phenoData$Age)
data_matrix = data_matrix[,ageOrder]
phenoData = phenoData[ageOrder,]

colnames(data_matrix ) == rownames(phenoData)

head(data_matrix )
head(phenoData)

test =  data_matrix
new =ifelse(test == "missing", 1, 1)
new =ifelse(test == "one2many", 10, new)
new =ifelse(test == "one2one",  10, new)
new =ifelse(test == "many2many",20,new)
new =ifelse(test == "many2one", 20,new)

unique(test[,1])
filter = apply(new, 1, function(x) var(x) != 0)
keep = new[filter,]

# Then calculate how far it is off from the truth for each sample in SD. (Probably ranked)
# Sort by those who have the least fail of all of them.


#########################################################################
## Reversed approach
#########################################################################

# truePositives = intersect(rownames(keep) ,ageingData[,1])
# 
# heatmap.2(keep[truePositives,], 
#           main = "Known ageing genes",
#           labCol = paste(phenoData$Organism, phenoData$Age, sep = ' '), 
#           Rowv = F, Colv = T, 
#           trace = "none",
#           #colCol = phenoData$AgeColor,
#           #RowSideColors = detected,
#           ColSideColors = phenoData$colors)
# 
# head(test)
# 
# intersect(names(keep),ageingData[,1])
# 
# heatmap.2(keep[names(first100),mostPowerfullColumns], 
#           main = localHeader,
#           labCol = paste(phenoData$Organism[mostPowerfullColumns], phenoData$Age[mostPowerfullColumns], sep = ' '), 
#           Rowv = F, Colv = F, 
#           trace = "none",
#           #colCol = phenoData$AgeColor,
#           RowSideColors = detected,
#           ColSideColors = phenoData$colors[mostPowerfullColumns])


#########################################################################
## Average age calculations
#########################################################################
test =  data_matrix
new =ifelse(test == "missing", 1, 0)
new =ifelse(test == "one2many", 2, new)
new =ifelse(test == "one2one",  2, new)
new =ifelse(test == "many2many",3,new)
new =ifelse(test == "many2one", 3,new)

unique(test[,1])
filter = apply(new, 1, function(x) var(x) != 0)
relations = new[filter,]

# 1. Take the average age of each gene relation, within species.
meanAgeFish = mean(phenoData$Age[phenoData$Group == "Fish"])
medianAgeFish = median(phenoData$Age[phenoData$Group == "Fish"])

meanAgeWhale = mean(phenoData$Age[phenoData$Group == "Whale"])
medianAgeWhale = median(phenoData$Age[phenoData$Group == "Whale"])

meanAgeBirds = mean(phenoData$Age[phenoData$Group == "Birds"])
medianAgeBirds = median(phenoData$Age[phenoData$Group == "Birds"])

meanAgeBats = mean(phenoData$Age[phenoData$Group == "Bats"])
medianAgeBats = median(phenoData$Age[phenoData$Group == "Bats"])

meanAgeRodents = mean(phenoData$Age[phenoData$Group == "Rodents"])
medianAgeRodents = median(phenoData$Age[phenoData$Group == "Rodents"])

meanAgePrimates = mean(phenoData$Age[phenoData$Group == "Primates"])
medianAgePrimates = median(phenoData$Age[phenoData$Group == "Primates"])


FishFrame = relations[,1:3]
colnames(FishFrame) = c("missing","one", "many")

WhaleFrame = relations[,1:3]
colnames(WhaleFrame) = c("missing","one", "many")

BirdsFrame = relations[,1:3]
colnames(BirdsFrame) = c("missing","one", "many")

BatsFrame = relations[,1:3]
colnames(BatsFrame) = c("missing","one", "many")

RodentsFrame = relations[,1:3]
colnames(RodentsFrame) = c("missing","one", "many")

PrimatesFrame = relations[,1:3]
colnames(PrimatesFrame) = c("missing","one", "many")



for (i in 1:nrow(scoreFrame)){
  print(i)
  
  ############## Fish
  # Missing
  missingFish = mean(phenoData$Age[relations[i,] == 1 & phenoData$Group == "Fish" ])
  FishFrame[i,1] = missingFish - meanAgeFish
  # One
  oneFish = mean(phenoData$Age[relations[i,] == 2 & phenoData$Group == "Fish" ])
  FishFrame[i,2] = oneFish - meanAgeFish
  # Many
  manyFish = mean(phenoData$Age[relations[i,] == 3 & phenoData$Group == "Fish" ])
  FishFrame[i,3] = manyFish - meanAgeFish
  
  ############## Whale
  # Missing
  missingWhale = mean(phenoData$Age[relations[i,] == 1 & phenoData$Group == "Whale" ])
  WhaleFrame[i,1] = missingWhale - meanAgeWhale
  # One
  oneWhale = mean(phenoData$Age[relations[i,] == 2 & phenoData$Group == "Whale" ])
  WhaleFrame[i,2] = oneWhale - meanAgeWhale
  # Many
  manyWhale = mean(phenoData$Age[relations[i,] == 3 & phenoData$Group == "Whale" ])
  WhaleFrame[i,3] = manyWhale - meanAgeWhale
  
  ############## Birds
  # Missing
  missingBirds = mean(phenoData$Age[relations[i,] == 1 & phenoData$Group == "Birds" ])
  BirdsFrame[i,1] = missingBirds - meanAgeBirds
  # One
  oneBirds = mean(phenoData$Age[relations[i,] == 2 & phenoData$Group == "Birds" ])
  BirdsFrame[i,2] = oneBirds - meanAgeBirds
  # Many
  manyBirds = mean(phenoData$Age[relations[i,] == 3 & phenoData$Group == "Birds" ])
  BirdsFrame[i,3] = manyBirds - meanAgeBirds
  
  ############## Bats
  # Missing
  missingBats = mean(phenoData$Age[relations[i,] == 1 & phenoData$Group == "Bats" ])
  BatsFrame[i,1] = missingBats - meanAgeBats
  # One
  oneBats = mean(phenoData$Age[relations[i,] == 2 & phenoData$Group == "Bats" ])
  BatsFrame[i,2] = oneBats - meanAgeBats
  # Many
  manyBats = mean(phenoData$Age[relations[i,] == 3 & phenoData$Group == "Bats" ])
  BatsFrame[i,3] = manyBats - meanAgeBats
  
  ############## Rodents
  # Missing
  missingRodents = mean(phenoData$Age[relations[i,] == 1 & phenoData$Group == "Rodents" ])
  RodentsFrame[i,1] = missingRodents - meanAgeRodents
  # One
  oneRodents = mean(phenoData$Age[relations[i,] == 2 & phenoData$Group == "Rodents" ])
  RodentsFrame[i,2] = oneRodents - meanAgeRodents
  # Many
  manyRodents = mean(phenoData$Age[relations[i,] == 3 & phenoData$Group == "Rodents" ])
  RodentsFrame[i,3] = manyRodents - meanAgeRodents
  
  ############## Primates
  # Missing
  missingPrimates = mean(phenoData$Age[relations[i,] == 1 & phenoData$Group == "Primates" ])
  PrimatesFrame[i,1] = missingPrimates - meanAgePrimates
  # One
  onePrimates = mean(phenoData$Age[relations[i,] == 2 & phenoData$Group == "Primates" ])
  PrimatesFrame[i,2] = onePrimates - meanAgePrimates
  # Many
  manyPrimates = mean(phenoData$Age[relations[i,] == 3 & phenoData$Group == "Primates" ])
  PrimatesFrame[i,3] = manyPrimates - meanAgePrimates
}


############### Fish Plot ######################################################################
plot( density(na.omit(FishFrame[,1])), col = "red", main = "mean ages\n Fish", xlim = c(-5,5))
lines(density(na.omit(FishFrame[,2])), col = "black")
lines(density(na.omit(FishFrame[,3])), col = "forestgreen")
abline(v = 0, lty = 2)

############### Whale Plot ######################################################################
plot( density(na.omit(WhaleFrame[,1])), col = "red", main = "mean ages\n Whale", xlim = c(-35,35))
lines(density(na.omit(WhaleFrame[,2])), col = "black")
lines(density(na.omit(WhaleFrame[,3])), col = "forestgreen")
abline(v = 0, lty = 2)

############### Birds Plot ######################################################################
plot( density(na.omit(BirdsFrame[,1])), col = "red", main = "mean ages\n Birds", xlim = c(-15,15))
lines(density(na.omit(BirdsFrame[,2])), col = "black")
lines(density(na.omit(BirdsFrame[,3])), col = "forestgreen")
abline(v = 0, lty = 2)

############### Bats Plot ######################################################################
plot( density(na.omit(BatsFrame[,1])), col = "red", main = "mean ages\n Bats", xlim = c(-15,15))
lines(density(na.omit(BatsFrame[,2])), col = "black")
lines(density(na.omit(BatsFrame[,3])), col = "forestgreen")
abline(v = 0, lty = 2)

############### Rodents Plot ######################################################################
plot( density(na.omit(RodentsFrame[,1])), col = "red", main = "mean ages\n Rodents", xlim = c(-20,20))
lines(density(na.omit(RodentsFrame[,2])), col = "black")
lines(density(na.omit(RodentsFrame[,3])), col = "forestgreen")
abline(v = 0, lty = 2)

############### Primates Plot ######################################################################
plot( density(na.omit(PrimatesFrame[,1])), col = "red", main = "mean ages\n Primates", xlim = c(-25,25))
lines(density(na.omit(PrimatesFrame[,2])), col = "black")
lines(density(na.omit(PrimatesFrame[,3])), col = "forestgreen")
abline(v = 0, lty = 2)



#####################################################################################################

# 2. Take the average age of each gene relation, overal.

#####################################################################################################
meanAge = mean(phenoData$Age)
medianAge = median(phenoData$Age)

scoreFrame = relations[,1:3]
colnames(scoreFrame) = c("missing","one", "many")

for (i in 1:nrow(scoreFrame)){
  print(i)
  # Missing
  Missing = mean(phenoData$Age[relations[i,] == 1])
  scoreFrame[i,1] = Missing - meanAge
  
  # One
  One     = mean(phenoData$Age[relations[i,] == 2])
  One - meanAge
  scoreFrame[i,2] = One - meanAge
  
  # Many
  Many= mean(phenoData$Age[relations[i,] == 3])
  scoreFrame[i,3] = Many - meanAge
}


plot( density(na.omit(scoreFrame[,1])), col = "red", main = "All ages\n All species",xlim = c(-50,50)  )
lines(density(na.omit(scoreFrame[,2])), col = "black")
lines(density(na.omit(scoreFrame[,3])), col = "forestgreen")
abline(v = 0, lty = 2)



# 3. How much is the average age greater than the other relations?

# 3b. How much is the average greater than the overal average?

# 4. Give gene + relation a score.

# 5. Rank by information score?



