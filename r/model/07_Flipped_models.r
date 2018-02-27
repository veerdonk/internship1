Pcut = 0.0001
#library(KEGGREST)

#install.packages("https://www.jstatsoft.org/index.php/jss/article/downloadSuppFile/v033i01/glmnet_1.1-4.tar.gz", repos=NULL, method="libcurl")
#install.packages("https://cran.r-project.org/src/contrib/doMC_1.3.5.tar.gz", repos=NULL, method="libcurl")

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

#write.csv(rownames(allData), "allgenes.csv", quote = TRUE, row.names = FALSE)

# Phenotype data
phenoData2 =  read.table("/home/dylan/Documents/david/scripts_git/internship1/r/model/animalGroups2.csv",sep = ",", head = T)
rownames(phenoData2) = phenoData2$Organism

# Remove that tarsier
# allData =    allData[,  -which(colnames(allData)=="soe"   | colnames(allData)=="tsy")]
# phenoData2 = phenoData2[-which(rownames(phenoData2)=="soe"| rownames(phenoData2)=="tsy"),]

#####################################################################################
# Sorting everything out

phenoData2 = phenoData2[order(phenoData2$Age, decreasing = FALSE),]
newAge =c(as.character(phenoData2$Organism[phenoData2$Group == "Fish"]),
          as.character(phenoData2$Organism[phenoData2$Group == "Birds"]),
          as.character(phenoData2$Organism[phenoData2$Group == "Whale"]),
          as.character(phenoData2$Organism[phenoData2$Group == "Bats"]),
          as.character(phenoData2$Organism[phenoData2$Group == "Rodents"]),
          as.character(phenoData2$Organism[phenoData2$Group == "Primates"]))

phenoData2 = phenoData2[newAge,]

phenoData2$colors = ifelse(phenoData2$Group == "Fish",     "lightblue", "grey")
phenoData2$colors = ifelse(phenoData2$Group == "Birds",    "green",  phenoData2$colors)
phenoData2$colors = ifelse(phenoData2$Group == "Whale",    "blue",   phenoData2$colors)
phenoData2$colors = ifelse(phenoData2$Group == "Bats",     "black",  phenoData2$colors)
phenoData2$colors = ifelse(phenoData2$Group == "Rodents",  "red",    phenoData2$colors)
phenoData2$colors = ifelse(phenoData2$Group == "Primates", "Orange", phenoData2$colors)


#####################################################################################
# Creating Ranks!

phenoData2$localRank = phenoData2$Age
phenoData2$localRank[phenoData2$Group == "Fish"] =     rank(phenoData2$Age[phenoData2$Group == "Fish"])
phenoData2$localRank[phenoData2$Group == "Whale"] =    rank(phenoData2$Age[phenoData2$Group == "Whale"])
phenoData2$localRank[phenoData2$Group == "Birds"] =    rank(phenoData2$Age[phenoData2$Group == "Birds"])
phenoData2$localRank[phenoData2$Group == "Bats"] =     rank(phenoData2$Age[phenoData2$Group == "Bats"])
phenoData2$localRank[phenoData2$Group == "Rodents"] =  rank(phenoData2$Age[phenoData2$Group == "Rodents"])
phenoData2$localRank[phenoData2$Group == "Primates"] = rank(phenoData2$Age[phenoData2$Group == "Primates"])


#####################################################################################
# Retrieving the new organisms into the mix!
## IF all are against the same reference, how can I have a 2 and a 5 in the same row?!
new = ifelse(allData == "missing",   1, 6)
new = ifelse(allData == "one2one",   2, new)
new = ifelse(allData == "one2many",  3, new)
new = ifelse(allData == "many2one",  4, new)
new = ifelse(allData == "many2many", 5, new)
new = new[,-1]

# 1. Retrieve the Fish reference
mostFish = rownames(phenoData2[ phenoData2[,2] == "Fish" & phenoData2[,5] == "no",])
refFish = rownames(phenoData2[ phenoData2[,2] == "Fish" & phenoData2[,5] == "yes",])
selected = colnames(new) %in% mostFish
newFish = ifelse (apply(new[,selected], 1, max) == 2, 2, 1)  
newFish = ifelse (apply(new[,selected], 1, max) == 3, 2, newFish)
newFish = ifelse (apply(new[,selected], 1, max) == 4, 5, newFish)
newFish = ifelse (apply(new[,selected], 1, max) == 5, 5, newFish)
new = cbind(new, newFish)
colnames(new)[ncol(new)] = refFish

# 2. Retrieve the Birds reference
mostBirds = rownames(phenoData2[ phenoData2[,2] == "Birds" & phenoData2[,5] == "no",])
refBirds = rownames(phenoData2[ phenoData2[,2] == "Birds" & phenoData2[,5] == "yes",])
selected = colnames(new) %in% mostBirds
newBirds = ifelse (apply(new[,selected], 1, max) == 2, 2, 1)  
newBirds = ifelse (apply(new[,selected], 1, max) == 3, 2, newBirds)
newBirds = ifelse (apply(new[,selected], 1, max) == 4, 5, newBirds)
newBirds = ifelse (apply(new[,selected], 1, max) == 5, 5, newBirds)
new = cbind(new, newBirds)
colnames(new)[ncol(new)] = refBirds

# 3. Retrieve the Whale reference
mostWhale = rownames(phenoData2[ phenoData2[,2] == "Whale" & phenoData2[,5] == "no",])
refWhale = rownames(phenoData2[ phenoData2[,2] == "Whale" & phenoData2[,5] == "yes",])
selected = colnames(new) %in% mostWhale
newWhale = ifelse (apply(new[,selected], 1, max) == 2, 2, 1)  
newWhale = ifelse (apply(new[,selected], 1, max) == 3, 2, newWhale)
newWhale = ifelse (apply(new[,selected], 1, max) == 4, 5, newWhale)
newWhale = ifelse (apply(new[,selected], 1, max) == 5, 5, newWhale)
new = cbind(new, newWhale)
colnames(new)[ncol(new)] = refWhale

# 4. Retrieve the Bats reference
mostBats = rownames(phenoData2[ phenoData2[,2] == "Bats" & phenoData2[,5] == "no",])
refBats = rownames(phenoData2[ phenoData2[,2] == "Bats" & phenoData2[,5] == "yes",])
selected = colnames(new) %in% mostBats
newBats = ifelse (apply(new[,selected], 1, max) == 2, 2, 1)  
newBats = ifelse (apply(new[,selected], 1, max) == 3, 2, newBats)
newBats = ifelse (apply(new[,selected], 1, max) == 4, 5, newBats)
newBats = ifelse (apply(new[,selected], 1, max) == 5, 5, newBats)
new = cbind(new, newBats)
colnames(new)[ncol(new)] = refBats

# 5. Retrieve the Rodents reference
mostRodents = rownames(phenoData2[ phenoData2[,2] == "Rodents" & phenoData2[,5] == "no",])
refRodents = rownames(phenoData2[ phenoData2[,2] == "Rodents" & phenoData2[,5] == "yes",])
selected = colnames(new) %in% mostRodents
newRodents = ifelse (apply(new[,selected], 1, max) == 2, 2, 1)  
newRodents = ifelse (apply(new[,selected], 1, max) == 3, 2, newRodents)
newRodents = ifelse (apply(new[,selected], 1, max) == 4, 5, newRodents)
newRodents = ifelse (apply(new[,selected], 1, max) == 5, 5, newRodents)
new = cbind(new, newRodents)
colnames(new)[ncol(new)] = refRodents

# 6. Retrieve the Primates reference
mostPrimates = rownames(phenoData2[ phenoData2[,2] == "Primates" & phenoData2[,5] == "no",])
refPrimates = rownames(phenoData2[ phenoData2[,2] == "Primates" & phenoData2[,5] == "yes",])
selected = colnames(new) %in% mostPrimates
newPrimates = ifelse (apply(new[,selected], 1, max) == 2, 2, 1)  
newPrimates = ifelse (apply(new[,selected], 1, max) == 3, 2, newPrimates)
newPrimates = ifelse (apply(new[,selected], 1, max) == 4, 5, newPrimates)
newPrimates = ifelse (apply(new[,selected], 1, max) == 5, 5, newPrimates)
new = cbind(new, newPrimates)
colnames(new)[ncol(new)] = refPrimates

data_matrix = new

colnames(data_matrix)[23] = "mph"

#Align names
data_matrix = data_matrix[, rownames(phenoData2)]
## I hereby see that the phenoData matches data_matrix
colnames(data_matrix ) == rownames(phenoData2)

phenoData2$youngOld = as.character(phenoData2$Full.Name)

phenoData2$youngOld[phenoData2$Group == "Fish"] = c(rep("young",5),rep("old",4))
phenoData2$youngOld[phenoData2$Group == "Birds"] = c(rep("young",6),rep("old",7))
phenoData2$youngOld[phenoData2$Group == "Whale"] = c(rep("young",3),rep("old",2))
phenoData2$youngOld[phenoData2$Group == "Bats"] = c(rep("young",3),rep("old",2))
phenoData2$youngOld[phenoData2$Group == "Rodents"] = c(rep("young",6),rep("old",7))
phenoData2$youngOld[phenoData2$Group == "Primates"] = c(rep("young",3),rep("old",4))


#####################################################################################
## Simplify the relations in the matrix

new = ifelse(allData == "missing",   1, 6)
new = ifelse(allData == "one2one",   2, new)
new = ifelse(allData == "one2many",  3, new)
new = ifelse(allData == "many2one",  4, new)
new = ifelse(allData == "many2many", 5, new)

simple =ifelse(data_matrix == 1, 1, 0)
simple =ifelse(data_matrix == 2, 2, simple)
simple =ifelse(data_matrix == 3, 3, simple)
simple =ifelse(data_matrix == 4, 2, simple)
simple =ifelse(data_matrix == 5, 3, simple)

filter = apply(simple, 1, function(x) var(x) != 0)
relations = simple[filter,]

#Ageing data
ageingData = read.table("/home/dylan/Documents/david/resources/ageing_genes/ageing_genes.csv", sep = ",", head = T, stringsAsFactors  = T)


#########################################################################
## Average age calculations
#########################################################################

FishFrame = relations[,1:2]
colnames(FishFrame) = c("Pval", "Adj-Pval")

WhaleFrame = relations[,1:2]
colnames(WhaleFrame) = c("Pval", "Adj-Pval")

BirdsFrame = relations[,1:2]
colnames(BirdsFrame) = c("Pval", "Adj-Pval")

BatsFrame = relations[,1:2]
colnames(BatsFrame) = c("Pval", "Adj-Pval")

RodentsFrame = relations[,1:2]
colnames(RodentsFrame) = c("Pval", "Adj-Pval")

PrimatesFrame = relations[,1:2]
colnames(PrimatesFrame) = c("Pval", "Adj-Pval")

AllFrame = relations[,1:2]
colnames(AllFrame) = c("Pval", "Adj-Pval")

NotfishFrame = relations[,1:2]
colnames(NotfishFrame) = c("Pval", "Adj-Pval")

NotfishbirdsFrame = relations[,1:2]
colnames(NotfishbirdsFrame) = c("Pval", "Adj-Pval")



for (i in 1:nrow(FishFrame)){
    print(i)
  
    ############## Fish ##############################################################
    youngFish = relations[i,phenoData2$Group == "Fish" & phenoData2$youngOld == "young"]
    oldFish = relations[i,phenoData2$Group == "Fish" & phenoData2$youngOld == "old"]
    
    Pval = wilcox.test(x = youngFish,
                       y = oldFish, alternative = "two.sided", 
                       paired =FALSE,
                       exact = FALSE)
    FishFrame[i,1] = Pval$p.value
    FishFrame[i,2] = Pval$p.value*nrow(relations)
    ############## Whale ##############################################################
    youngWhale = relations[i,phenoData2$Group == "Whale" & phenoData2$youngOld == "young"]
    oldWhale = relations[i,phenoData2$Group == "Whale" & phenoData2$youngOld == "old"]
    
    Pval = wilcox.test(x = youngWhale,
                       y = oldWhale, alternative = "two.sided", 
                       paired =FALSE,
                       exact = FALSE)
    WhaleFrame[i,1] = Pval$p.value
    WhaleFrame[i,2] = Pval$p.value*nrow(relations)
  

    ############## Birds ##############################################################
    youngBirds = relations[i,phenoData2$Group == "Birds" & phenoData2$youngOld == "young"]
    oldBirds = relations[i,phenoData2$Group == "Birds" & phenoData2$youngOld == "old"]
    
    Pval = wilcox.test(x = youngBirds,
                       y = oldBirds, alternative = "two.sided", 
                       paired =FALSE,
                       exact = FALSE)
    BirdsFrame[i,1] = Pval$p.value
    BirdsFrame[i,2] = Pval$p.value*nrow(relations)
    ############## Bats ##############################################################
    youngBats = relations[i,phenoData2$Group == "Bats" & phenoData2$youngOld == "young"]
    oldBats = relations[i,phenoData2$Group == "Bats" & phenoData2$youngOld == "old"]
    
    Pval = wilcox.test(x = youngBats,
                       y = oldBats, alternative = "two.sided", 
                       paired =FALSE,
                       exact = FALSE)
    BatsFrame[i,1] = Pval$p.value
    BatsFrame[i,2] = Pval$p.value*nrow(relations)
    ############## Rodents ##############################################################
    youngRodents = relations[i,phenoData2$Group == "Rodents" & phenoData2$youngOld == "young"]
    oldRodents = relations[i,phenoData2$Group == "Rodents" & phenoData2$youngOld == "old"]
    
    Pval = wilcox.test(x = youngRodents,
                       y = oldRodents, alternative = "two.sided", 
                       paired =FALSE,
                       exact = FALSE)
    RodentsFrame[i,1] = Pval$p.value
    RodentsFrame[i,2] = Pval$p.value*nrow(relations)
    ############## Primates ##############################################################
    youngPrimates = relations[i,phenoData2$Group == "Primates" & phenoData2$youngOld == "young"]
    oldPrimates = relations[i,phenoData2$Group == "Primates" & phenoData2$youngOld == "old"]
    
    Pval = wilcox.test(x = youngPrimates,
                       y = oldPrimates, alternative = "two.sided", 
                       paired =FALSE,
                       exact = FALSE)
    PrimatesFrame[i,1] = Pval$p.value
    PrimatesFrame[i,2] = Pval$p.value*nrow(relations)
    ############## All ##############################################################
    youngAll = relations[i, phenoData2$youngOld == "young"]
    oldAll = relations[i, phenoData2$youngOld == "old"]
    
    Pval = wilcox.test(x = youngAll,
                       y = oldAll, alternative = "two.sided", 
                       paired =FALSE,
                       exact = TRUE)
    AllFrame[i,1] = Pval$p.value
    AllFrame[i,2] = Pval$p.value*nrow(relations)
    
}

############### Creating all plots ######################################################################



goodFish =      FishFrame[FishFrame[,1] < 0.05 & !is.na(FishFrame[,1]),]
goodWhale =     WhaleFrame[WhaleFrame[,1] < 0.05 & !is.na(WhaleFrame[,1]),]
goodBirds =     BirdsFrame[BirdsFrame[,1] < 0.05 & !is.na(BirdsFrame[,1]),]
goodBats =      BatsFrame[BatsFrame[,1] < 0.05 & !is.na(BatsFrame[,1]),]
goodRodents =   RodentsFrame[RodentsFrame[,1] < 0.05 & !is.na(RodentsFrame[,1]),]
goodPrimates =  PrimatesFrame[PrimatesFrame[,1] < 0.05 & !is.na(PrimatesFrame[,1]),]
goodAll =       AllFrame[AllFrame[,1] < 0.05 & !is.na(AllFrame[,1]),]

length(goodFish)    # 218
length(goodWhale)    # 0
length(goodBirds)   # 110
length(goodBats)   # 0
length(goodRodents)  # 722
length(goodPrimates) # 642

length(goodAll) # 1228

head(goodAll)

# I should check per species!

###############################################################################################

newAge =c(as.character(phenoData2$Organism[phenoData2$Group == "Fish"]),
          as.character(phenoData2$Organism[phenoData2$Group == "Birds"]),
          as.character(phenoData2$Organism[phenoData2$Group == "Whale"]),
          as.character(phenoData2$Organism[phenoData2$Group == "Bats"]),
          as.character(phenoData2$Organism[phenoData2$Group == "Rodents"]),
          as.character(phenoData2$Organism[phenoData2$Group == "Primates"]))

phenoData2 = phenoData2[newAge,]
sortedData = relations[,newAge]

colnames(sortedData) == rownames(phenoData2)

########## ########################## ########################## ##########################
#0. Testing All

makePlots = function (goodData, nameOfSet){
    
    if(class(goodData) == "character"){
      allGood = goodData
    }else{
        allGood = rownames(goodData)  #[1:4000]
    }
    
    currentData = sortedData[allGood,]
    geneOrder =   order(rowSums(sortedData[allGood,]))
    sampleOrder = order(colSums(sortedData[allGood,]))
    
    # Age vs score rank
    plot(phenoData2$Age[sampleOrder],
         main = paste(nameOfSet,"\n", "n=", length(allGood) ),
         xlab = "Gene presence score (Rank)",
         ylab = "Max age in years",
         col = phenoData2$colors[sampleOrder], pch = 16)
    correlation = round( cor(     1:length(phenoData2$Age[sampleOrder]),
                                  phenoData2$Age[sampleOrder]),2)
    legend(x = "topleft", legend = paste("cor(x,y) = ",correlation))
    
   
    ## Age and Gene presence Ranks
    plot( x = colSums(currentData[,sampleOrder]),
          y = phenoData2$Age[sampleOrder],
          main = paste(nameOfSet,"\n", "n=", length(allGood) ),
          xlab = "Gene presence score (Not Ranked)",
          ylab = "Max age in years",
          col = phenoData2$colors[sampleOrder], pch = 16)
      correlation = round( cor(     sort(colSums(sortedData[allGood,])),
                                    phenoData2$Age[sampleOrder]),2)
      legend(x = "topleft", legend = paste("cor(x,y) = ",correlation))
    ##
     

    ## Age and Gene presence Ranks
    plot( x = colSums(currentData[,sampleOrder]),
          y = rank(phenoData2$Age[sampleOrder]),
          main = paste(nameOfSet,"\n", "n=", length(allGood) ),
          xlab = "Gene presence score",
          ylab = "Max age in years (Rank)",
          col = phenoData2$colors[sampleOrder], pch = 16)
      correlation = round( cor(     sort(colSums(sortedData[allGood,])),
                                    rank(phenoData2$Age[sampleOrder])),2)
      legend(x = "topleft", legend = paste("cor(x,y) = ",correlation))
    
    
    detected = ifelse(allGood %in%  ageingData[,1], "forestgreen","red")
    heatmap.2(currentData[geneOrder,sampleOrder],
              main = paste(nameOfSet,"\n", "n=", length(allGood) ),
              labCol = paste(phenoData2$Organism[sampleOrder], phenoData2$Age[sampleOrder],  phenoData2$youngOld[sampleOrder], sep = ' '),
              #labCol = paste(phenoData2$Organism[heaviestGene], phenoData2$localRank[heaviestGene], sep = ' '),
              Rowv = F, Colv = F,
              dendro = "none",
              trace = "none",
              #colCol = phenoData2$AgeColor,
              RowSideColors = detected[geneOrder],
              #colsep = c(9,22,27,32,45), sepcolor = "black",
              ColSideColors = phenoData2$colors[sampleOrder])
    
    heatmap.2(currentData,
              main = paste(nameOfSet,"\n", "n=", length(allGood) ),
              labCol = paste(phenoData2$Organism, phenoData2$Age,  phenoData2$youngOld, sep = ' '),
              #labCol = paste(phenoData2$Organism[heaviestGene], phenoData2$localRank[heaviestGene], sep = ' '),
              Rowv = T, Colv = F,
              dendro = "none",
              trace = "none",
              #colCol = phenoData2$AgeColor,
              RowSideColors = detected,
              colsep = c(9,22,27,32,45), sepcolor = "black",
              ColSideColors = phenoData2$colors)

    
      write.csv(x = goodData, paste("results/B_", nameOfSet,".csv", sep = "")  )
}

########################## ########################## ########################## ##########################
#0. Testing All

pdf("results/allPlots2.pdf")

makePlots(goodAll, "All Data")

makePlots(goodFish, "Fish")
makePlots(goodBirds, "Birds")
# makePlots(goodWhale, "Whale")
# makePlots(goodBats, "Bats")
makePlots(goodRodents, "Rodents")
makePlots(goodPrimates, "Primates")

goodPrimatesandRodents = Reduce(intersect, list(rownames(goodPrimates),rownames(goodRodents)))
makePlots(goodPrimatesandRodents, "Intersect primates and rodents")

goodBirdsRodents = Reduce(intersect, list(rownames(goodBirds),rownames(goodRodents)))
makePlots(goodPrimatesandRodents, "Intersect birds and rodents")


dev.off()

