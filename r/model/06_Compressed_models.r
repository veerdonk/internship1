Pcut = 0.0001
#library(KEGGREST)

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

#Align names
data_matrix = data_matrix[,rownames(phenoData2)]
## I hereby see that the phenoData matches data_matrix
colnames(data_matrix ) == rownames(phenoData2)


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
ageingData = read.table("data/ageing_genes.csv", sep = ",", head = T, stringsAsFactors  = T)


#########################################################################
## Average age calculations
#########################################################################

FishFrame = relations[,1:4]
colnames(FishFrame) = c("none2one", "none2many","none2more","one2many")

WhaleFrame = relations[,1:4]
colnames(WhaleFrame) = c("none2one", "none2many","none2more","one2many")

BirdsFrame = relations[,1:4]
colnames(BirdsFrame) = c("none2one", "none2many","none2more","one2many")

BatsFrame = relations[,1:4]
colnames(BatsFrame) = c("none2one", "none2many","none2more","one2many")

RodentsFrame = relations[,1:4]
colnames(RodentsFrame) = c("none2one", "none2many","none2more","one2many")

PrimatesFrame = relations[,1:4]
colnames(PrimatesFrame) = c("none2one", "none2many","none2more","one2many")

AllFrame = relations[,1:4]
colnames(AllFrame) = c("none2one", "none2many","none2more","one2many")

NotfishFrame = relations[,1:4]
colnames(NotfishFrame) = c("none2one", "none2many","none2more","one2many")

NotfishbirdsFrame = relations[,1:4]
colnames(NotfishbirdsFrame) = c("none2one", "none2many","none2more","one2many")



for (i in 1:nrow(FishFrame)){
    print(i)
    
    ############## Fish ##############################################################
    # Missing
    missingFish =phenoData2$Age[relations[i,] == 1 & phenoData2$Group == "Fish" ]
    # One
    oneFish = phenoData2$Age[relations[i,] == 2 & phenoData2$Group == "Fish" ]
    # Many
    manyFish = phenoData2$Age[relations[i,] == 3 & phenoData2$Group == "Fish" ]
    
    # none to one
    if ( length(missingFish) != 0 & length(oneFish) != 0 ){
      none2one = wilcox.test(x = missingFish,
                             y = oneFish, alternative = "two.sided", 
                             paired =FALSE,
                             exact = FALSE)
      none2one = none2one$p.value
    }else{ none2one = 1}
    FishFrame[i,1] = none2one
    
    # none to many
    if ( length(missingFish) != 0 & length(manyFish) != 0 ){
      none2many = wilcox.test(x = missingFish,
                              y = manyFish, alternative = "two.sided", 
                              paired =FALSE,
                              exact = FALSE)
      none2many = none2many$p.value
    }else{ none2many = 1}
    FishFrame[i,2]  = none2many
    
    # none to more
    if ( length(missingFish) != 0 & length(oneFish) != 0 & length(manyFish) != 0){
      none2more = wilcox.test(x = missingFish,
                              y = c(oneFish,manyFish), alternative = "two.sided", 
                              paired =FALSE,
                              exact = FALSE)
      none2more = none2more$p.value
    }else{ none2more = 1}
    FishFrame[i,3]  = none2more
    
    # one to many
    if ( length(oneFish) != 0 & length(manyFish) != 0){
      one2many = wilcox.test(x = oneFish,
                            y = manyFish, alternative = "two.sided", 
                            paired =FALSE,
                            exact = FALSE)
      one2many = one2many$p.value
    }else{ one2many = 1}
    FishFrame[i,4] = one2many
  
    ############## Whale ##############################################################
    # Missing
    missingWhale = phenoData2$localRank[relations[i,] == 1 & phenoData2$Group == "Whale" ]
    # One
    oneWhale = phenoData2$localRank[relations[i,] == 2 & phenoData2$Group == "Whale" ]
    # Many
    manyWhale = phenoData2$localRank[relations[i,] == 3 & phenoData2$Group == "Whale" ]
    
    # none to one
    if ( length(missingWhale) != 0 & length(oneWhale) != 0 & length(manyWhale) == 0 ){
      customTest = (max(missingWhale) < min(oneWhale) |
                      min(missingWhale) > max(oneWhale))
      none2one = ifelse(customTest, 0.04, 1)
      
    }else{ none2one = 1}
    WhaleFrame[i,1] = none2one
    
    # none to many
    if ( length(missingWhale) != 0 & length(oneWhale) == 0 & length(manyWhale) != 0 ){
      customTest = (max(missingWhale) < min(manyWhale) |
                      min(missingWhale) > max(manyWhale))
      none2many = ifelse(customTest, 0.04, 1)
    }else{ none2many = 1}
    WhaleFrame[i,2]  = none2many
    
    # none to more
    if ( length(missingWhale) != 0 & length(oneWhale) != 0 & length(manyWhale) != 0){
      customTest = (max(missingWhale) < min(oneWhale) & min(oneWhale) < min( manyWhale) |
                      min(missingWhale) > max(oneWhale) & max(oneWhale) > max(manyWhale))
      none2more = ifelse(customTest, 0.04, 1)
    }else{ none2more = 1}
    WhaleFrame[i,3]  = none2more
    
    # one to many
    if ( length(missingWhale) == 0 & length(oneWhale) != 0 & length(manyWhale) != 0){
      customTest = (max(oneWhale) < min(manyWhale) |
                      min(oneWhale) > max(manyWhale))
      one2many = ifelse(customTest, 0.04, 1)
    }else{ one2many = 1}
    WhaleFrame[i,4] = one2many
    
    
    
    ############## Birds ##############################################################
    # Missing
    missingBirds =phenoData2$Age[relations[i,] == 1 & phenoData2$Group == "Birds" ]
    # One
    oneBirds = phenoData2$Age[relations[i,] == 2 & phenoData2$Group == "Birds" ]
    # Many
    manyBirds = phenoData2$Age[relations[i,] == 3 & phenoData2$Group == "Birds" ]
    
    # none to one
    if ( length(missingBirds) != 0 & length(oneBirds) != 0 ){
      none2one = wilcox.test(x = missingBirds,
                             y = oneBirds, alternative = "two.sided", 
                             paired =FALSE,
                             exact = FALSE)
      none2one = none2one$p.value
    }else{ none2one = 1}
    BirdsFrame[i,1] = none2one
    
    # none to many
    if ( length(missingBirds) != 0 & length(manyBirds) != 0 ){
      none2many = wilcox.test(x = missingBirds,
                              y = manyBirds, alternative = "two.sided", 
                              paired =FALSE,
                              exact = FALSE)
      none2many = none2many$p.value
    }else{ none2many = 1}
    BirdsFrame[i,2]  = none2many
    
    # none to more
    if ( length(missingBirds) != 0 & length(oneBirds) != 0 & length(manyBirds) != 0){
      none2more = wilcox.test(x = missingBirds,
                              y = c(oneBirds,manyBirds), alternative = "two.sided", 
                              paired =FALSE,
                              exact = FALSE)
      none2more = none2more$p.value
    }else{ none2more = 1}
    BirdsFrame[i,3]  = none2more
    
    # one to many
    if ( length(oneBirds) != 0 & length(manyBirds) != 0){
      one2many = wilcox.test(x = oneBirds,
                             y = manyBirds, alternative = "two.sided", 
                             paired =FALSE,
                             exact = FALSE)
      one2many = one2many$p.value
    }else{ one2many = 1}
    BirdsFrame[i,4] = one2many
    
    ############## Bats ##############################################################
    # Missing
    missingBats = phenoData2$localRank[relations[i,] == 1 & phenoData2$Group == "Bats" ]
    # One
    oneBats = phenoData2$localRank[relations[i,] == 2 & phenoData2$Group == "Bats" ]
    # Many
    manyBats = phenoData2$localRank[relations[i,] == 3 & phenoData2$Group == "Bats" ]
    
    # none to one
    if ( length(missingBats) != 0 & length(oneBats) != 0 & length(manyBats) == 0 ){
      customTest = (max(missingBats) < min(oneBats) |
                      min(missingBats) > max(oneBats))
      none2one = ifelse(customTest, 0.04, 1)
      
    }else{ none2one = 1}
    BatsFrame[i,1] = none2one
    
    # none to many
    if ( length(missingBats) != 0 & length(oneBats) == 0 & length(manyBats) != 0 ){
      customTest = (max(missingBats) < min(manyBats) |
                      min(missingBats) > max(manyBats))
      none2many = ifelse(customTest, 0.04, 1)
    }else{ none2many = 1}
    BatsFrame[i,2]  = none2many
    
    # none to more
    if ( length(missingBats) != 0 & length(oneBats) != 0 & length(manyBats) != 0){
      customTest = (max(missingBats) < min(oneBats) & min(oneBats) < min( manyBats) |
                      min(missingBats) > max(oneBats) & max(oneBats) > max(manyBats))
      none2more = ifelse(customTest, 0.04, 1)
    }else{ none2more = 1}
    BatsFrame[i,3]  = none2more
    
    # one to many
    if ( length(missingBats) == 0 & length(oneBats) != 0 & length(manyBats) != 0){
      customTest = (max(oneBats) < min(manyBats) |
                      min(oneBats) > max(manyBats))
      one2many = ifelse(customTest, 0.04, 1)
    }else{ one2many = 1}
    BatsFrame[i,4] = one2many
    
    
    
    ############## Rodents ##############################################################
    # Missing
    missingRodents =phenoData2$Age[relations[i,] == 1 & phenoData2$Group == "Rodents" ]
    # One
    oneRodents = phenoData2$Age[relations[i,] == 2 & phenoData2$Group == "Rodents" ]
    # Many
    manyRodents = phenoData2$Age[relations[i,] == 3 & phenoData2$Group == "Rodents" ]
    
    # none to one
    if ( length(missingRodents) != 0 & length(oneRodents) != 0 & length(manyRodents) == 0 ){
      none2one = wilcox.test(x = missingRodents,
                             y = oneRodents, alternative = "two.sided", 
                             paired =FALSE,
                             exact = FALSE)
      none2one = none2one$p.value
    }else{ none2one = 1}
    RodentsFrame[i,1] = none2one
    
    # none to many
    if ( length(missingRodents) != 0 & length(oneRodents) == 0 &  length(manyRodents) != 0 ){
      none2many = wilcox.test(x = missingRodents,
                              y = manyRodents, alternative = "two.sided", 
                              paired =FALSE,
                              exact = FALSE)
      none2many = none2many$p.value
    }else{ none2many = 1}
    RodentsFrame[i,2]  = none2many
    
    # none to more
    if ( length(missingRodents) != 0 & length(oneRodents) != 0 & length(manyRodents) != 0){
      none2more = wilcox.test(x = missingRodents,
                              y = c(oneRodents,manyRodents), alternative = "two.sided", 
                              paired =FALSE,
                              exact = FALSE)
      none2more = none2more$p.value
    }else{ none2more = 1}
    RodentsFrame[i,3]  = none2more
    
    # one to many
    if ( length(missingRodents) == 0 & length(oneRodents) != 0 & length(manyRodents) != 0){
      one2many = wilcox.test(x = oneRodents,
                             y = manyRodents, alternative = "two.sided", 
                             paired =FALSE,
                             exact = FALSE)
      one2many = one2many$p.value
    }else{ one2many = 1}
    RodentsFrame[i,4] = one2many
    
    ############## Primates ##############################################################
    # Missing
    missingPrimates = phenoData2$localRank[relations[i,] == 1 & phenoData2$Group == "Primates" ]
    # One
    onePrimates = phenoData2$localRank[relations[i,] == 2 & phenoData2$Group == "Primates" ]
    # Many
    manyPrimates = phenoData2$localRank[relations[i,] == 3 & phenoData2$Group == "Primates" ]
    
    # none to one
    if ( length(missingPrimates) != 0 & length(onePrimates) != 0 & length(manyPrimates) == 0 ){
      customTest =   (max(missingPrimates) < min(onePrimates) |
                      min(missingPrimates) > max(onePrimates))
      none2one = ifelse(customTest, 0.04, 1)
      
    }else{ none2one = 1}
    PrimatesFrame[i,1] = none2one
    
    # none to many
    if ( length(missingPrimates) != 0 & length(onePrimates) == 0 & length(manyPrimates) != 0 ){
      customTest = (max(missingPrimates) < min(manyPrimates) |
                      min(missingPrimates) > max(manyPrimates))
      none2many = ifelse(customTest, 0.04, 1)
    }else{ none2many = 1}
    PrimatesFrame[i,2]  = none2many
    
    # none to more
    if ( length(missingPrimates) != 0 & length(onePrimates) != 0 & length(manyPrimates) != 0){
      customTest =   (max(missingPrimates) < min(onePrimates) & max(onePrimates) < min(manyPrimates) |
                      min(missingPrimates) > max(onePrimates) & min(onePrimates) > max(manyPrimates))
      none2more = ifelse(customTest, 0.04, 1)
    }else{ none2more = 1}
    PrimatesFrame[i,3]  = none2more
    
    # one to many
    if ( length(missingPrimates) == 0 & length(onePrimates) != 0 & length(manyPrimates) != 0){
      customTest = (max(onePrimates) < min(manyPrimates) |
                      min(onePrimates) > max(manyPrimates))
      one2many = ifelse(customTest, 0.04, 1)
    }else{ one2many = 1}
    PrimatesFrame[i,4] = one2many
    
    
    
    ############## All ##############################################################
    # Missing
    missingAll =phenoData2$Age[relations[i,] == 1]
    # One
    oneAll = phenoData2$Age[relations[i,] == 2]
    # Many
    manyAll = phenoData2$Age[relations[i,] == 3]
    
    # none to one
    if ( length(missingAll) != 0 & length(oneAll) != 0 ){
      none2one = wilcox.test(x = missingAll,
                             y = oneAll, alternative = "two.sided", 
                             paired =FALSE,
                             exact = FALSE)
      none2one = none2one$p.value
    }else{ none2one = 1}
    AllFrame[i,1] = none2one
    
    # none to many
    if ( length(missingAll) != 0 & length(manyAll) != 0 ){
      none2many = wilcox.test(x = missingAll,
                              y = manyAll, alternative = "two.sided", 
                              paired =FALSE,
                              exact = FALSE)
      none2many = none2many$p.value
    }else{ none2many = 1}
    AllFrame[i,2]  = none2many
    
    # none to more
    if ( length(missingAll) != 0 & length(oneAll) != 0 & length(manyAll) != 0){
      none2more = wilcox.test(x = missingAll,
                              y = c(oneAll,manyAll), alternative = "two.sided", 
                              paired =FALSE,
                              exact = FALSE)
      none2more = none2more$p.value
    }else{ none2more = 1}
    AllFrame[i,3]  = none2more
    
    # one to many
    if ( length(oneAll) != 0 & length(manyAll) != 0){
      one2many = wilcox.test(x = oneAll,
                             y = manyAll, alternative = "two.sided", 
                             paired =FALSE,
                             exact = FALSE)
      one2many = one2many$p.value
    }else{ one2many = 1}
    AllFrame[i,4] = one2many
    
    
    ############## Notfish ##############################################################
    # Missing
    missingNotfish =phenoData2$Age[relations[i,] == 1 & phenoData2$Group != "Fish" ]
    # One
    oneNotfish = phenoData2$Age[relations[i,] == 2 & phenoData2$Group != "Fish" ]
    # Many
    manyNotfish = phenoData2$Age[relations[i,] == 3 & phenoData2$Group != "Fish" ]
    
    # none to one
    if ( length(missingNotfish) != 0 & length(oneNotfish) != 0 ){
      none2one = wilcox.test(x = missingNotfish,
                             y = oneNotfish, alternative = "two.sided", 
                             paired =FALSE,
                             exact = FALSE)
      none2one = none2one$p.value
    }else{ none2one = 1}
    NotfishFrame[i,1] = none2one
    
    # none to many
    if ( length(missingNotfish) != 0 & length(manyNotfish) != 0 ){
      none2many = wilcox.test(x = missingNotfish,
                              y = manyNotfish, alternative = "two.sided", 
                              paired =FALSE,
                              exact = FALSE)
      none2many = none2many$p.value
    }else{ none2many = 1}
    NotfishFrame[i,2]  = none2many
    
    # none to more
    if ( length(missingNotfish) != 0 & length(oneNotfish) != 0 & length(manyNotfish) != 0){
      none2more = wilcox.test(x = missingNotfish,
                              y = c(oneNotfish,manyNotfish), alternative = "two.sided", 
                              paired =FALSE,
                              exact = FALSE)
      none2more = none2more$p.value
    }else{ none2more = 1}
    NotfishFrame[i,3]  = none2more
    
    # one to many
    if ( length(oneNotfish) != 0 & length(manyNotfish) != 0){
      one2many = wilcox.test(x = oneNotfish,
                             y = manyNotfish, alternative = "two.sided", 
                             paired =FALSE,
                             exact = FALSE)
      one2many = one2many$p.value
    }else{ one2many = 1}
    NotfishFrame[i,4] = one2many
    
    ############## Notfish AND notBirds ##############################################################
    # Missing
    missingNotfishbirds =phenoData2$Age[relations[i,] == 1 & phenoData2$Group != "Fish" & phenoData2$Group != "Birds" ]
    # One
    oneNotfishbirds = phenoData2$Age[relations[i,] == 2 & phenoData2$Group != "Fish" & phenoData2$Group != "Birds" ]
    # Many
    manyNotfishbirds = phenoData2$Age[relations[i,] == 3 & phenoData2$Group != "Fish" & phenoData2$Group != "Birds" ]
    
    # none to one
    if ( length(missingNotfishbirds) != 0 & length(oneNotfishbirds) != 0 ){
      none2one = wilcox.test(x = missingNotfishbirds,
                             y = oneNotfishbirds, alternative = "two.sided", 
                             paired =FALSE,
                             exact = FALSE)
      none2one = none2one$p.value
    }else{ none2one = 1}
    NotfishbirdsFrame[i,1] = none2one
    
    # none to many
    if ( length(missingNotfishbirds) != 0 & length(manyNotfishbirds) != 0 ){
      none2many = wilcox.test(x = missingNotfishbirds,
                              y = manyNotfishbirds, alternative = "two.sided", 
                              paired =FALSE,
                              exact = FALSE)
      none2many = none2many$p.value
    }else{ none2many = 1}
    NotfishbirdsFrame[i,2]  = none2many
    
    # none to more
    if ( length(missingNotfishbirds) != 0 & length(oneNotfishbirds) != 0 & length(manyNotfishbirds) != 0){
      none2more = wilcox.test(x = missingNotfishbirds,
                              y = c(oneNotfishbirds,manyNotfishbirds), alternative = "two.sided", 
                              paired =FALSE,
                              exact = FALSE)
      none2more = none2more$p.value
    }else{ none2more = 1}
    NotfishbirdsFrame[i,3]  = none2more
    
    # one to many
    if ( length(oneNotfishbirds) != 0 & length(manyNotfishbirds) != 0){
      one2many = wilcox.test(x = oneNotfishbirds,
                             y = manyNotfishbirds, alternative = "two.sided", 
                             paired =FALSE,
                             exact = FALSE)
      one2many = one2many$p.value
    }else{ one2many = 1}
    NotfishFrame[i,4] = one2many
    
    
    
}


############### Creating all plots ######################################################################

goodFish =     na.omit(FishFrame[apply(FishFrame, 1, FUN = min) < 0.05,])
goodWhale =    na.omit(WhaleFrame[apply(WhaleFrame, 1, FUN = min) < 0.05,])
goodBirds =    na.omit(BirdsFrame[apply(BirdsFrame, 1, FUN = min) < 0.05,])
goodBats =     na.omit(BatsFrame[apply(BatsFrame, 1, FUN = min) < 0.05,])
goodRodents =  na.omit(RodentsFrame[apply(RodentsFrame, 1, FUN = min) < 0.05,])
goodPrimates = na.omit(PrimatesFrame[apply(PrimatesFrame, 1, FUN = min) < 0.05,])

NotfishFrame[,3] = ifelse(NotfishFrame[,3] < Pcut & NotfishFrame[,4] > Pcut, 1, NotfishFrame[,3] )
goodNotfish = na.omit(NotfishFrame[apply(NotfishFrame, 1, FUN = min) < Pcut,])

NotfishbirdsFrame[,3] = ifelse(NotfishbirdsFrame[,3] < Pcut & NotfishbirdsFrame[,4] > Pcut, 1, NotfishbirdsFrame[,3] )
goodNotfishbirds = na.omit(NotfishbirdsFrame[apply(NotfishbirdsFrame, 1, FUN = min) < Pcut,])

AllFrame[,3] = ifelse(AllFrame[,3] < Pcut & AllFrame[,4] > Pcut, 1, AllFrame[,3] )
goodAll = na.omit(AllFrame[apply(AllFrame, 1, FUN = min) < Pcut,])



length(rownames(goodFish))     # 88
length(rownames(goodWhale))    # 1362
length(rownames(goodBirds))    # 192
length(rownames(goodBats))    # 1378
length(rownames(goodRodents))  # 3093
length(rownames(goodPrimates)) # 1042

length(rownames(goodAll)) # 1027
length(rownames(goodNotfish)) #701

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
              labCol = paste(phenoData2$Organism[sampleOrder], phenoData2$Age[sampleOrder], sep = ' '),
              #labCol = paste(phenoData2$Organism[heaviestGene], phenoData2$localRank[heaviestGene], sep = ' '),
              Rowv = T, Colv = F,
              dendro = "none",
              trace = "none",
              #colCol = phenoData2$AgeColor,
              RowSideColors = detected[geneOrder],
              #colsep = c(9,22,27,32,45), sepcolor = "black",
              ColSideColors = phenoData2$colors[sampleOrder])
    
    heatmap.2(currentData,
              main = paste(nameOfSet,"\n", "n=", length(allGood) ),
              labCol = paste(phenoData2$Organism, phenoData2$Age, sep = ' '),
              #labCol = paste(phenoData2$Organism[heaviestGene], phenoData2$localRank[heaviestGene], sep = ' '),
              Rowv = T, Colv = F,
              dendro = "none",
              trace = "none",
              #colCol = phenoData2$AgeColor,
              RowSideColors = detected,
              colsep = c(9,22,27,32,45), sepcolor = "black",
              ColSideColors = phenoData2$colors)

    
      write.csv(x = goodData, paste("results/", nameOfSet,".csv", sep = "")  )
}

########################## ########################## ########################## ##########################
#0. Testing All

pdf("results/allPlots.pdf")

makePlots(goodAll, "All Data")

makePlots(goodNotfish, "Not Fish")
makePlots(goodNotfishbirds, "Not Fish nor Birds")

makePlots(goodFish, "Fish")
makePlots(goodBirds, "Birds")
makePlots(goodWhale, "Whale")
makePlots(goodBats, "Bats")
makePlots(goodRodents, "Rodents")
makePlots(goodPrimates, "Primates")

goodPrimatesandRodents = Reduce(intersect, list(rownames(goodPrimates),rownames(goodRodents)))
makePlots(goodPrimatesandRodents, "Intersect primates and rodents")

goodBirdsRodents = Reduce(intersect, list(rownames(goodBirds),rownames(goodRodents)))
makePlots(goodPrimatesandRodents, "Intersect birds and rodents")


dev.off()

