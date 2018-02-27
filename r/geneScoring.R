# birds = c("Turkey_Vulture",
# "Dalmation_pelican",
# "Barn_owl",
# "Bald_eagle",
# "Cicken",
# "Emperor_penguin",
# "Mallard",
# "Sunbittern",
# "Turkey",
# "Common_cuckoo",
# "Adelie_Penguin",
# "Collared_flycatcher",
# "Annas_Hummingbird")
# 
# lifespan = c(3,3,3,2,2,2,2,1,1,1,1,1,1)
# 
# birdsdf = data.frame(birds,lifespan)
# 
# longevity = factor(c("short","medium","long"))

genesFile = read.csv("~/Documents/david/resources/bats/brandtsNames.tsv", sep = "\t")
geneList = data.frame(genesFile$identifier)
colnames(geneList) = "reference"

allfiles = dir("/home/dylan/Documents/david/scripts_git/out/names/all", full.names = TRUE)
# allfiles = allfiles[44]

orthologyFiles = dir("/home/dylan/data/orthologs/allOrthologyTypes/", full.names = TRUE)
# orthologyFiles = orthologyFiles[44]
test = data.frame(allfiles, orthologyFiles)

allgenes = c()
# idToNameMap = data.frame(character(0), character(0))
# colnames(idToNameMap) = c("id", "genename")
for(file in allfiles){
  print(file)
  namefile = read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  names = data.frame(namefile$V2, namefile$V1)
  colnames(names) = c("id", "genename")
  # idToNameMap = rbind(idToNameMap, names)
  allgenes = c(allgenes, na.omit(namefile[!grepl(".*?_.*?", namefile$V1),1]))
}
allgenes <- data.frame(as.character(unique(toupper(allgenes))))
colnames(allgenes) = c("name")

allgenes <- apply(test, 1, function(x){
  orgCode = strsplit(x[1], "_", fixed = TRUE)[[1]][3]
  print(orgCode)
  namesFile = na.omit(read.delim(x[1], header = FALSE, stringsAsFactors = F))
  orthologyFile = read.delim(x[2], header = FALSE, stringsAsFactors = F)
  orthologyFile[,1] = apply(orthologyFile, 1, splitName)
  allgenes[,orgCode] <- apply(allgenes, 1, function(x){
    allgenes[,orgCode] <- orthologyFile[orthologyFile$V1 == namesFile[namesFile$V1 == x[1],2],2][1]
  })
  return(allgenes)
})
allgenesTable = data.frame(allgenes)
allgenesTable = allgenesTable[, -grep("name.", colnames(allgenesTable))]
allgenesTable <- allgenesTable[-grep("ENS", allgenesTable$name),]

allgenesTable[is.na(test)] <- "missing"
rownames(allgenesTable) <- allgenesTable$name
write.csv(allgenesTable, "~/data/orthologs/allOrganismsOrthologyTable.csv")


orthologyFile = read.delim(as.character(test[,2]), header = FALSE, stringsAsFactors = F)
namesFile = read.delim(as.character(test[,1]), header = FALSE, stringsAsFactors = F)
orthologyFile[,1] = apply(orthologyFile, 1, splitName)
orgCode = "min"
allgenes[,orgCode] = rep("missing", nrow(allgenes))
allgenes <- apply(orthologyFile, 1, function(x){
  #print(allgenes[allgenes$name == namesFile[namesFile[,2] == x[1],1],])
  allgenes[,orgCode] <- orthologyFile[orthologyFile$V1 == namesFile[namesFile$V1 == x[1],2],2][1]
})
allgenes[allgenes$name == namesFile[namesFile$V2 == x[1],1],2] <- x[2]

allgenes[,orgCode] <- apply(allgenes, 1, function(x){
  
  allgenes[,orgCode] <- orthologyFile[orthologyFile$V1 == namesFile[namesFile$V1 == x[1],2],2][1]
})


df <- data.frame(1:3,5:7)

apply(df, 1, function(x){
  df[df$X1.3 == 2,2] <- 100
})

df[df$X1.3 == 2,2] <- 100





presenceCheck = function(x){
  if (length(ortData[ortData$reference == x[1],1]) == 0){
    return("missing")
  }
  else{
    if(length(ortData[ortData$reference == x[1],1]) == 1 && ortData[ortData$reference == x[1],2] == "one2one"){
      return("one2one")
    }
    else{
      if(length(ortData[ortData$reference == x[1],1]) == 1 && ortData[ortData$reference == x[1],2] == "many2one"){
        return("deletion")
      }
      else{
        return("duplication")
      }
    }
  }
}
# orgCode = strsplit(first, "_", fixed = TRUE)[[1]][2]
# geneList[,orgCode] = apply(geneList, 1, presenceCheck)
# 
# geneList[geneList$reference == "Cau_R000992",]

splitName = function(x){
  splitname = strsplit(x, split="|", fixed = T)[[1]][2]
  
  if(is.na(splitname)){
    return(as.character(x[1]))
  }
  else{
    return(as.character(splitname))
  }
}


findGeneName = function(id){
  splitid = strsplit(id, split="|", fixed = T)[[1]][2]
  if(is.na(splitid)){
    return(idToNameMap[idToNameMap$id == id,2])
  }
  else{
    return(as.character(idToNameMap[idToNameMap$id == splitid,2]))
  }
}

idstonames = apply(ortData, 1, findGeneName)

for (file in orthologyFiles){
  orgCode = strsplit(file, "_", fixed = TRUE)[[1]][3]
  print(orgCode)
  ortData = read.csv(file, sep="\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(ortData) = c("reference", "type", "ortholog")
  ortData$reference = lapply(ortData$reference, splitName)
  geneList[,orgCode] = apply(geneList, 1, presenceCheck)
}

write.csv(geneList, "~/data/orthologs/bats/genescores_duplication.csv")
bats = read.csv("~/data/orthologs/bats/genescores_duplication.csv")
row.names(bats) = bats$X
bats$X = NULL
#################### read dnds file(s) ########################
test = data.frame(geneList$reference, geneList$omy)
omy_dnds = read.delim("~/Documents/david/scripts_git/out/dnds/ungapped/fish/loc_omy_dnds.tsv")
omyFiltered = omy_dnds$ref[omy_dnds$dN.dS > 0.0001 & omy_dnds$dN.dS < 10]

factor(x = geneList$efu, labels = )

changecol = function(colName){
  ifelse(colName == 0, "missing", ifelse(colName == 1, "deletion", ifelse(colName == 2, "one2one", "duplication")))
}
# bats$min = as.factor(changecol(bats$min))
# write.csv(bats, "~/data/batsOrthologyTable.csv")
# 
# geneList$myd = as.factor(changecol(geneList$myd))

geneSubset = subset(geneList, !(min == myd & myd == pale & pale == efu))



clean_interaction_set = subset(data.frame(clean_interaction_set), !(efu == min & min == myd & myd == pale))
pivoted_clean_interaction_set = data.frame(t(clean_interaction_set))
lifespan_orthology_regression = lm(lifespan_category ~ 0 + . , data = pivoted_clean_interaction_set)