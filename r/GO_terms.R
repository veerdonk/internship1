##########################################################################
## 4. GO-terms
##########################################################################
#install.packages("gProfileR")
library("gProfileR")
###########################################################################
# # # Gprofile background
# background = rownames(dat)
# localGprofile   =   gprofiler(background, organism = "mmusculus", ordered_query = F,
#                               significant = F, exclude_iea = F, underrep = F, evcodes = F,
#                               region_query = F, max_p_value = 1, min_set_size = 0, max_set_size = 0,
#                               min_isect_size = 0, correction_method = "analytical",
#                               hier_filtering = "none", domain_size = "annotated", custom_bg = "",
#                               numeric_ns = "", png_fn = NULL, include_graph = F, src_filter = c("GO:BP", "GO:MF", "GO:CC"))
# 
# # src_filter = c("GO:BP", "GO:MF", "GO:CC")
# # src_filter = NULL
# #localGprofile = localGprofile[,1:13]
# write.table(localGprofile, file =  "data/GObackground2.csv", sep = "\t", row.names = FALSE)
###########################################################################


#query = background#YoungWT_OldWT


# For background use: significant = F, hier_filtering = "none"
goProfileThis = function(query, outname){
  localGprofile = gprofiler(query, organism = "hsapiens", ordered_query = F,
                                significant = F, exclude_iea = F, underrep = F, evcodes = F,
                                region_query = F, max_p_value = 1, min_set_size = 0, max_set_size = 0,
                                min_isect_size = 0, correction_method = "analytical",
                                hier_filtering = "none", #("none","moderate, strong")
                                domain_size = "annotated", custom_bg = "",
                                numeric_ns = "", png_fn = NULL, include_graph = F, src_filter = c("GO:BP", "GO:MF", "GO:CC"))     # src_filter = c("GO:BP", "GO:MF", "GO:CC")
  ## (GO:BP, GO:MF, GO:CC to select a particular GO branch), KEGG, REAC, TF, MI, CORUM, HP, HPA, OMIM.
  ## src_filter = NULL
  localGprofile = localGprofile[,1:13]
  rownames(localGprofile) =  localGprofile$term.name

  outFolder = "~/data/stopCodons/goTerms/gprofiled/"
  outFile = paste(outname,".csv",sep = "")
  write.csv(localGprofile, file =  paste(outFolder,outFile, sep = ""), row.names = FALSE)
  return(localGprofile)
}

plotGO = function(inSet, plotLabel){
  cutoff = 1.0
  a = reference[rownames(inSet),6]/reference[rownames(inSet),5] # GO terms for the reference
  b = inSet[,6]/inSet[,5]                                       # GO terms for my querie     

  Interest = rownames(inSet[which(b/a > cutoff),])
  localCol = ifelse (b/a > cutoff, "red","black")
  range = c(min(c(a,b), na.rm = T), max(c(a,b), na.rm = T))

  
  if (length(Interest) >0){
    print(paste( plotLabel,length(Interest)))
    myBar = barplot(( inSet[Interest,6]/inSet[Interest,5]) /
                      (reference[Interest,6]/reference[Interest,5]), las = 2 ,
                    horiz = TRUE,
                    col = c("lightgrey", "darkgrey"))#,myColors))
    axis(2, at=myBar[,1], 
         labels=inSet[Interest,12], 
         las = 2,
         cex.axis = ifelse(length(Interest)>10,0.6,0.8))
    
    text(inSet[Interest,6], y = myBar[,1], x = 0.4, cex = 0.5)
    
    title(main = paste(plotLabel,"\n","overrepresented"),
          cex.main = 0.75)
  }
  #abline(v = 1, lty = 2)
  #abline(v = 1.5, lty = 2)
  return(Interest)
}

# Load the background
# reference = read.csv("data/GObackground.csv", sep = ",")
# rownames(reference) = reference$term.name
# #pdf("results/GO_terms.pdf")
# 
# 
# 
# pdf("GO_Terms.pdf")
# 
# par(mar=c(5.1,18.1,4.1,2.1), xpd = TRUE,mfrow=c(1,1)) #  bottom, left, top and right
# 
# 
# # 1. Consequences of ageing without mutation.
# GO  = goProfileThis(YoungWT_OldWT, "YoungWT_OldWT")
# plotGO(GO, "YoungWT_OldWT")

####################################################################################
################################# running gprofiler ################################
### Background
# allGenes = read.csv("~/data/stopCodons/goTerms/allGenes.csv", stringsAsFactors = F)
# background = allGenes$x
# GObackground = goProfileThis(background, "background")
# reference = GObackground
# write.csv(reference, "~/data/stopCodons/goTerms/gprofiled/background.csv")

### conserved stop = TAA
# conservedTaa = read.csv("~/data/stopCodons/goTerms/conserved_taa.csv", stringsAsFactors = F)
# taaquery = conservedTaa$x
# tqnf = goProfileThis(taaquery, "taaqueryNoFilter")
# goProfileThis(taaquery, "taaqueryModerateFilter")
# goProfileThis(taaquery, "taaqueryStrongFilter")
# 
# ### conserved stop = TGA
# conservedTga = read.csv("~/data/stopCodons/goTerms/conservedTga.csv", stringsAsFactors = F)
# tgaquery = conservedTga$x
# goProfileThis(tgaquery, "tgaqueryNoFilter")
# goProfileThis(tgaquery, "tgaqueryModerateFilter")
# goProfileThis(tgaquery, "tgaqueryStrongFilter")
# 
# ### conserved stop = TAG
# conservedTag = read.csv("~/data/stopCodons/goTerms/conservedTag.csv", stringsAsFactors = F)
# tagquery = conservedTag$x
# goProfileThis(tagquery, "tagqueryNoFilter")
# goProfileThis(tagquery, "tagqueryModerateFilter")
# goProfileThis(tagquery, "tagqueryStrongFilter")
# 
# ### non conseved stops
# nonConserved = read.csv("~/data/stopCodons/goTerms/nonConserved.csv", stringsAsFactors = F)
# nonConservedQuery = nonConserved$x
# goProfileThis(nonConservedQuery, "nonConservedQueryNoFilter")
# goProfileThis(nonConservedQuery, "nonConservedQueryModerateFilter")
# goProfileThis(nonConservedQuery, "nonConservedQueryStrongFilter")

################################# load previously run data ################################
### loading TAA with 3 levels of hierarchical filtering
taaGOqueryNoFilter = read.csv("~/data/stopCodons/goTerms/gprofiled/taaqueryNoFilter.csv", stringsAsFactors = F)
taaGOqueryModerateFilter = read.csv("~/data/stopCodons/goTerms/gprofiled/taaqueryModerateFilter.csv", stringsAsFactors = F)
taaGOqueryStrongFilter = read.csv("~/data/stopCodons/goTerms/gprofiled/taaqueryStrongFilter.csv", stringsAsFactors = F)
rownames(taaGOqueryNoFilter) = taaGOqueryNoFilter$term.name
rownames(taaGOqueryModerateFilter) = taaGOqueryModerateFilter$term.name
rownames(taaGOqueryStrongFilter) = taaGOqueryStrongFilter$term.name

### loading TGA with 3 levels of hierarchical filtering
tgaGOqueryNoFilter = read.csv("~/data/stopCodons/goTerms/gprofiled/tgaqueryNoFilter.csv", stringsAsFactors = F)
tgaGOqueryModerateFilter = read.csv("~/data/stopCodons/goTerms/gprofiled/tgaqueryModerateFilter.csv", stringsAsFactors = F)
tgaGOqueryStrongFilter = read.csv("~/data/stopCodons/goTerms/gprofiled/tgaqueryStrongFilter.csv", stringsAsFactors = F)
rownames(tgaGOqueryNoFilter) = tgaGOqueryNoFilter$term.name
rownames(tgaGOqueryModerateFilter) = tgaGOqueryModerateFilter$term.name
rownames(tgaGOqueryStrongFilter) = tgaGOqueryStrongFilter$term.name

### loading TAG with 3 levels of hierarchical filtering
tagGOqueryNoFilter = read.csv("~/data/stopCodons/goTerms/gprofiled/tagqueryNoFilter.csv", stringsAsFactors = F)
tagGOqueryModerateFilter =  read.csv("~/data/stopCodons/goTerms/gprofiled/tagqueryModerateFilter.csv", stringsAsFactors = F)
tagGOqueryStrongFilter = read.csv("~/data/stopCodons/goTerms/gprofiled/tagqueryStrongFilter.csv", stringsAsFactors = F)
rownames(tagGOqueryNoFilter) = tagGOqueryNoFilter$term.name
rownames(tagGOqueryModerateFilter) = tagGOqueryModerateFilter$term.name
rownames(tagGOqueryStrongFilter) = tagGOqueryStrongFilter$term.name

### loading non conserved stop codons
nonConservedGOqueryNoFilter = read.csv("~/data/stopCodons/goTerms/gprofiled/nonConservedQueryNoFilter.csv", stringsAsFactors = F)
nonConservedGOqueryModerateFilter = read.csv("~/data/stopCodons/goTerms/gprofiled/nonConservedQueryModerateFilter.csv", stringsAsFactors = F)
nonConservedGOqueryStrongFilter = read.csv("~/data/stopCodons/goTerms/gprofiled/nonConservedQueryStrongFilter.csv", stringsAsFactors = F)
rownames(nonConservedGOqueryNoFilter) = nonConservedGOqueryNoFilter$term.name
rownames(nonConservedGOqueryModerateFilter) = nonConservedGOqueryModerateFilter$term.name
rownames(nonConservedGOqueryStrongFilter) = nonConservedGOqueryStrongFilter$term.name

### loading background -no filtering -no significance
reference = read.csv("~/data/stopCodons/goTerms/gprofiled/background.csv")
rownames(reference) = reference$term.name


### plotting bar plots
pdf("GOBarPlots.pdf")
## set margins for plot so go terms can be read
par(mar=c(5,16,4,1)+.1)

## plot no filter
plotGO(taaGOqueryNoFilter, "conserved TAA stopcodons\n-no filtering -cutoff 1.8")
plotGO(tgaGOqueryNoFilter, "conserved TGA stopcodons\n-no filtering -cutoff 1.8")
plotGO(tagGOqueryNoFilter, "conserved TAG stopcodons\n-no filtering -cutoff 1.8")
plotGO(nonConservedGOqueryNoFilter, "non conserved stopcodons\n-no filtering -cutoff 1.0")

## plot moderate filter
plotGO(taaGOqueryModerateFilter, "conserved TAA stopcodons\n-moderate filtering -cutoff 1.0")
plotGO(tgaGOqueryModerateFilter, "conserved TGA stopcodons\n-moderate filtering -cutoff 1.0")
plotGO(tagGOqueryModerateFilter, "conserved TAG stopcodons\n-moderate filtering -cutoff 1.0")
plotGO(nonConservedGOqueryModerateFilter, "non conserved stopcodons\n-no filtering -cutoff 1.0")

## plot strong filter
plotGO(taaGOqueryStrongFilter, "conserved TAA stopcodons\n-strong filtering -cutoff 1.0")
plotGO(tgaGOqueryStrongFilter, "conserved TGA stopcodons\n-strong filtering -cutoff 1.0")
plotGO(tagGOqueryStrongFilter, "conserved TAG stopcodons\n-strong filtering -cutoff 1.0")
plotGO(nonConservedGOqueryStrongFilter, "non conserved stopcodons\n-no filtering -cutofff 1.0")

dev.off()



