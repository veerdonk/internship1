library(dplyr)
library(stringr)
# install.packages("tidytext")
library(tidytext)

data("stop_words")

setwd("~/Documents/david/scripts_git/internship1/r/model/")
genelist <- read.delim("./results/pantherGeneList.txt", header = FALSE, stringsAsFactors = FALSE)
sigGenes <- read.csv("/home/dylan/Documents/david/scripts_git/internship1/r/model/results/B_All Data.csv")
colnames(sigGenes) <- c("gene", "pval", "adj_pval")

ageing_genes <- read.csv("~/Documents/david/resources/ageing_genes/ageing_genes.csv")

tidyGenes <- genelist %>%
  unnest_tokens(word, V3)

tidyDescriptions <- genelist%>%
  unnest_tokens(word, V4)

tidyGenes <- tidyGenes %>%
  anti_join(stop_words)

tidyDescriptions <- tidyDescriptions %>%
  anti_join(stop_words)

bio_terms <- data.frame(c("ortholog", "protein", "1", "2", "3", "4", "domain", "cell", "5", "6", "7", "8", "type"))
colnames(bio_terms) <- "word"
tidyGenes <- tidyGenes %>%
  anti_join(bio_terms)
tidyDescriptions <- tidyDescriptions %>%
  anti_join(bio_terms)

most_common <- tidyGenes %>%
  count(word, sort = TRUE) 

wanted_terms <- c("[cC]ancer","[oO]nco", "[tT]umor", "[tT]elom", "[Aa]ggregat","[Rr]epair")#, "[Rr]egulat")#"[eE]xtracellular matrix", 

genelist[genelist$V2 %in% tidyGenes[tidyGenes$word %in% wanted_terms,2],2] 
genelist[genelist$V2 %in% tidyDescriptions[tidyDescriptions$word %in% wanted_terms,2],2]

genesOfInterest <- c()
for(term in wanted_terms){
  genesOfInterest <- c(genesOfInterest, genelist[grep(term, genelist$V3),2])
  genesOfInterest <- c(genesOfInterest, genelist[grep(term, genelist$V4),2])
}

genesOfInterestInfo <- genelist[genelist$V2 %in% genesOfInterest,]
genesInAgeing <- genelist[genelist$V2 %in% ageing_genes$symbol,]

sigGenes[sigGenes$gene %in% genesOfInterestInfo$V2,]












