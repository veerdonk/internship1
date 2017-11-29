
# # # TEST # # #
upsetdata <- read.csv("/home/dylan/Documents/david/scripts_git/test.csv", sep=",", header = TRUE)

rodentNames <- c("gene", "known ageing genes", "hamster", "damara mole-rat", "guinea pig", "shrew mouse", "mouse", "deer mouse", "degu", "squirrel", "rat", "kangaroo mouse", "chinchilla")
fishNames <- c("gene", "known ageing gene", "Amazon molly", "Cave fish", "Nile tilapia", "Medaka", "Southern platyfish")
primateNames <- c("gene", "ageing gene", "olive baboon", "squirrel monkey", "green monkey", "gorilla", "tarsier", "macaque")


colnames(upsetdata) <- primateNames
ageingGenes <- upsetdata$gene[upsetdata$ageing_genes]

#tail(names, -1)
#install.packages("UpSetR")
library("UpSetR")
require(ggplot2); require(plyr); require(gridExtra); require(grid);
# , queries = list(list(query = intersects, params = list("ageing.genes", "csa"), active = T))
upset(upsetdata, sets = c("ageing_genes", "csa", "soe", "pan", "ggo", "tsy", "mmu"), keep.order = T, order.by = "freq", nintersects = NA, queries = list(
  list(query=intersects, params= list("ageing_genes", "pan", "ggo"), active=T),
  list(query=intersects, params= list("pan", "ggo", "mmu"),active=T),
  list(query=intersects, params= list("soe", "csa", "tsy"),active=T),
  list(query= elements, params = list("gene", "ADA2", "TEP1", 'FAS', 'FMC1', 'HAP1', 'MRPS16', 'PLAU', 'SWI5'))
))


whale_kg <-  read.csv("/home/dylan/Documents/david/scripts_git/out/upset/whales_kg.csv", sep=",", header = TRUE, row.names = "gene")
upset(whale_kg, nsets = 12, queries = list(
  list(query = intersects, params = list("oro", "phm")),
  list(query = intersects, params = list("tut", "bacu"))
)) # , queries = list(list(query = elements, params = list("gene", "FOXO3")))

whale_all <-  read.csv("/home/dylan/Documents/david/scripts_git/out/upset/whales_kg_allPSG.csv", sep=",", header = TRUE, row.names = "gene")
upset(whale_all, nsets = 12, queries = list(
  list(query = intersects, params = list("oro", "phm")),
  list(query = intersects, params = list("tut", "bacu"))
))

primate_kg <- read.csv("/home/dylan/Documents/david/scripts_git/out/upset/primates_kg.csv", sep = ",", header = TRUE)
upset(primate_kg, nsets = 12, order.by = "freq",queries = list(
  list(query = elements, params = list("ggo", "pan"))
))
primate_all <- read.csv("/home/dylan/Documents/david/scripts_git/out/upset/primates_kg_allPSG.csv", sep = ",", header = TRUE, row.names = "gene")
upset(primate_all, nsets = 12, nintersects = NA, order.by = "freq", decreasing = FALSE, queries = list(
  list(query = elements, params = list("ggo", "pan", "mmu")),
  list(query = elements, params = list("soe", "csa", "tsy"))
))

bats_kg <- read.csv("/home/dylan/Documents/david/scripts_git/out/upset/bats_kg.csv", sep = ",", header = TRUE)
upset(bats_kg, nsets = 12, order.by = "freq",queries = list(
  list(query = elements, params = list("myd"))
))

birds_kg <- read.csv("/home/dylan/Documents/david/scripts_git/out/upset/birds_kg.csv", sep = ",", header = TRUE, row.names = "gene")
upset(birds_kg, nsets = 12)
birds_all <- read.csv("/home/dylan/Documents/david/scripts_git/out/upset/birds_kg_allPSG.csv", sep = ",", header = TRUE)
upset(birds_all, nsets = 12, queries = list(
  list(query = intersects, params = list("gene","krt5"))
))

fish_kg <- read.csv("/home/dylan/Documents/david/scripts_git/out/upset/fish_kg.csv", sep = ",", header = TRUE, row.names = "gene")
upset(fish_kg, nsets = 12,order.by = "freq", decreasing = FALSE)