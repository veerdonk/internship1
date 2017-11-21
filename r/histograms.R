setwd("~/Documents/david/scripts/genenames/out")
pdf("histograms.pdf")
rodentFiles <- dir("./rodents", full.names = TRUE)
fishFiles <- dir("./fish_gar", full.names = TRUE)
whaleFiles <- dir("./whales", full.names = TRUE)
birdFiles <- dir("./birds", full.names = TRUE)
batFiles <- dir("./bats", full.names = TRUE)
primateFiles <- dir("./primates", full.names = TRUE)
# species <- c(rodentFiles, fishFiles, whaleFiles, birdFiles, batFiles, primateFiles)
species <- c(primateFiles)
# # # files ###
hsa_mmu <- read.csv("/home/dylan/Documents/david/scripts/genenames/out/primates/hsa_mmu_dndsGeneNames.tsv", sep="\t", header = FALSE)
hsa_pan <- read.csv("/home/dylan/Documents/david/scripts/genenames/out/primates/hsa_pan_dndsGeneNames.tsv", sep="\t", header = FALSE)
ensembl_hsa_pan <- read.csv("/home/dylan/Documents/david/resources/primates/hsa_pan_dn_ds.txt", sep="\t")
human_pan_nostop <- read.csv("../human_pan_dndsGeneNames.tsv", sep="\t", header = FALSE)
ensemblHuman <- read.csv("/home/dylan/Documents/david/resources/primates/human_macaque_dn_ds.txt", sep="\t")

hgl_mus <- read.csv("/home/dylan/Documents/david/resources/rodents/dnds/hgl_mus.txt", sep="\t")
hgl_mus_paml <- read.csv("/home/dylan/Documents/david/scripts/genenames/out/rodents/hgl_mus_dndsGeneNames.tsv", sep="\t", header = FALSE)

hgl_mus_omega5 <- read.csv("/home/dylan/Documents/david/scripts/genenames/test/hgl_mus_dnds_high_omegaGeneNames.tsv", sep = "\t", header = FALSE)
hgl_mus_omega1 <- read.csv("/home/dylan/Documents/david/scripts/genenames/test/hgl_mus_dnds_omega1GeneNames.tsv", sep = "\t", header = FALSE)

### test with dnds ratios from ensembl ###
# human - macaque ensembl
ensembldnds <- ensemblHuman$dN.with.Macaque/ensemblHuman$dS.with.Macaque
hist(log(ensembldnds),breaks = 100, freq = FALSE, main="human & macaque dN/dS (ensembl)", xlab = "log dnds", ylab = "density")
lines(density(log(na.omit(ensembldnds))), lwd = 1, col = "red")
#human macaque calculated
hist(log(hsa_mmu[hsa_mmu$V4 > 0.001 & hsa_mmu$V4 < 99,4]),breaks = 100, freq = FALSE, main="human & macaque dN/dS (PAML)", xlab = "log dnds", ylab = "density")
lines(density(log(na.omit(hsa_mmu[hsa_mmu$V4 > 0.001 & hsa_mmu$V4 < 99, 4]))), lwd = 1, col = "red")

#human babboon calculated
hist(log(hsa_mmu[hsa_mmu$V4 > 0.001 & hsa_mmu$V4 < 99,4]),breaks = 100, freq = FALSE, main="human & baboon dN/dS (PAML)", xlab = "log dnds", ylab = "density")
lines(density(log(na.omit(hsa_mmu[hsa_mmu$V4 > 0.001 & hsa_mmu$V4 < 99, 4]))), lwd = 1, col = "red")
#human-babboon ensembl
ensembl_hsa_pan_dnds <- ensembl_hsa_pan$dN.with.Olive.baboon/ensembl_hsa_pan$dS.with.Olive.baboon
hist(log(ensembl_hsa_pan_dnds),breaks = 100, freq = FALSE, main="human & baboon dN/dS (ensembl)", xlab = "log dnds", ylab = "density")
lines(density(log(na.omit(ensembl_hsa_pan_dnds))), lwd = 1, col = "red")
#human baoon calculated without stop/gap
hist(log(human_pan_nostop[human_pan_nostop$V4 > 0.001 & human_pan_nostop$V4 < 99,4]),breaks = 100, freq = FALSE, main="human & baboon dN/dS (PAML-no stop/gap)", xlab = "log dnds", ylab = "density")
lines(density(log(na.omit(human_pan_nostop[human_pan_nostop$V4 > 0.001 & human_pan_nostop$V4 < 99, 4]))), lwd = 1, col = "red")

#mole rat - mouse calculated
hist(log(hgl_mus_paml[hgl_mus_paml$V4>0.001&hgl_mus_paml$V4<99,4]), breaks = 100, freq = FALSE, main="mole rat & mouse dN/dS (PAML)", xlab = "log dnds", ylab = "density")
lines(density(log(hgl_mus_paml[hgl_mus_paml$V4>0.001&hgl_mus_paml$V4<99,4])), lwd=1, col="red")
range(hgl_mus_paml[hgl_mus_paml$V4>0.001&hgl_mus_paml$V4<99,4])

#mole rat - mouse ensembl
ensembl_hgl_mus <- hgl_mus$dN.with.Mouse/hgl_mus$dS.with.Mouse
hist(log(ensembl_hgl_mus),breaks = 100, freq = FALSE, main="mole rat & mouse dN/dS (Ensembl)", xlab = "log dnds", ylab = "density")
lines(density(log(na.omit(ensembl_hgl_mus))), lwd = 1, col = "red")
mean(log(na.omit(ensembl_hgl_mus)))

omittedNas <- na.omit(ensembl_hgl_mus)
na.omit(ensembl_hgl_mus[ensembl_hgl_mus>1])

#molerat - mouse without stopcodons/gaps
hgl_mus_nostop <- read.csv("/home/dylan/Documents/david/scripts/genenames/moleRat_mouse_dndsGeneNames.tsv", sep="\t", header = FALSE)
hist(log(hgl_mus_nostop[hgl_mus_nostop$V4>0.001 & hgl_mus_nostop$V4<99,4]), breaks = 100, freq = FALSE, main = "naked molerat & mouse (PAML - no stop/gap)", xlab = "log dnds", ylab = "density")
lines(density(log(hgl_mus_nostop[hgl_mus_nostop$V4>0.001 & hgl_mus_nostop$V4<99,4])), lwd=1, col="red")
range(hgl_mus_nostop[hgl_mus_nostop$V4>0.001 & hgl_mus_nostop$V4<99,4])
mean(log(hgl_mus_nostop[hgl_mus_nostop$V4>0.001 & hgl_mus_nostop$V4<99,4]))
hgl_mus_nostop[hgl_mus_nostop$V4>8 & hgl_mus_nostop$V4<99,4]

### fish with loc as reference ###
#loc - oni
loc_oni <- read.csv("/home/dylan/Documents/david/scripts/genenames/out/fish_gar/loc_oniGeneNames.tsv", sep="\t", header = FALSE)
loc_oni_dnds <- loc_oni[loc_oni$V4>0.001 & loc_oni$V4<99, 4]
hist(log(loc_oni_dnds), breaks = 100, freq = FALSE, main = "Spotted gar - Nile tilapia", xlab = "log dnds", ylab = "density")
lines(density(log(loc_oni_dnds)), lwd=1, col="red")
#loc - pfo
loc_pfo <- read.csv("/home/dylan/Documents/david/scripts/genenames/out/fish_gar/loc_pfoGeneNames.tsv", sep = "\t", header = FALSE)
loc_pfo_dnds <- loc_pfo[loc_pfo$V4>0.001 & loc_pfo$V4<99, 4]
hist(log(loc_pfo_dnds), breaks = 100, freq = FALSE, main = "Spotted gar - Amazon molly", xlab = "log dnds", ylab = "density")
lines(density(log(loc_pfo_dnds)), lwd=1, col="red")

# # # LOC = reference # # #
setwd("~/Documents/david/scripts_git/out")
fish <- dir("./names/fish", full.names = TRUE)
fishTitles <- c("Gar - Tetra (cave fish)", "Gar - Channel catfish", "Gar - Killifish", "Gar - Medaka", "Gar - Rainbow trout", "Gar - Nile tilapia", "Gar - Amazon molly", "Gar - Southern platyfish")



for (specie in fish){
  lapply(specie, function(x){
    df <- read.csv(x, sep="\t", header = FALSE)
    dndsRatios = df[df$V4 >0.001 & df$V4<99,4]
    hist(log(dndsRatios),breaks = 100, freq = FALSE, main=fishTitles[which(fish == specie)], xlab = "log dnds", ylab = "density")
    lines(density(log(dndsRatios)), lwd = 1, col = "red")
  })
}

# # # Birds -gaps -stopcodons # # #
setwd("~/Documents/david/scripts_git/out")
birds <- dir("./names/birds", full.names = TRUE)
birdTitles <- c("Anna's Hummingbird", "Emperor penguin", "Mallard", "Cuckoo", "Sunbittern", "Collared flycatcher", "Chicken", "Bald eagle", "Turkey", "Adelie penguin", "Dalmation pelican", "Barn owl")


for (bird in birds){
  lapply(bird, function(x){
    df <- read.csv(x, sep="\t", header = FALSE)
    dndsRatios = df[df$V4 >0.001 & df$V4<99,4]
    hist(log(dndsRatios),breaks = 100, freq = FALSE, main=paste(birdTitles[which(birds == bird)]," - Turkey vulture"), xlab = "log dnds", ylab = "density")
    lines(density(log(dndsRatios)), lwd = 1, col = "red")
  })
}

# # # TEST # # #
upsetdata <- read.csv("/home/dylan/Documents/david/scripts/genenames/test.csv", sep=",", header = TRUE)

rodentNames <- c("gene", "known ageing genes", "hamster", "damara mole-rat", "guinea pig", "shrew mouse", "mouse", "deer mouse", "degu", "squirrel", "rat", "kangaroo mouse", "chinchilla")
fishNames <- c("gene", "known ageing gene", "Amazon molly", "Cave fish", "Nile tilapia", "Medaka", "Southern platyfish")
colnames(upsetdata) <- fishNames


#tail(names, -1)
#install.packages("UpSetR")
library("UpSetR")
require(ggplot2); require(plyr); require(gridExtra); require(grid);
upset(upsetdata, sets = fishNames[2:7], nintersects = NA, order.by = "freq")


movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=T, sep=";" )
upset(movies, sets = c("Action", "Adventure", "Children", "War", "Noir"),
      queries = list(list(query = intersects, params = list("War"), active = T),
                     list(query = intersects, params = list("Adventure", "Action"))))


dev.off()





