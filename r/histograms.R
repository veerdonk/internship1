pdf("histograms2.pdf")

# # # files ###
hsa_mmu <- read.csv("/home/dylan/Documents/david/scripts/genenames/out/primates/hsa_mmu_dndsGeneNames.tsv", sep="\t", header = FALSE)
hsa_pan <- read.csv("/home/dylan/Documents/david/scripts_git/out/names/primates/hsa_pan_dndsGeneNames.tsv", sep="\t", header = FALSE)
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
hist(log(hsa_pan[hsa_pan$V4 > 0.001 & hsa_pan$V4 < 99,4]),breaks = 50, freq = FALSE, main="human & baboon dN/dS (PAML)", xlab = "log dnds", ylab = "density")
lines(density(log(na.omit(hsa_pan[hsa_pan$V4 > 0.001 & hsa_pan$V4 < 99, 4]))), lwd = 1, col = "red")
#human-babboon ensembl
ensembl_hsa_pan_dnds <- ensembl_hsa_pan$dN.with.Olive.baboon/ensembl_hsa_pan$dS.with.Olive.baboon
hist(log(ensembl_hsa_pan_dnds),breaks = 50, freq = FALSE, main="human & baboon dN/dS (ensembl)", xlab = "log dnds", ylab = "density")
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

# # # Fish -gaps -stopcodons LOC = reference # # #
par(mfrow=c(2,2))
setwd("~/Documents/david/scripts_git/out")
fish <- dir("./names/fish", full.names = TRUE)
fishTitles <- c("Tetra (cave fish)", "Channel catfish", "Killifish", "Medaka", "Rainbow trout", "Nile tilapia", "Amazon molly", "Southern platyfish")

for (specie in fish){
  lapply(specie, function(x){
    df <- read.csv(x, sep="\t", header = FALSE)
    dndsRatios = df[df$V4 >0.001 & df$V4<99,4]
    hist(log(dndsRatios),breaks = 100, freq = FALSE, main=paste(fishTitles[which(fish == specie)], " - Gar"), xlab = "log dnds", ylab = "density")
    lines(density(log(dndsRatios)), lwd = 1, col = "red")
  })
}

# # # Birds -gaps -stopcodons # # #
setwd("~/Documents/david/scripts_git/out")
birds <- dir("./names/birds/both", full.names = TRUE)
birdTitles <- c("Anna's Hummingbird", "Emperor penguin", "Mallard", "Cuckoo", "Sunbittern", "Collared flycatcher", "Chicken", "Bald eagle", "Turkey", "Adelie penguin", "Dalmation pelican", "Barn owl")


for (bird in birds){
  lapply(bird, function(x){
    df <- read.csv(x, sep="\t", header = FALSE)
    dndsRatios = df[df$V4 >0.001 & df$V4<99,4]
    hist(log(dndsRatios),breaks = 100, freq = FALSE, main=paste(birdTitles[which(birds == bird)]," - Turkey vulture"), xlab = "log dnds", ylab = "density")
    lines(density(log(dndsRatios)), lwd = 1, col = "red")
  })
}

# # # Rodents -gaps -stopcodons # # #
setwd("~/Documents/david/scripts_git/out")
rodents <- dir("./names/rodents/both", full.names = TRUE)
rodentTitles <- c("Beaver", "Brazillian guinea pig", "Long-tailed chinchilla", "Ord's kangaroo rat", "Damara mole-rat", "Golden hamster", "Shrew mouse", "Mouse", "Degu", "NA deer mouse", "Rat", "Squirrel")

for (rodent in rodents){
  lapply(rodent, function(x){
    df <- read.csv(x, sep="\t", header = FALSE)
    dndsRatios = df[df$V4 >0.001 & df$V4<99,4]
    hist(log(dndsRatios),breaks = 100, freq = FALSE, main=paste(rodentTitles[which(rodents == rodent)]," - Naked mole-rat"), xlab = "log dnds", ylab = "density")
    abline(h = 0, v = 0, col = "gray60")
    lines(density(log(dndsRatios)), lwd = 1, col = "red")
  })
}

# # # Whales -gaps -stopcodons # # #
setwd("~/Documents/david/scripts_git/out")
whales <- dir("./names/whales/", full.names = TRUE)
whaleTitles <- c("Minke whale", "Killer whale", "Sperm whale", "Dolphin")

for (whale in whales){
  lapply(whale, function(x){
    df <- read.csv(x, sep="\t", header = FALSE)
    dndsRatios <- df[df$V4 > 0.001 & df$V4 < 99, 4]
    hist(log(dndsRatios), breaks = 100, freq = FALSE, main = paste(whaleTitles[which(whales == whale)], " - Bowhead whale"), xlab = "log dNdS", ylab = "density")
    abline(h = 0, v = 0, col = "gray60")
    lines(density(log(dndsRatios)), lwd = 1, col = "red")
  })
}

# # # Bats -gaps -stopcodons # # #
setwd("~/Documents/david/scripts_git/out")
bats <- dir("./names/bats/", full.names = TRUE)
batsTitles <- c("Big brown bat", "Natal long-fingered bat", "David's myotis", "Black flying fox")

for (bat in bats){
  lapply(bat, function(x){
    df <- read.csv(x, sep="\t", header = FALSE)
    dndsRatios <- df[df$V4 > 0.001 & df$V4 < 99, 4]
    hist(log(dndsRatios), breaks = 100, freq = FALSE, main = paste(batsTitles[which(bats == bat)], " - Brandt's bat"), xlab = "log dNdS", ylab = "density")
    abline(h = 0, v = 0, col = "gray60")
    lines(density(log(dndsRatios)), lwd = 1, col = "red")
  })
}

# # # Primates -gaps -stopcodons # # #
primates <- dir("./names/primates/", full.names = TRUE)
primateTitles <- c("Green monkey", "Gorilla", "Macaque", "Olive baboon", "Squirrel monkey", "Phillipine tarsier")

for (primate in primates){
  lapply(primate, function(x){
    df <- read.csv(x, sep = "\t", header = FALSE)
    dndsRatios <- df[df$V4 > 0.001 & df$V4 < 99, 4]
    hist(log(dndsRatios), breaks = 100, freq = FALSE, main = paste(primateTitles[which(primates == primate)], " - Human"))
    abline(h = 0, v = 0, col = "gray60")
    lines(density(log(dndsRatios)), lwd = 1, col = "red")
  })
}

dev.off()





