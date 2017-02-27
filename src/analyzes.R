#Convert the data to GenABEL format
require(GenABEL)
load("../data/genotypes.RData")
load("../data/phenotypes.RData")

#ped file
gdata.ped <- gdata
gdata.ped[gdata.ped == 1] <- 2
gdata.ped[gdata.ped == -1] <- 1
gdata.ped <- gdata.ped[, rep(x = colnames(gdata.ped), times = rep(x = 2, times = ncol(gdata.ped)))] #duplicate every column in the genotype matrix
ids <- rownames(gdata.ped)
gdata.ped <- data.frame(ids, ids, "0 0 0 -9", gdata.ped)
write.table(x = gdata.ped, file = "../data/gdata.ped", sep = "\t", quote = F, row.names = F, col.names = F)

#map file
markers <- colnames(gdata)
tmp <- strsplit(markers, split = "_")
map_SacCer3 <- as.numeric(sapply(X = tmp, FUN = function(x){x[3]}))
chr <- sapply(X = tmp, FUN = function(x){x[2]})
tmp <- unique(chr)
for(i in 1:length(tmp)){
  chr[chr == tmp[i]] <- i
}
chr <- as.numeric(chr)
gdata.map <- data.frame(chrom = chr, name = markers, position = map_SacCer3)
write.table(x = gdata.map, file = "../data/gdata.map", sep = "\t", quote = F, row.names = F, col.names = T)

#phenotype file
pheno <- data.frame(id = ids, sex = 1, matrix(nrow = length(pheno_raw[[1]]), ncol = 60))
names <- c() #the names of the traits
names[seq(from = 1, to = 58, by = 3)] <- paste(paste("trait", 1:20, sep = "_"), 1, sep = ".")
names[seq(from = 2, to = 59, by = 3)] <- paste(paste("trait", 1:20, sep = "_"), 2, sep = ".")
names[seq(from = 3, to = 60, by = 3)] <- paste(paste("trait", 1:20, sep = "_"), "mean", sep = ".")
names(pheno)[3:62] <- names
j <- 3
for(i in 1:length(pheno_raw)){
  trait <- pheno_raw[[i]]
  trait1 <- sapply(X = trait, FUN = function(x){x[1]})
  trait2 <- sapply(X = trait, FUN = function(x){x[2]})
  pheno[,j] <- trait1
  j <- j+1
  pheno[,j] <- trait2
  j <- j+1
  pheno[,j] <- sapply(X = trait, FUN = function(x){ mean(c(x[1], x[2]), na.rm = T) })
  j <- j+1
}
write.table(x = pheno, file = "../data/phenotypes_meanAndRep", quote = F, sep = "\t", row.names = F)

#Create raw file and load it as a GenABEL object
convert.snp.ped(pedfile = "../data/gdata.ped", mapfile = "../data/gdata.map", outfile = "../data/gdata.raw")
data <- load.gwaa.data(phenofile = "../data/phenotypes_meanAndRep", genofile = "../data/gdata.raw")
save(list = "data", file = "../data/data_GenABEL.RData")



#Build networks
require(igraph)
bloom.int <- read.table(file = "../data/bloom2015_int.txt", sep = "\t", header = T) #The pairwise interactions mapped in Bloom2015
n <- 50 #consider SNPs within +-n to be the same
traits <- unique(as.character(bloom.int$trait))
networks.index <- as.list(1:length(traits))
networks.ids <- as.list(1:length(traits))
names(networks.index) <- traits
names(networks.ids) <- traits
for(trait in traits){
  tmp <- bloom.int[bloom.int$trait == trait, ]
  snps <- c(tmp$Q1_index, tmp$Q2_index)
  pos <- c(tmp$Q1_pos, tmp$Q2_pos)
  diffs <- abs(outer(snps, snps, FUN = "-")) #all pairwise differences in index
  diffs[lower.tri(diffs)] <- 1000
  diag(diffs) <- 1000
  clump <- which(diffs < n) #differences smaller than n
  diffs.col <- ceiling(clump/ncol(diffs)) #differences smaller than n, columns in the distance matrix
  diffs.row <- clump - (diffs.col - 1)*nrow(diffs) #differences smaller than n, rows in the distance matrix
  diffs.graph <- graph_from_data_frame(data.frame(diffs.col, diffs.row))
  diffs.clust <- clusters(diffs.graph)
  for(i in 1:diffs.clust$no){
    nodesInCluster <- as.numeric(names(diffs.clust$membership[diffs.clust$membership == i]))
    snps[nodesInCluster] <- snps[nodesInCluster[1]]
    pos[nodesInCluster] <- pos[nodesInCluster[1]]
  }
  tmp.clumped <- tmp
  tmp.clumped$Q1_index <- snps[1:nrow(tmp)]
  tmp.clumped$Q2_index <- snps[(nrow(tmp) + 1):(2*nrow(tmp))]
  tmp.clumped$Q1_pos <- pos[1:nrow(tmp)]
  tmp.clumped$Q2_pos <- pos[(nrow(tmp) + 1):(2*nrow(tmp))]
  
  relations <- data.frame(from = paste(tmp.clumped$Q1_chr, tmp.clumped$Q1_pos), to = paste(tmp.clumped$Q2_chr, tmp.clumped$Q2_pos), LOD = tmp$LOD)
  relations.graph <- graph_from_data_frame(relations, directed = F)
  networks.ids[[trait]] <- relations.graph
  networks.index[[trait]] <- graph_from_data_frame(data.frame(tmp.clumped$Q1_index, tmp.clumped$Q2_index, LOD = tmp$LOD), directed = F)
}
save(list = c("networks.ids", "networks.index"), file = "../data/networks.RData")



#Count the number of interactors of every epistatic QTL
require(igraph)
data.traitNames <- names(pheno_raw)
nrVertices <- sum(sapply(networks.ids, FUN = function(x){length(V(x))}))
networkNeighbours <- data.frame(snpName = rep(NA, nrVertices), snpIndex = NA, nrNeighbours = NA, trait = NA)
#networkNeighbours_effectSize <- data.frame(snpName = rep(NA, nrVertices), snpIndex = NA, nrNeighbours = NA, vGWA_dGLM_eff = NA, vGWA_dGLM_p = NA, trait = NA)
for(trait in names(networks.index)){
  #Count the number of neighbours of every loci in the network
  network <- networks.index[[trait]]
  network.vertices <- as.numeric(V(network)$name)
  network.vertices_neighbours <- sapply(1:length(network.vertices), FUN = function(x){ length(network.vertices[neighbors(graph = network, v = x, mode = "total")]) })
  
  from <- min(which(is.na(networkNeighbours$snpName)))
  to <- from + length(network.vertices) - 1

  networkNeighbours$snpName[from:to] <- data@gtdata@snpnames[network.vertices]
  networkNeighbours$snpIndex[from:to] <- network.vertices
  networkNeighbours$nrNeighbours[from:to] <- network.vertices_neighbours
  networkNeighbours$trait[from:to] <- trait
}
save(list = "networkNeighbours", file = "../data/networkInfo.RData")



#Check the capacitating effect of every hub-QTL through the following procedure:
#1. Divide the population in two based on the genotype at the hub
#2. Calculate the narrow sense heritability h2 in each group
load("../data/data_GenABEL.RData")
load("../data/networkInfo.RData")
load("../data/phenotypes.RData")
networkNeighbours.hubs <- networkNeighbours[networkNeighbours$nrNeighbours > 4,]

library(doSNOW)
library(foreach)
cl <- makeCluster(11)
registerDoSNOW(cl)
ptm <- proc.time()
h2Split_polygenic_hubs <- foreach(i = 1:nrow(networkNeighbours.hubs)) %dopar% {
  library(GenABEL)
  snp <- networkNeighbours.hubs$snpName[i]
  trait <- paste("trait_", grep(pattern = networkNeighbours.hubs$trait[i], x = names(pheno_raw)), ".mean", sep = "")
  geno <- as.double.gwaa.data(data[, snp])
  data.group1 <- data[geno[,1] == 0,]
  data.group2 <- data[geno[,1] == 2,]
  
  #mixed model in group 1
  relmat <- ibs(data.group1)
  data.group1_polygenic <- try(polygenic(formula = data.group1@phdata[,trait], kinship.matrix = relmat, data = data.group1))
  if(class(data.group1_polygenic) == "try-error")
    group1_h2 <- NA
  else
    group1_h2 <- data.group1_polygenic$esth2
  
  #mixed model in group 2
  relmat <- ibs(data.group2)
  data.group2_polygenic <- try(polygenic(formula = data.group2@phdata[,trait], kinship.matrix = relmat, data = data.group2))
  if(class(data.group2_polygenic) == "try-error")
    group2_h2 <- NA
  else
    group2_h2 <- data.group2_polygenic$esth2
  
  return(list(group1_h2, group2_h2))
}
proc.time() - ptm
stopCluster(cl)
save(list = "h2Split_polygenic_hubs", file = "../data/hubSplits_h2Polygenic.RData")



#Same procedure for a random subset of radial QTL
load("../data/data_GenABEL.RData")
load("../data/networkInfo.RData")
load("../data/phenotypes.RData")
tmp <- which(networkNeighbours$nrNeighbours == 1)
radialSet <- sample(x = tmp, size = 40)
library(doSNOW)
library(foreach)
cl <- makeCluster(30)
registerDoSNOW(cl)
ptm <- proc.time()
h2Split_polygenic_radial40 <- foreach(i = 1:length(radialSet)) %dopar% {
  library(GenABEL)
  snp <- networkNeighbours$snpName[ radialSet[i] ]
  trait <- paste("trait_", grep(pattern = networkNeighbours$trait[ radialSet[i] ], x = names(pheno_raw)), ".mean", sep = "")
  geno <- as.double.gwaa.data(data[, snp])
  data.group1 <- data[geno[,1] == 0,]
  data.group2 <- data[geno[,1] == 2,]
  
  #mixed model in group 1
  relmat <- ibs(data.group1)
  data.group1_polygenic <- try(polygenic(formula = data.group1@phdata[,trait], kinship.matrix = relmat, data = data.group1))
  if(class(data.group1_polygenic) == "try-error")
    group1_h2 <- NA
  else
    group1_h2 <- data.group1_polygenic$esth2
  
  #mixed model in group 2
  relmat <- ibs(data.group2)
  data.group2_polygenic <- try(polygenic(formula = data.group2@phdata[,trait], kinship.matrix = relmat, data = data.group2))
  if(class(data.group2_polygenic) == "try-error")
    group2_h2 <- NA
  else
    group2_h2 <- data.group2_polygenic$esth2
  
  return(list(group1_h2, group2_h2, snp, trait))
}
proc.time() - ptm
stopCluster(cl)
save(list = "h2Split_polygenic_radial40", file = "../data/40radialSplits_h2Polygenic.RData")



#Permutation test for the capacitation analysis
#This example gives an empirical NULL distribution for the first trait. The same procedure needs to be repeated to create a NULL distribution for each trait (change the trait variable accordingly)
n <- 1000 #number of permutations per trait
trait <- "trait_1.mean"
library(doSNOW)
library(foreach)
cl <- makeCluster(33)
registerDoSNOW(cl)

h2SplitPermut_polygenic_trait1 <- foreach(i = 1:n) %dopar% {
  require(GenABEL)
  load("data_GenABEL.RData")
  
  split <- sample(x = 1:2, size = data@gtdata@nids, replace = T)
  
  group1 <- data@gtdata@idnames[split == 1]
  group2 <- data@gtdata@idnames[split == 2]
  trait.notMissing <- data@gtdata@idnames[!is.na(data@phdata[, trait])]
  
  #mixed model in group 1
  group1.notMissing <- group1[group1 %in% trait.notMissing]
  data.group1 <- data[group1.notMissing, ]
  relmat <- ibs(data.group1)
  data.group1_polygenic <- try(polygenic(formula = data.group1@phdata[,trait], kinship.matrix = relmat, data = data.group1))
  if(class(data.group1_polygenic) == "try-error")
    group1_h2 <- NA
  else
    group1_h2 <- data.group1_polygenic$esth2
  
  #mixed model in group 2
  group2.notMissing <- group2[group2 %in% trait.notMissing]
  data.group2 <- data[group2.notMissing, ]
  relmat <- ibs(data.group2)
  data.group2_polygenic <- try(polygenic(formula = data.group2@phdata[,trait], kinship.matrix = relmat, data = data.group2))
  if(class(data.group2_polygenic) == "try-error")
    group2_h2 <- NA
  else
    group2_h2 <- data.group2_polygenic$esth2
  
  return(list(group1_h2, group2_h2))
}
stopCluster(cl)
save(list = "h2SplitPermut_polygenic_trait1", file = "../data/randomSplitTrait1_polygenic_n1000.RData")



#Test significance of the capacitation
for(i in 1:20){
  load(paste("../data/randomSplitTrait", i, "_polygenic_n1000.RData", sep = ""))
}
load("../data/networkInfo.RData")
load("../data/hubSplits_h2Polygenic.RData")
load("../data/phenotypes.RData")
data.traitNames <- names(pheno_raw)
networkNeighbours.hubs <- networkNeighbours[networkNeighbours$nrNeighbours > 4,]

#Skip the hubs that give the completely overlapping CopperSulfate networks
networkNeighbours.hubs <- networkNeighbours.hubs[c(1:2, 7:19), ] 
h2Split_polygenic_hubs <- h2Split_polygenic_hubs[c(1:2, 7:19)]

hubs.capaciPvals_n1000 <- rep(NA, nrow(networkNeighbours.hubs))
hubs.h2Diff <- rep(NA, nrow(networkNeighbours.hubs))
for(i in 1:nrow(networkNeighbours.hubs)){
  trait <- networkNeighbours.hubs$trait[i]
  data.traitNr <- grep(pattern = trait, x = data.traitNames)
  
  #The observed difference in h2
  tmp <- h2Split_polygenic_hubs[[i]]
  hub.h2Diff <- abs(tmp[[1]] - tmp[[2]])
  hubs.h2Diff[i] <- hub.h2Diff
  
  #The permutation differences
  tmp <- get(paste("h2SplitPermut_polygenic_trait", data.traitNr, sep = ""))
  null.h2Diff <- unlist(lapply(X = tmp, FUN = function(x){abs(x[[1]] - x[[2]])}))
  
  #p-value
  missing <- is.na(null.h2Diff)
  hubs.capaciPvals_n1000[i] <- sum(hub.h2Diff < null.h2Diff[!missing])/sum(!missing)
}
sum(hubs.capaciPvals_n1000 < .05/15)
hubs.summary <- data.frame(networkNeighbours.hubs, h2Diff = hubs.h2Diff, h2Diff_p = hubs.capaciPvals_n1000)
save(list = "hubs.summary", file = "../data/hubs_capacEff.RData")



#Perform IAA GWAS in the segregants with BY and RM alleles at the hub separately
#Compare the network implicated by this to the one obtained in the epistatic scan
load("../data/data_GenABEL.RData")
load("../data/networks.RData")
iaaRadialNetwork <- as.numeric(names(neighbors(graph = networks.index$IndolaceticAcid, v = "11146")))
iaaRadialNetwork.pos <- as.numeric(gsub(pattern = "(.*?)_.*", replacement = "\\1", x = data@gtdata@snpnames[iaaRadialNetwork]))
geno <- as.double.gwaa.data(data[, 11146])
iaaScan_capac <- qtscore(trait_8.mean, data[geno[,1] == 0, ])
iaaScan_canal <- qtscore(trait_8.mean, data[geno[,1] == 2, ])

geno.coVar <- as.double.gwaa.data(data[geno[,1] == 0, c("2358650_chrIV_998628_A_T", "9714664_chrXIV_467028_A_G", "1242017_chrIII_198615_T_G", "8733525_chrXIII_410320_T_C")])
iaaScan_capacCovar <- qtscore(data@phdata[geno[,1] == 0, "trait_8.mean"] ~ geno.coVar[,1] + geno.coVar[,2] + geno.coVar[,3] + geno.coVar[,4], data[geno[,1] == 0, ])

geno.coVar <- as.double.gwaa.data(data[geno[,1] == 2, c("2358650_chrIV_998628_A_T", "9714664_chrXIV_467028_A_G", "1242017_chrIII_198615_T_G", "8733525_chrXIII_410320_T_C")])
iaaScan_canalCovar <- qtscore(data@phdata[geno[,1] == 2, "trait_8.mean"] ~ geno.coVar[,1] + geno.coVar[,2] + geno.coVar[,3] + geno.coVar[,4], data[geno[,1] == 2, ])

pos <- as.numeric(gsub(pattern = "(.*?)_.*", replacement = "\\1", x = data@gtdata@snpnames))
ymax <- max(c(-log10(iaaScan_capac[,"P1df"]), -log10(iaaScan_capacCovar[,"P1df"])))
par(mfrow = c(3,1), mar = c(5,4,2,2)+.1)
col <- rep("black", data@gtdata@nsnps)
col[iaaScan_canal[,"P1df"] < .05/data@gtdata@nsnps] <- "red"
plot(pos, -log10(iaaScan_canal[,"P1df"]), pch = 19, cex = .5, bty = "n", ylim = c(0, ymax), main = "Canalized", col = col)
abline(v = iaaRadialNetwork.pos, lty = 1, lwd = 2)
abline(h = -log10(.05/data@gtdata@nsnps), lty = 2, col = "red")

col <- rep("black", data@gtdata@nsnps)
col[iaaScan_capac[,"P1df"] < .05/data@gtdata@nsnps] <- "red"
plot(pos, -log10(iaaScan_capac[,"P1df"]), pch = 19, cex = .5, bty = "n", ylim = c(0, ymax), main = "Capacitated", col = col)
abline(v = iaaRadialNetwork.pos, lty = 1, lwd = 2)
abline(h = -log10(.05/data@gtdata@nsnps), lty = 2, col = "red")

col <- rep("black", data@gtdata@nsnps)
col[iaaScan_capacCovar[,"P1df"] < .05/data@gtdata@nsnps] <- "red"
plot(pos, -log10(iaaScan_capacCovar[,"P1df"]), pch = 19, cex = .5, bty = "n", ylim = c(0, ymax), main = "Capac. 4 covar", col = col)
abline(v = iaaRadialNetwork.pos, lty = 1, lwd = 2)
abline(h = -log10(.05/data@gtdata@nsnps), lty = 2, col = "red")



#Quantify the synergistic effect of radial QTLs, when combined with the capacitating allele at the corresponding hub QTL
#This is done by comparing the fit of additive and exponential models, fitted to the 5 radial QTLs in the capacitated group of segregants
networkNeighbours.hubs <- networkNeighbours[networkNeighbours$nrNeighbours > 4,]
#Skip the hubs that give the completely overlapping CopperSulfate networks
networkNeighbours.hubs <- networkNeighbours.hubs[c(1:2, 7:19), ]

#Copper Sulfate
hub <- networkNeighbours.hubs$snpIndex[networkNeighbours.hubs$trait == "CopperSulfate"]
network <- networks.index$CopperSulfate

#Pick the top 5 interactions around the hub (Really messy. There should be a simpler way if you understand igraph)
hubNeighbors <- E(network)[from(as.character(hub))] #the neighboring edges
tmp <- order(hubNeighbors$LOD, decreasing=TRUE)[1:5] #sort according to LOD
hubNeighbors.top <- E(network)[ as.vector(hubNeighbors)[tmp] ] #sorted edge list
tmp2 <- get.edges(graph = network, hubNeighbors.top) #indices of the nodes in hubNeighbors.top
hubNeighbors.names <- as.numeric(V(network)[unique(c(tmp2))]$name) #take out these nodes from network

pheno <- data@phdata$trait_2.mean
names(pheno) <- data@gtdata@idnames
geno <- as.double.gwaa.data(data[,hubNeighbors.names])
geno[geno == 2] <- 1
#Code the genotypes at the radial QTLs so that 1 is the allele with a lower phenotype
for(i in 2:ncol(geno)){ 
  pheno.split <- split(pheno, geno[,i])
  if(mean(pheno.split$`1`, na.rm = T) > mean(pheno.split$`0`, na.rm = T)){
    geno[,i] <- geno[,i]+1
    geno[geno[,i] == 2, i] <- 0
  }
}

group.RM <- rownames(geno)[geno[,1] == 1]
group.BY <- rownames(geno)[geno[,1] == 0]
#Make allele-count variables
tmp <- geno[group.RM, 2:6]
group.RM.count <- apply(X = tmp, MARGIN = 1, FUN = sum)
tmp <- geno[group.BY, 2:6]
group.BY.count <- apply(X = tmp, MARGIN = 1, FUN = sum)

#Fit additive and exponential models to the allele count
group.RM.count_lmAdd <- lm(pheno[group.RM] ~ group.RM.count) #additive
tmp <- data.frame(pheno = pheno[group.RM], count = group.RM.count)
tmp <- tmp[!is.na(tmp$pheno), ]
group.RM.count_nlmExp <- summary(nls(pheno ~ a + b*exp(c*count), start = list(a = 0, b = 1, c = 1), data = tmp))

#plot the fit of the two models
ylim <- c(min(pheno, na.rm = T), max(pheno, na.rm = T))
box.groupRM <- boxplot(pheno[group.RM] ~ group.RM.count, boxwex = .15, at = seq(.9, 5.9), col = "green", ylim = ylim, xaxt = "n", frame = F, yaxt = "n")
box.groupBY <- boxplot(pheno[group.BY] ~ group.BY.count, boxwex = .15, at = seq(1.1, 6.1), add = T, col = "grey", ylim = ylim, xaxt = "n", frame = F, yaxt = "n")
axis(side = 1, at = 1:6, labels = 0:5, cex.axis = 2, padj = .2)
axis(side = 2, cex.axis = 2)
mtext("Number of growth decreasing alleles", side = 1, cex = 2, line = 3)
mtext("Colony Growth CopperSulfate", side = 2, cex = 2, line = 2.5)
legend("topright", c("BY", "RM"), col = c("grey", "green"), pch = 19, title = "Hub allele")

#Add the two regression lines to the plot
a <- group.RM.count_nlmExp$coefficients[1,1]
b <- group.RM.count_nlmExp$coefficients[2,1]
c <- group.RM.count_nlmExp$coefficients[3,1]
f_exp <- function(x){a + b*exp(c*x) }
lines(x = seq(.9, 5.9), y = f_exp(0:5), lwd = 2, col = "blue", lty = 2)
lines(x = c(.9, 5.9), y = c(group.RM.count_lmAdd$coefficients[1], group.RM.count_lmAdd$coefficients[1] + group.RM.count_lmAdd$coefficients[2]*5), lwd = 2)
R2_exp <- 1 - var(pheno[group.RM] - f_exp(group.RM.count), na.rm = T)/var(pheno[group.RM], na.rm = T)
R2_add <- 1 - var(group.RM.count_lmAdd$residuals, na.rm = T)/var(pheno[group.RM], na.rm = T)
legend("bottomleft", c(paste("Additive model       (R2 = ", round(R2_add, digits = 3), ")", sep = ""), paste("Exponential model (R2 = ", round(R2_exp, digits = 3), ")", sep = "")), col = c("black", "blue"), title = "Regression Lines", lty = "solid", lwd = 3)

#IAA
pheno <- data@phdata$trait_8.mean
names(pheno) <- data@gtdata@idnames
geno <- as.double.gwaa.data(data[, c("4944074_chrVIII_114144_A_G", "9716860_chrXIV_469224_A_G", "1242017_chrIII_198615_T_G", "2358650_chrIV_998628_A_T", "8733525_chrXIII_410320_T_C", "7890567_chrXII_645539_A_G")])
geno[geno == 2] <- 1
#Code the genotypes at the radial QTLs so that 1 is the allele with a lower phenotype
for(i in 2:ncol(geno)){ 
  pheno.split <- split(pheno, geno[,i])
  if(mean(pheno.split$`1`, na.rm = T) > mean(pheno.split$`0`, na.rm = T)){
    geno[,i] <- geno[,i]+1
    geno[geno[,i] == 2, i] <- 0
  }
}

group.RM <- rownames(geno)[geno[,1] == 0]
group.BY <- rownames(geno)[geno[,1] == 1]
#Make allele-count variables
tmp <- geno[group.RM, 2:6]
group.RM.count <- apply(X = tmp, MARGIN = 1, FUN = sum)
tmp <- geno[group.BY, 2:6]
group.BY.count <- apply(X = tmp, MARGIN = 1, FUN = sum)

#Fit additive and exponential models to the allele count
group.RM.count_lmAdd <- lm(pheno[group.RM] ~ group.RM.count) #additive
tmp <- data.frame(pheno = pheno[group.RM], count = group.RM.count)
tmp <- tmp[!is.na(tmp$pheno), ]
group.RM.count_nlmExp <- summary(nls(pheno ~ a + b*exp(c*count), start = list(a = 0, b = 1, c = 1), data = tmp))

#plot the fit of the two models
ylim <- c(min(pheno, na.rm = T), max(pheno, na.rm = T))
box.groupRM <- boxplot(pheno[group.RM] ~ group.RM.count, boxwex = .15, at = seq(.9, 5.9), col = "green", ylim = ylim, xaxt = "n", frame = F, yaxt = "n")
box.groupBY <- boxplot(pheno[group.BY] ~ group.BY.count, boxwex = .15, at = seq(1.1, 6.1), add = T, col = "grey", ylim = ylim, xaxt = "n", frame = F, yaxt = "n")
axis(side = 1, at = 1:6, labels = 0:5, cex.axis = 2, padj = .2)
axis(side = 2, cex.axis = 2)
mtext("Number of growth decreasing alleles", side = 1, cex = 2, line = 3)
mtext("Colony Growth IAA", side = 2, cex = 2, line = 2.5)
legend("topright", c("BY", "RM"), col = c("grey", "green"), pch = 19, title = "Hub allele")

#Add the two regression lines to the plot
a <- group.RM.count_nlmExp$coefficients[1,1]
b <- group.RM.count_nlmExp$coefficients[2,1]
c <- group.RM.count_nlmExp$coefficients[3,1]
f_exp <- function(x){a + b*exp(c*x) }
lines(x = seq(.9, 5.9), y = f_exp(0:5), lwd = 2, col = "blue", lty = 2)
lines(x = c(.9, 5.9), y = c(group.RM.count_lmAdd$coefficients[1], group.RM.count_lmAdd$coefficients[1] + group.RM.count_lmAdd$coefficients[2]*5), lwd = 2)
R2_exp <- 1 - var(pheno[group.RM] - f_exp(group.RM.count), na.rm = T)/var(pheno[group.RM], na.rm = T)
R2_add <- 1 - var(group.RM.count_lmAdd$residuals, na.rm = T)/var(pheno[group.RM], na.rm = T)
legend("bottomleft", c(paste("Additive model       (R2 = ", round(R2_add, digits = 3), ")", sep = ""), paste("Exponential model (R2 = ", round(R2_exp, digits = 3), ")", sep = "")), col = c("black", "blue"), title = "Regression Lines", lty = "solid", lwd = 3)



#Simultaneously analyze the mean and variance effect of all epistatic loci
require(dglm)
require(igraph)
load("../data/data_GenABEL.RData")
load("../data/networks.RData")
load("../data/phenotypes.RData")
data.traitNames <- names(pheno_raw)
nrVertices <- sum(sapply(networks.index, FUN = function(x){length(V(x))}))
networkNeighbours_dGLM <- data.frame(snpIndex = rep(NA, nrVertices), nrNeighbours = NA, meanEff = NA, varEff = NA, meanEff_p = NA, varEff_p = NA, trait = NA)

for(trait in names(networks.index)){
  #Count the number of neighbours of every loci in the network
  network <- networks.index[[trait]]
  network.vertices <- as.numeric(V(network)$name)
  network.vertices_neighbours <- sapply(1:length(network.vertices), FUN = function(x){ length(network.vertices[neighbors(graph = network, v = x, mode = "total")]) })
  
  #Put it into the data frame
  from <- min(which(is.na(networkNeighbours_dGLM$snpIndex)))
  to <- from + length(network.vertices) - 1
  networkNeighbours_dGLM$snpIndex[from:to] <- network.vertices
  networkNeighbours_dGLM$nrNeighbours[from:to] <- network.vertices_neighbours
  networkNeighbours_dGLM$trait[from:to] <- trait
  
  #fit dGLM
  data.traitNr <- grep(pattern = trait, x = data.traitNames)
  dGLM.effects <- matrix(nrow = length(network.vertices), ncol = 4)
  for(i in 1:length(network.vertices)){
    geno <- as.double.gwaa.data(data[,network.vertices[i]])
    geno[geno == 2] <- 1
    pheno <- data@phdata[,paste("trait_", data.traitNr, ".mean", sep = "")]
    dGLM <- dglm(formula = pheno~geno, dformula = pheno~geno, data = data.frame(pheno = pheno, geno = geno[,1]))
    dGLM.effects[i, 1] <- dGLM$dispersion.fit$coefficients[2] #The dispersion effect
    dGLM.effects[i, 2] <- dGLM$coefficients[2] #The mean effect
    
    dGLM.effects[i, 3] <- summary(dGLM$dispersion.fit)$coefficients["geno", 4] #The p-value of the dispersion effect
    dGLM.effects[i, 4] <- summary(dGLM)$coefficients["geno", 4] #The p-value of the mean effect
  } 
  
  #Put it into the data frame
  networkNeighbours_dGLM$varEff[from:to] <- dGLM.effects[,1]
  networkNeighbours_dGLM$meanEff[from:to] <- dGLM.effects[,2]
  networkNeighbours_dGLM$varEff_p[from:to] <- dGLM.effects[,3]
  networkNeighbours_dGLM$meanEff_p[from:to] <- dGLM.effects[,4]
}
save(list = "networkNeighbours_dGLM", file = "160513_epiLoci_dGLM.RData")
#Correlation test
cor.test(networkNeighbours_dGLM$nrNeighbours, exp(abs(networkNeighbours_dGLM$varEff))) #Taking exp(b) where b is the dispersion effect from the dglm model gives the fold difference in variance between the genotypes 



#Additive vs epistatic models fitted to the whole population
require(igraph)
load("../data/data_GenABEL.RData")
load("../data/networkInfo.RData")
load("../data/phenotypes.RData")
load("../data/networks.RData")
networkNeighbours.hubs <- networkNeighbours[networkNeighbours$nrNeighbours > 4,]
#Skip the hubs that give the completely overlapping CopperSulfate networks
networkNeighbours.hubs <- networkNeighbours.hubs[c(1:2, 7:19), ]
networkNeighbours.hubs_R2 <- cbind(networkNeighbours.hubs, add.R2 = NA, epi2.R2 = NA, epi3.R2 = NA, epiFull.R2 = NA)
for(i in 1:nrow(networkNeighbours.hubs)){
  hub <- networkNeighbours.hubs$snpIndex[i]
  trait <- networkNeighbours.hubs$trait[i]
  trait.MyName <- paste("trait_", grep(pattern = trait, x = names(pheno_raw)), ".mean", sep = "") #the trait name used in data
  network <- networks.index[[trait]]
  
  #Pick the top 5 interactions around the hub (Really messy. There should be a simpler way if you understand igraph)
  hubNeighbors <- E(network)[from(as.character(hub))] #the neighboring edges
  tmp <- order(hubNeighbors$LOD, decreasing=TRUE)[1:5] #sort according to LOD
  hubNeighbors.top <- E(network)[ as.vector(hubNeighbors)[tmp] ] #sorted edge list
  tmp2 <- get.edges(graph = network, hubNeighbors.top) #indices of the nodes in hubNeighbors.top
  hubNeighbors.names <- as.numeric(V(network)[unique(c(tmp2))]$name) #take out these nodes from network
  
  #construct data-frame
  geno <- as.double.gwaa.data(data[,hubNeighbors.names])
  geno[geno == 2] <- 1
  #create dummy variable that defines the 64 genotypic classes
  data.lm <- data.frame(pheno = data@phdata[,trait.MyName], geno)
  #data.lm <- data.cv[!is.na(data.cv$pheno), ]
  colnames(data.lm)[2:7] <- paste("snp", 1:6, sep = "")
  
  #compare model fits
  add.lm <- lm(pheno ~ snp1+snp2+snp3+snp4+snp5+snp6, data = data.lm)
  epi2.lm <- lm(pheno ~ (snp1+snp2+snp3+snp4+snp5+snp6)^2, data = data.lm)
  epi3.lm <- lm(pheno ~ (snp1+snp2+snp3+snp4+snp5+snp6)^3, data = data.lm)
  epiFull.lm <- lm(pheno ~ snp1*snp2*snp3*snp4*snp5*snp6, data = data.lm)
  
  networkNeighbours.hubs_R2$add.R2[i] <- summary(add.lm)$r.squared
  networkNeighbours.hubs_R2$epi2.R2[i] <- summary(epi2.lm)$r.squared
  networkNeighbours.hubs_R2$epi3.R2[i] <- summary(epi3.lm)$r.squared
  networkNeighbours.hubs_R2$epiFull.R2[i] <- summary(epiFull.lm)$r.squared
}
mean(networkNeighbours.hubs_R2$epiFull.R2 - networkNeighbours.hubs_R2$add.R2) #The average increase in R2 of a fully parameterized 6-way epistatic model compared to an additive model



#Cross validations
m <- 10 #number of folds
hubs_CV_add_epi <- list()
for(i in 1:nrow(networkNeighbours.hubs)){
  hub <- networkNeighbours.hubs$snpIndex[i]
  trait <- networkNeighbours.hubs$trait[i]
  trait.MyName <- paste("trait_", grep(pattern = trait, x = names(pheno_raw)), ".mean", sep = "") #the trait name used in data
  network <- networks.index[[trait]]
  
  #Pick the top 5 interactions around the hub (Really messy. There should be a simpler way if you understand igraph)
  hubNeighbors <- E(network)[from(as.character(hub))] #the neighboring edges
  tmp <- order(hubNeighbors$LOD, decreasing=TRUE)[1:5] #sort according to LOD
  hubNeighbors.top <- E(network)[ as.vector(hubNeighbors)[tmp] ] #sorted edge list
  tmp2 <- get.edges(graph = network, hubNeighbors.top) #indices of the nodes in hubNeighbors.top
  hubNeighbors.names <- as.numeric(V(network)[unique(c(tmp2))]$name) #take out these nodes from network
  
  #construct data-frame
  geno <- as.double.gwaa.data(data[,hubNeighbors.names])
  geno[geno == 2] <- 1
  #create dummy variable that defines the 64 genotypic classes
  class <- apply(X = geno, FUN = paste, MARGIN = 1, collapse = "_")
  class <- factor(class)
  data.cv <- data.frame(pheno = data@phdata[,trait.MyName], geno, class = class, cvpredAdd = NA, cvpredEpi2 = NA, cvpredEpi3 = NA)
  data.cv <- data.cv[!is.na(data.cv$pheno), ]
  colnames(data.cv)[2:7] <- paste("snp", 1:6, sep = "")
  
  #Perform the CV
  n <- dim(data.cv)[1]
  rand <- sample(n)%%m + 1
  foldnum <- sort(unique(rand))
  for (j in foldnum) {
    rows.in <- rand != j
    rows.out <- rand == j
    
    #variable selection
    model.all2way <- lm(pheno ~ (snp1+snp2+snp3+snp4+snp5+snp6)^2, data = data.cv[rows.in, ]) #model including all pairwise interactions
    keep <- step(model.all2way, direction = "backward", k = 2) #backwards elimination based on AIC
    tmp <- names(keep$coefficients)
    keep2way.formula <- formula(paste("pheno ~", paste(tmp[2:length(tmp)], collapse = "+"), collapse = ""))
    
    model.all3way <- lm(pheno ~ (snp1+snp2+snp3+snp4+snp5+snp6)^3, data = data.cv[rows.in, ]) #model including all pairwise and threeway interactions
    keep <- step(model.all3way, direction = "backward", k = 2) #backwards elimination based on AIC
    tmp <- names(keep$coefficients)
    keep3way.formula <- formula(paste("pheno ~", paste(tmp[2:length(tmp)], collapse = "+"), collapse = ""))
    
    #Fit models to training data
    subsAdd.lm <- lm(pheno ~ snp1+snp2+snp3+snp4+snp5+snp6, data = data.cv[rows.in, ])
    subsEpi2.lm <- lm(keep2way.formula, data = data.cv[rows.in, ])
    subsEpi3.lm <- lm(keep3way.formula, data = data.cv[rows.in, ])
    
    #Predict in test data
    data.cv[rows.out, "cvpredAdd"] <- predict(subsAdd.lm, newdata = data.cv[rows.out, ])
    data.cv[rows.out, "cvpredEpi2"] <- predict(subsEpi2.lm, newdata = data.cv[rows.out, ])
    data.cv[rows.out, "cvpredEpi3"] <- predict(subsEpi3.lm, newdata = data.cv[rows.out, ])
  }
  
  #wrap it up
  hubs_CV_add_epi[[i]] <- data.cv
}
save(list = "hubs_CV_add_epi", file = "../data/addVSepi_CV.RData")



#Analyze cross validation results
biasTest_tTest <- as.list(rep(NA, 15))
for(i in 1:15){
  CVs <- hubs_CV_add_epi[[i]]
  #CV residuals
  add.resid <- CVs$pheno - CVs$cvpredAdd
  epi2.resid <- CVs$pheno - CVs$cvpredEpi2
  epi3.resid <- CVs$pheno - CVs$cvpredEpi3
  
  #Split according to genotype
  add.resid_split <- split(x = add.resid, f = CVs$class)
  epi2.resid_split <- split(x = epi2.resid, f = CVs$class)
  epi3.resid_split <- split(x = epi3.resid, f = CVs$class)
  
  #Test the bias in every genotype
  add.resid_bias <- sapply(add.resid_split, function(x){ t.test(x)$p.value })
  epi2.resid_bias <- sapply(epi2.resid_split, function(x){ t.test(x)$p.value })
  epi3.resid_bias <- sapply(epi3.resid_split, function(x){ t.test(x)$p.value })
  biasTest_tTest[[i]] <- list(add = add.resid_bias, epi2 = epi2.resid_bias, epi3 = epi3.resid_bias)
}
#Count the number of biased genotypes
sapply(biasTest_tTest, function(x){sum(x$add < .05/64)})
sapply(biasTest_tTest, function(x){sum(x$epi2 < .05/64)})
sapply(biasTest_tTest, function(x){sum(x$epi3 < .05/64)})



#Figure 1
#The code generates all subfigures of figure 1
load("../data/networks.RData")
load("../data/networkInfo.RData")
for(i in 1:length(networks.index)){
  network <- networks.index[[i]]
  trait <- names(networks.index)[i]
  network_nrNeighbors <- sapply(X = V(network), FUN = function(x){length(neighbors(graph = network, v = x))})
  V(network)$label <- ""
  #   V(network)$label <- network_nrNeighbors
  #   V(network)$label.cex <- 1.6
  E(network)$width <- 3
  hubs <- names(V(network))[network_nrNeighbors > 4]
  V(network)$color <- "lightblue"
  V(network)$color[names(V(network)) %in% hubs] <- "red"
  pdf(file = paste("../figures/fig1/", trait, ".pdf", sep = ""), width = 9, height = 6)
  plot(network)
  dev.off()
}
#histogram
pdf(file = "../figures/fig1/nrIntHisto.pdf", width = 9, height = 6)
par(mar = c(5,4,2,2) + .1)
hist(networkNeighbours$nrNeighbours, main = "", xlab = "", ylab = "", cex.axis = 1.6, col = c(rep("lightblue", 4), rep("red", 8)))
abline(v = 5, lty = 2)
mtext(text = "Number of loci", side = 2, cex = 1.6, line = 2.5)
mtext(text = "Number of genetic interactions", side = 1, cex = 1.6, line = 3)
dev.off()



#Figure 2 panel a
load("../data/data_GenABEL.RData")
snps <- c("4944074_chrVIII_114144_A_G", "9716860_chrXIV_469224_A_G", "1242017_chrIII_198615_T_G", "2358650_chrIV_998628_A_T", "8733525_chrXIII_410320_T_C", "7890567_chrXII_645539_A_G")
havePheno <- data@gtdata@idnames[!is.na(data@phdata$trait_8.mean)]
geno <- as.double.gwaa.data(data[havePheno, snps])
geno[geno == 2] <- 1
pheno <- data@phdata[havePheno, "trait_8.mean"]
#Code the genotypes so that 1 is the allele with a lower phenotype. The purpose is to later on add a x-axis with the genotypes
for(i in 2:ncol(geno)){ 
  pheno.split <- split(pheno, geno[,i])
  if(mean(pheno.split$`1`, na.rm = T) > mean(pheno.split$`0`, na.rm = T)){
    geno[,i] <- geno[,i]+1
    geno[geno[,i] == 2, i] <- 0
  }
}

#Re-label the genotypes
geno[geno == 1] <- "L"
geno[geno == "0"] <- "H"
geno[geno[,1] == "L", 1] <- "RM"
geno[geno[,1] == "H", 1] <- "BY"
#Keep track of which allele is from which radial QTL
tmp <- gsub(pattern = ".*_(chr.*?)_.*", replacement = "\\1", x = colnames(geno))
geno[,2] <- paste(geno[,2], tmp[2], sep = "_")
geno[,3] <- paste(geno[,3], tmp[3], sep = "_")
geno[,4] <- paste(geno[,4], tmp[4], sep = "_")
geno[,5] <- paste(geno[,5], tmp[5], sep = "_")
geno[,6] <- paste(geno[,6], tmp[6], sep = "_")

#Sort the genotypes according to the predictions from an additive model
names(pheno) <- data@gtdata@idnames[havePheno]
iaa.lm <- lm(pheno ~ geno[,2]+geno[,3]+geno[,4]+geno[,5]+geno[,6] + geno[,1])
box <- boxplot(iaa.lm$fitted.values ~ geno[,2]+geno[,3]+geno[,4]+geno[,5]+geno[,6] + geno[,1], plot = F)
newOrder.tmp <- c(order(box$stats[3,1:32]))
newOrder <- c(newOrder.tmp, newOrder.tmp + 32)
#sort the pred
box.pred <- box
box.pred$stats <- box$stats[,newOrder]
box.pred$n <- box$n[newOrder]
box.pred$names <- box$names[newOrder]
box.pred$conf <- box$conf[,newOrder]
for(i in 1:length(newOrder)){
  box.pred$group[box$group == newOrder[i]] <- i
}
#sort the GP-map
box <- boxplot(pheno ~ geno[,2]+geno[,3]+geno[,4]+geno[,5]+geno[,6] + geno[,1], plot = F)
box.pheno <- box
box.pheno$stats <- box$stats[,newOrder]
box.pheno$n <- box$n[newOrder]
box.pheno$names <- box$names[newOrder]
box.pheno$conf <- box$conf[,newOrder]
for(i in 1:length(newOrder)){
  box.pheno$group[box$group == newOrder[i]] <- i
}

pdf(file = "../figures/fig2_a.pdf", width = 9, height = 6)
par(mar = c(4.5, 3, 3, 2) + .1)
bxp(box.pheno, ylim = c(-3.2, 1.9), xaxt = "n", boxfill = c(rep("green", 32), rep("grey", 32)), frame = F, yaxt = "n", outcex = .6)
axis(side = 2, cex.axis = 1.6, line = -1)

#Add points as x-labels to indicate the genotype in every class
axis(side = 1, at = 1:64, labels = F)
genoClass <- gsub(pattern = "(.*\\.).*", replacement = "\\1", x = box.pheno$names) #skip info about the hub-genotype
tmp <- sapply(X = genoClass, FUN = function(x){strsplit(x = x, split = "\\.") })
genoClass <- matrix(unlist(tmp), ncol = 64)
for(i in 1:5){
  col <- rep("orange", 64)
  #col[genoClass[i,] == "L"] <- "orange"
  col[grep(pattern = "L", x = genoClass[i,])] <- "blue"
  y <- -3.5 - .1*i
  points(x = 1:64, y = rep(y, 64), col = col, pch = 19, xpd = T, cex = .6)
}

axis(side = 1, at = 1:64, labels = box.pheno$n, cex.axis = .75, line = -2, tick = F, las = 2) #Write sample sizes
lines(x = 1:32, y = box.pred$stats[3, 1:32], col = "black", lwd = 2.5)
lines(x = 33:64, y = box.pred$stats[3, 33:64], col = "black", lwd = 2.5)
mtext("Colony growth on IAA containing medium", side = 2, cex = 1.6, line = 1.5, adj = 0)
mtext("Six-locus genotype-class", side = 1, cex = 1.6, line = 3)
legend(x = 44, y = -2.15, c("BY hub-QTL allele", "RM hub-QTL allele"), col = c("green", "grey"), cex = .75, pch = 15, bty = "n") #title = "hub genotype"
legend(x = 44, y = -2.6, c("Growth-increasing allele radial-QTL", "Growth-decreasing allele radial-QTL"), col = c("orange", "blue"), cex = .75, pch = 19, bty = "n")
dev.off()



#Figure 2 panel b
snps <- c("4944074_chrVIII_114144_A_G", "9716860_chrXIV_469224_A_G", "1242017_chrIII_198615_T_G", "2358650_chrIV_998628_A_T", "8733525_chrXIII_410320_T_C", "7890567_chrXII_645539_A_G")
geno <- as.double.gwaa.data(data[, snps])
geno[geno == 2] <- 1
pheno <- data@phdata[, "trait_8.mean"]
names(pheno) <- data@gtdata@idnames
#Code the genotypes so that 1 is the allele with a lower phenotype. The purpose is to later on add a x-axis with the genotypes
for(i in 2:ncol(geno)){ 
  pheno.split <- split(pheno, geno[,i])
  if(mean(pheno.split$`1`, na.rm = T) > mean(pheno.split$`0`, na.rm = T)){
    geno[,i] <- geno[,i]+1
    geno[geno[,i] == 2, i] <- 0
  }
}
group1 <- rownames(geno)[geno[,1] == 1]
group2 <- rownames(geno)[geno[,1] == 0]

#Make allele-count variables
tmp <- geno[group1, 2:6]
group1.count <- apply(X = tmp, MARGIN = 1, FUN = sum)
tmp <- geno[group2, 2:6]
group2.count <- apply(X = tmp, MARGIN = 1, FUN = sum)

#Fit additive and exponential models to the allele count
group2.count_lmAdd <- lm(pheno[group2] ~ group2.count) #additive
data.exp <- data.frame(pheno = pheno[group2], count = group2.count)
data.exp <- data.exp[!is.na(data.exp$pheno), ]
group2.count_nlmExp <- summary(nls(pheno ~ a + b*exp(c*count), start = list(a = 1, b = 1, c = 1), data = data.exp))

pdf(file = "../figures/fig2_b.pdf", width = 9, height = 6)
par(mar = c(5,2,2,2) + .1)
ylim <- c(min(pheno, na.rm = T), max(pheno, na.rm = T))
cols <- c("grey", "green")
box.group1 <- boxplot(pheno[group1] ~ group1.count, boxwex = .15, at = seq(.9, 5.9), col = cols[1], ylim = ylim, xaxt = "n", frame = F, yaxt = "n", outcex = .6)
box.group2 <- boxplot(pheno[group2] ~ group2.count, boxwex = .15, at = seq(1.1, 6.1), add = T, col = cols[2], ylim = ylim, xaxt = "n", frame = F, yaxt = "n", outcex = .6)
axis(side = 1, at = 1:6, labels = 0:5, cex.axis = 1.6, padj = .5, line = .8)
axis(side = 2, cex.axis = 1.6, line = -2)
mtext("Number of growth decreasing alleles at the radial QTL", side = 1.5, cex = 1.6, line = 3.5)
mtext("Colony growth on IAA containing medium", side = 2, cex = 1.6, line = .5)
legend(x = .6, y = -2.5, c("BY hub-QTL allele", "RM hub-QTL allele"), col = c(cols[2], cols[1]), pch = 15, cex = .75, bty = "n")
legend(x = .6, y = -2, c("Additive model fit", "Exponential model fit"), col = c("black", "blue"), lty = "solid", lwd = 3, cex = .75, bty = "n")

#Write sample sizes
axis(side = 1, at = seq(.9, 5.9), labels = box.group1$n, tick = F, line = -2, las = 2, cex.axis = 1)
axis(side = 1, at = seq(1.1, 6.1), labels = box.group2$n, tick = F, line = -2, las = 2, cex.axis = 1)

#Add the two regression lines to the plot
a <- group2.count_nlmExp$coefficients[1,1]
b <- group2.count_nlmExp$coefficients[2,1]
c <- group2.count_nlmExp$coefficients[3,1]
f_exp <- function(x){a + b*exp(c*x) }
lines(x = seq(1.1, 6.1, by = .1), y = f_exp(seq(0, 5, by = .1)), lwd = 2, col = "blue", lty = 2) #smother line
lines(x = c(1.1, 6.1), y = c(group2.count_lmAdd$coefficients[1], group2.count_lmAdd$coefficients[1] + group2.count_lmAdd$coefficients[2]*5), lwd = 2, col = "black")
dev.off()



#Figure 3
#The code generates figures like 2b for all 15 hub networks. The figures corresponding to the traits E6-berbamine, Manganese sulfate, Neomycin and Copper Sulfate are the 4 panels in fig 3
load("../data/networkInfo.RData")
load("../data/networks.RData")
load("../data/phenotypes.RData")
load("../data/genotypes.RData")
load("../data/data_GenABEL.RData")
networkNeighbours.hubs <- networkNeighbours[networkNeighbours$nrNeighbours > 4,]
#Skip the hubs that give the completely overlapping CopperSulfate networks
networkNeighbours.hubs <- networkNeighbours.hubs[c(1:2, 7:19), ]
for(i in 1:15){
  hub <- networkNeighbours.hubs$snpIndex[i]
  trait <- networkNeighbours.hubs$trait[i]
  trait.MyName <- paste("trait_", grep(pattern = trait, x = names(pheno_raw)), ".mean", sep = "") #the trait name used in data
  network <- networks.index[[trait]]
  
  #Pick the top 5 interactions around the hub (Really messy. There should be a simpler way if you understand igraph)
  hubNeighbors <- E(network)[from(as.character(hub))] #the neighboring edges
  tmp <- order(hubNeighbors$LOD, decreasing=TRUE)[1:5] #sort according to LOD
  hubNeighbors.top <- E(network)[ as.vector(hubNeighbors)[tmp] ] #sorted edge list
  tmp2 <- get.edges(graph = network, hubNeighbors.top) #indices of the nodes in hubNeighbors.top
  hubNeighbors.names <- as.numeric(V(network)[unique(c(tmp2))]$name) #take out these nodes from network
  
  geno <- as.double.gwaa.data(data[, c(hub, hubNeighbors.names[hubNeighbors.names != hub])])
  geno[geno == 2] <- 1
  pheno <- data@phdata[, trait.MyName]
  #Code the genotypes so that 1 is the allele with a lower phenotype. The purpose is to later on add a x-axis with the genotypes
  for(j in 2:ncol(geno)){ 
    pheno.split <- split(pheno, geno[,j])
    if(mean(pheno.split$`1`, na.rm = T) > mean(pheno.split$`0`, na.rm = T)){
      geno[,j] <- geno[,j]+1
      geno[geno[,j] == 2, j] <- 0
    }
  }
  
  #Check the line origin of the hub-alleles
  #Coding in gdata: BY = -1, RM = 1
  hub.count1 <- table(geno[,1])
  hub.count2 <- table(gdata[,hub])
  if(sum(hub.count1 %in% hub.count2) != 2)
    print(paste("Trubbel iteration:", i, "Trait:", trait))
  if(hub.count1[names(hub.count1) == 0] == hub.count2[names(hub.count2) == -1]){
    group.RM <- rownames(geno)[geno[,1] == 1]
    group.BY <- rownames(geno)[geno[,1] == 0]
  }else{
    group.RM <- rownames(geno)[geno[,1] == 0]
    group.BY <- rownames(geno)[geno[,1] == 1]
  }
  
  #Make allele-count variables
  tmp <- geno[group.RM, 2:6]
  group.RM.count <- apply(X = tmp, MARGIN = 1, FUN = sum)
  tmp <- geno[group.BY, 2:6]
  group.BY.count <- apply(X = tmp, MARGIN = 1, FUN = sum)
  
  #Fit additive models to the allele count
  group.RM.count_lmAdd <- lm(data@phdata[group.RM, trait.MyName] ~ group.RM.count) 
  group.BY.count_lmAdd <- lm(data@phdata[group.BY, trait.MyName] ~ group.BY.count) 
  
  pdf(file = paste("../figures/alleleCount/", trait, i, ".pdf", sep = ""), width = 9, height = 6)
  par(mar = c(4,2,2,1) + .1)
  ylim <- c(min(data@phdata[,trait.MyName], na.rm = T), max(data@phdata[,trait.MyName], na.rm = T))
  tmp <- which.max(c(max(abs(data@phdata[group.RM, trait.MyName]), na.rm = T), max(abs(data@phdata[group.BY, trait.MyName]), na.rm = T)))
  box.group.RM <- boxplot(data@phdata[group.RM, trait.MyName] ~ group.RM.count, boxwex = .15, at = seq(.9, 5.9), col = "grey", ylim = ylim, xaxt = "n", frame = F, yaxt = "n")
  box.group.BY <- boxplot(data@phdata[group.BY, trait.MyName] ~ group.BY.count, boxwex = .15, at = seq(1.1, 6.1), add = T, col = "green", ylim = ylim, xaxt = "n", frame = F, yaxt = "n")
  #title(main = )
  axis(side = 1, at = 1:6, labels = 0:5, cex.axis = 1.6, padj = .5)
  axis(side = 2, cex.axis = 1.6, line = -3)
  mtext("Number of growth decreasing alleles at the radial QTL", side = 1.5, cex = 1.3, line = 2.8)
  mtext(paste("Colony Growth on", trait, "containing medium"), side = 2, cex = 1.3, line = 0)
  
  #Add the two regression lines to the plot
  lines(x = c(.9, 5.9), y = c(group.RM.count_lmAdd$coefficients[1], group.RM.count_lmAdd$coefficients[1] + group.RM.count_lmAdd$coefficients[2]*5), lwd = 2, col = "black")
  lines(x = c(1.1, 6.1), y = c(group.BY.count_lmAdd$coefficients[1], group.BY.count_lmAdd$coefficients[1] + group.BY.count_lmAdd$coefficients[2]*5), lwd = 2, col = "black")
  dev.off()
}



#figure 4
require(gtools)
load("../data/addVSepi_CV.RData")
load("../data/networkInfo.RData")
load("../data/genotypes.RData")
load("../data/hubSplits_h2Polygenic.RData")
load("../data/hubs_capacEff.RData")
load("../data/data_GenABEL.RData")
networkNeighbours.hubs <- networkNeighbours[networkNeighbours$nrNeighbours > 4,]
#Skip the hubs that give the completely overlapping CopperSulfate networks
networkNeighbours.hubs <- networkNeighbours.hubs[c(1:2, 7:19), ]
h2Split_polygenic_hubs <- h2Split_polygenic_hubs[c(1:2, 7:19)]
cols <- rbind(c(0, 0, 0),   #the additive points (rgb)
              c(0.00392, 0.65, 0.44314), #epi2 points       0.847, 0.702, 0.396
              c(0.651, 0.380, 0.102)) #epi3 points    0.00392, 0.52157, 0.44314  olivedrab: 0.420, 0.557, 0.137

pch <- c(19, 2, 4) #The style of the points (additive, epi2, epi3)
alpha <- .5 #The transparency of all points
#Create the plotting frame that the lines will be added to
pdf(file = "../figures/figure4.pdf", width = 9, height = 6)
par(mar = c(3.5,4,4,2) + .1)
plot(1:67, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(-1.2, 1.1)) #Empty plot
lines(x = c(0,32), y = c(0,0), lty = 2)
lines(x = c(36, 67), y = c(0,0), lty = 2)
mtext(text = "a", side = 3, at = 2, cex = 2)
mtext(text = "b", side = 3, at = 37, cex = 2)
axis(side = 2, cex.axis = 1.6, line = -.5)
mtext(text = "Average Estimation Error", side = 2, line = 2.5, cex = 1.6)
axis(side = 1, at = 1:32, labels = F, line = -.8)
axis(side = 1, at = 36:67, labels = F, line = -.8)
#Dots on the x-axis
genoClass <- permutations(n = 2, r = 5, v = 0:1, repeats.allowed = T)
tmp <- apply(genoClass, MARGIN = 1, sum)
genoClass <- genoClass[order(tmp),]
for(i in 1:5){
  col <- rep("blue", 32)
  col[genoClass[,i] == 1] <- "orange"
  y <- -1.24 - .05*i
  points(x = 1:32, y = rep(y, 32), col = col, pch = 19, xpd = T, cex = .6)
  points(x = 36:67, y = rep(y, 32), col = col, pch = 19, xpd = T, cex = .6)
}
mtext("Six-locus genotype-class", side = 1, cex = 1.6, line = 2.5)

legend("topright", c("Additive", "Epistatic"), col = c(rgb(red = cols[1,1], green = cols[1,2], blue = cols[1,3]), rgb(red = cols[2,1], green = cols[2,2], blue = cols[2,3])), title = expression(bold("Model")), pch = pch, cex = .75, bty = "n")
legend(x = 48, y = -.9, c("Growth-increasing allele radial-QTL", "Growth-decreasing allele radial-QTL"), col = c("orange", "blue"), cex = .75, pch = 19, bty = "n")
capac.misspred_add <- matrix(NA, nrow = 10, ncol = 32)
canal.misspred_add <- matrix(NA, nrow = 10, ncol = 32)
capac.misspred_epi2 <- matrix(NA, nrow = 10, ncol = 32)
canal.misspred_epi2 <- matrix(NA, nrow = 10, ncol = 32)
capac.misspred_epi3 <- matrix(NA, nrow = 10, ncol = 32)
canal.misspred_epi3 <- matrix(NA, nrow = 10, ncol = 32)

j <- 1
for(i in 1:length(hubs_CV_add_epi)){
  pheno <- hubs_CV_add_epi[[i]]$pheno
  names(pheno) <- rownames(hubs_CV_add_epi[[i]])
  pheno.norm <- (pheno - mean(pheno))/sd(pheno)
  misspred_add <- hubs_CV_add_epi[[i]]$pheno - hubs_CV_add_epi[[i]]$cvpredAdd
  names(misspred_add) <- rownames(hubs_CV_add_epi[[i]])
  misspred_epi2 <- hubs_CV_add_epi[[i]]$pheno - hubs_CV_add_epi[[i]]$cvpredEpi2
  names(misspred_epi2) <- rownames(hubs_CV_add_epi[[i]])
  misspred_epi3 <- hubs_CV_add_epi[[i]]$pheno - hubs_CV_add_epi[[i]]$cvpredEpi3
  names(misspred_epi3) <- rownames(hubs_CV_add_epi[[i]])
  trait <- networkNeighbours.hubs$trait[i]
  hub <- networkNeighbours.hubs$snpName[i]
  
  #Check which hub-genotype is the capacitated one
  #The structure of h2Split_polygenic_hubs (which contains h2s in the two hub-genotypes, see code line 121 - 162):
  #  h2Split_polygenic_hubs[[i]]
  #     [[1]] h2 of the hub-genotype coded as 0
  #     [[2]] h2 of the hub-genotype coded as 2
  hub.geno <- as.double.gwaa.data(data[,hub])
  hub.h2s <- h2Split_polygenic_hubs[[i]]
  if(hub.h2s[[1]] > hub.h2s[[2]]){ #So if this is true, the guys with geno == 0 are the capacitated ones
    group.capac <- rownames(hub.geno)[hub.geno[,1] == 0]
    group.canal <- rownames(hub.geno)[hub.geno[,1] == 2]
  }else{
    group.capac <- rownames(hub.geno)[hub.geno[,1] == 2]
    group.canal <- rownames(hub.geno)[hub.geno[,1] == 0]
  }
  
  #Pick the radial loci
  geno <- hubs_CV_add_epi[[i]][,2:7]
  tmp <- hub.geno[rownames(geno),] #only the segregants with a phenotype
  tmp[tmp == 2] <- 1
  radial <- which(!sapply(1:6, function(x){all(geno[,x] == tmp)}))
  if(length(radial) != 5){print(paste("Trubbel i =", i, "hub:", hub))} #Security check
  geno.radial <- geno[,radial]
  
  #Code the genotypes at the radial loci so that 1 is the allele with the lower phenotype
  for(k in 1:ncol(geno.radial)){
    pheno.split <- split(pheno, geno.radial[,k])
    if(mean(pheno.split$`1`, na.rm = T) > mean(pheno.split$`0`, na.rm = T)){
      geno.radial[,k] <- geno.radial[,k]+1
      geno.radial[geno.radial[,k] == 2, k] <- 0
    }
  }

  #Resolve ties in allelecount based on the effect sizes
  geno.radial.lm <- lm(pheno ~ geno.radial[,1] + geno.radial[,2] + geno.radial[,3] + geno.radial[,4] + geno.radial[,5])
  allelesRank <- abs(geno.radial.lm$coefficients[2:6])/100
  geno.radial[geno.radial[,1] == 1,1] <- geno.radial[geno.radial[,1] == 1,1] + allelesRank[1]
  geno.radial[geno.radial[,2] == 1,2] <- geno.radial[geno.radial[,2] == 1,2] + allelesRank[2]
  geno.radial[geno.radial[,3] == 1,3] <- geno.radial[geno.radial[,3] == 1,3] + allelesRank[3]
  geno.radial[geno.radial[,4] == 1,4] <- geno.radial[geno.radial[,4] == 1,4] + allelesRank[4]
  geno.radial[geno.radial[,5] == 1,5] <- geno.radial[geno.radial[,5] == 1,5] + allelesRank[5]
  
  #The capacitated group
  geno.radial.capac <- na.omit(geno.radial[group.capac, ]) #Genotypes at the radial loci, capacitated group
  geno.radial.capac_genoClasses <- apply(X = geno.radial.capac, MARGIN = 1, FUN = function(x){paste(x, collapse = "_") })
  #Additive
  misspred_add.capac <- na.omit(misspred_add[group.capac]) #Additive misspredictions
  if(!identical(names(misspred_add.capac), names(geno.radial.capac_genoClasses))){print(paste("Trubbel: The additive predictions and the genotypes dont add up. i =", i))} #Security check
  misspred_add.capac_perClass <- split(x = misspred_add.capac, f = geno.radial.capac_genoClasses) #Additive misspredictions per genotype class
  misspred_add.capac_perClassMeans <- sapply(misspred_add.capac_perClass, mean) #Average missprediction per genotype class
  #Sort the classes according to the number of +/- radial alleles
  alleleCount <- sapply(names(misspred_add.capac_perClassMeans), FUN = function(x){ sum(as.numeric(strsplit(x, split = "_")[[1]] )) })
  misspred_add.capac_perClassMeans <- misspred_add.capac_perClassMeans[order(alleleCount, decreasing = T)] 
  #Epi2
  misspred_epi2.capac <- na.omit(misspred_epi2[group.capac]) #Epi2 misspredictions
  if(!identical(names(misspred_epi2.capac), names(geno.radial.capac_genoClasses))){print(paste("Trubbel: The epi2 predictions and the genotypes dont add up. i =", i))} #Security check
  misspred_epi2.capac_perClass <- split(x = misspred_epi2.capac, f = geno.radial.capac_genoClasses) #Additive misspredictions per genotype class
  misspred_epi2.capac_perClassMeans <- sapply(misspred_epi2.capac_perClass, mean) #Average missprediction per genotype class
  #Sort the classes according to the number of +/- radial alleles
  alleleCount <- sapply(names(misspred_epi2.capac_perClassMeans), FUN = function(x){ sum(as.numeric(strsplit(x, split = "_")[[1]] )) })
  misspred_epi2.capac_perClassMeans <- misspred_epi2.capac_perClassMeans[order(alleleCount, decreasing = T)] 
  #Epi3
  misspred_epi3.capac <- na.omit(misspred_epi3[group.capac]) #Epi3 misspredictions
  if(!identical(names(misspred_epi3.capac), names(geno.radial.capac_genoClasses))){print(paste("Trubbel: The epi3 predictions and the genotypes dont add up. i =", i))} #Security check
  misspred_epi3.capac_perClass <- split(x = misspred_epi3.capac, f = geno.radial.capac_genoClasses) #Additive misspredictions per genotype class
  misspred_epi3.capac_perClassMeans <- sapply(misspred_epi3.capac_perClass, mean) #Average missprediction per genotype class
  #Sort the classes according to the number of +/- radial alleles
  alleleCount <- sapply(names(misspred_epi3.capac_perClassMeans), FUN = function(x){ sum(as.numeric(strsplit(x, split = "_")[[1]] )) })
  misspred_epi3.capac_perClassMeans <- misspred_epi3.capac_perClassMeans[order(alleleCount, decreasing = T)] 
  
  #The canalized group
  geno.radial.canal <- na.omit(geno.radial[group.canal, ]) #Genotypes at the radial loci, canalized group
  geno.radial.canal_genoClasses <- apply(X = geno.radial.canal, MARGIN = 1, FUN = function(x){paste(x, collapse = "_") })
  #Additive
  misspred_add.canal <- na.omit(misspred_add[group.canal]) #Additive misspredictions
  if(!identical(names(misspred_add.canal), names(geno.radial.canal_genoClasses))){print(paste("Trubbel: The additive predictions and the genotypes dont add up. i =", i))} #Security check
  misspred_add.canal_perClass <- split(x = misspred_add.canal, f = geno.radial.canal_genoClasses) #Additive misspredictions per genotype class
  misspred_add.canal_perClassMeans <- sapply(misspred_add.canal_perClass, mean) #Average missprediction per genotype class
  #Sort the classes according to the number of +/- radial alleles
  alleleCount <- sapply(names(misspred_add.canal_perClassMeans), FUN = function(x){ sum(as.numeric(strsplit(x, split = "_")[[1]] )) })
  misspred_add.canal_perClassMeans <- misspred_add.canal_perClassMeans[order(alleleCount, decreasing = T)] 
  #Epi2
  misspred_epi2.canal <- na.omit(misspred_epi2[group.canal]) #Epi2 misspredictions
  if(!identical(names(misspred_epi2.canal), names(geno.radial.canal_genoClasses))){print(paste("Trubbel: The epi2 predictions and the genotypes dont add up. i =", i))} #Security check
  misspred_epi2.canal_perClass <- split(x = misspred_epi2.canal, f = geno.radial.canal_genoClasses) #Additive misspredictions per genotype class
  misspred_epi2.canal_perClassMeans <- sapply(misspred_epi2.canal_perClass, mean) #Average missprediction per genotype class
  #Sort the classes according to the number of +/- radial alleles
  alleleCount <- sapply(names(misspred_epi2.canal_perClassMeans), FUN = function(x){ sum(as.numeric(strsplit(x, split = "_")[[1]] )) })
  misspred_epi2.canal_perClassMeans <- misspred_epi2.canal_perClassMeans[order(alleleCount, decreasing = T)] 
  #Epi3
  misspred_epi3.canal <- na.omit(misspred_epi3[group.canal]) #Epi2 misspredictions
  if(!identical(names(misspred_epi3.canal), names(geno.radial.canal_genoClasses))){print(paste("Trubbel: The epi3 predictions and the genotypes dont add up. i =", i))} #Security check
  misspred_epi3.canal_perClass <- split(x = misspred_epi3.canal, f = geno.radial.canal_genoClasses) #Additive misspredictions per genotype class
  misspred_epi3.canal_perClassMeans <- sapply(misspred_epi3.canal_perClass, mean) #Average missprediction per genotype class
  #Sort the classes according to the number of +/- radial alleles
  alleleCount <- sapply(names(misspred_epi3.canal_perClassMeans), FUN = function(x){ sum(as.numeric(strsplit(x, split = "_")[[1]] )) })
  misspred_epi3.canal_perClassMeans <- misspred_epi3.canal_perClassMeans[order(alleleCount, decreasing = T)] 
  
  #Plot
  if(hubs.summary$h2Diff_p[i] < .05/15){ #plot either the capacitated (<) or non-capacitated (>) networks
    capac.misspred_add[j,] <- misspred_add.capac_perClassMeans
    canal.misspred_add[j,] <- misspred_add.canal_perClassMeans
    capac.misspred_epi2[j,] <- misspred_epi2.capac_perClassMeans
    canal.misspred_epi2[j,] <- misspred_epi2.canal_perClassMeans
    capac.misspred_epi3[j,] <- misspred_epi3.capac_perClassMeans
    canal.misspred_epi3[j,] <- misspred_epi3.canal_perClassMeans
    
    points(1:32, misspred_add.capac_perClassMeans, lwd = 1, col = rgb(red = cols[1,1], green = cols[1,2], blue = cols[1,3], alpha = alpha), pch = pch[1], cex = .4)
    points(36:67, misspred_add.canal_perClassMeans, lwd = 1, col = rgb(red = cols[1,1], green = cols[1,2], blue = cols[1,3], alpha = alpha), pch = pch[1], cex = .4)
    points(1:32, misspred_epi2.capac_perClassMeans, lwd = 1, col = rgb(red = cols[2,1], green = cols[2,2], blue = cols[2,3], alpha = alpha), pch = pch[2], cex = .4)
    points(36:67, misspred_epi2.canal_perClassMeans, lwd = 1, col = rgb(red = cols[2,1], green = cols[2,2], blue = cols[2,3], alpha = alpha), pch = pch[2], cex = .4)

    j<-j+1
  }
}
#Add regression lines
x <- c(sapply(X = 1:32, FUN = function(x){rep(x, times = 10)}))
lwd <- 3
lty <- c(1,2,5)
#Capac.
#Additive
capac.misspred_add_vec <- c(capac.misspred_add)
capac.misspred_add_lm <- lm(capac.misspred_add_vec ~ x)
regY <- c(capac.misspred_add_lm$coefficients[1] + capac.misspred_add_lm$coefficients[2], capac.misspred_add_lm$coefficients[1] + capac.misspred_add_lm$coefficients[2]*32)
lines(x = c(1, 32), y = regY, lwd = lwd, col = rgb(red = cols[1,1], green = cols[1,2], blue = cols[1,3]), lty = lty[1])
#Epi2
capac.misspred_epi2_vec <- c(capac.misspred_epi2)
capac.misspred_epi2_lm <- lm(capac.misspred_epi2_vec ~ x)
regY <- c(capac.misspred_epi2_lm$coefficients[1] + capac.misspred_epi2_lm$coefficients[2], capac.misspred_epi2_lm$coefficients[1] + capac.misspred_epi2_lm$coefficients[2]*32)
#lines(x = c(1, 32), y = regY, lwd = lwd, col = rgb(red = cols[2,1], green = cols[2,2], blue = cols[2,3]), lty = lty[2])
lines(x = c(1, 32), y = regY, lwd = lwd, col = "black", lty = lty[2])

#Canal.
#Additive
canal.misspred_add_vec <- c(canal.misspred_add)
canal.misspred_add_lm <- lm(canal.misspred_add_vec ~ x)
regY <- c(canal.misspred_add_lm$coefficients[1] + canal.misspred_add_lm$coefficients[2], canal.misspred_add_lm$coefficients[1] + canal.misspred_add_lm$coefficients[2]*32)
lines(x = c(36, 67), y = regY, lwd = lwd, col = rgb(red = cols[1,1], green = cols[1,2], blue = cols[1,3]), lty = lty[1])
#Epi2
canal.misspred_epi2_vec <- c(canal.misspred_epi2)
canal.misspred_epi2_lm <- lm(canal.misspred_epi2_vec ~ x)
regY <- c(canal.misspred_epi2_lm$coefficients[1] + canal.misspred_epi2_lm$coefficients[2], canal.misspred_epi2_lm$coefficients[1] + canal.misspred_epi2_lm$coefficients[2]*32)
lines(x = c(36, 67), y = regY, lwd = lwd, col = "black", lty = lty[2])
dev.off()



#figure 5
load("../data/simIAA.RData")
simR2_epi <- as.list(1:6)
simR2_add <- as.list(1:6)
for(i in 1:6){
  simR2_epi.tmp <- simIAA[[i]]$R2
  simR2_add.tmp <- simIAA_add[[i]]$R2
  
  tmp <- names(simR2_epi.tmp)
  tmp <- do.call(rbind, strsplit(x = tmp, split = "_"))
  simR2_epi[[i]] <- simR2_epi.tmp[tmp[, i] == .5 & simIAA[[i]]$p < 1e-4]
  
  tmp <- names(simR2_add.tmp)
  tmp <- do.call(rbind, strsplit(x = tmp, split = "_"))
  simR2_add[[i]] <- simR2_add.tmp[tmp[, i] == .5 & simIAA_add[[i]]$p < 1e-4]
}
labels <- c("chrVIII:114,144", "chrXIV:469,224", "chrIII:198,615", "chrIV:998,628", "chrXIII:410,320", "chrXII:645,539")
pdf(file = "../figures/figure5.pdf", width = 9, height = 6)
par(mar = c(4,4,2,1)+.1)
box1 <- boxplot(simR2_epi[[1]], simR2_epi[[2]], simR2_epi[[3]], simR2_epi[[4]], simR2_epi[[5]], simR2_epi[[6]], boxwex = .15, at = seq(.9, 5.9), xaxt = "n", col = "lightblue", frame = F, yaxt = "n")
box2 <- boxplot(simR2_add[[1]], simR2_add[[2]], simR2_add[[3]], simR2_add[[4]], simR2_add[[5]], simR2_add[[6]], boxwex = .15, at = seq(1.1, 6.1), xaxt = "n", add = T, col = "red", frame = F, yaxt = "n", border = "firebrick")
polygon(x = c(.83, .83, .97, .97), y = c(box1$stats[4,1], box1$stats[2,1], box1$stats[2,1], box1$stats[4,1]), density = 9, border = NA)
polygon(x = c(.83, .83, .97, .97) + .2, y = c(box2$stats[4,1], box2$stats[2,1], box2$stats[2,1], box2$stats[4,1]), density = 9, border = NA)
axis(1, at = 1:6, labels = rep("", 6))
text(x = 1:6, y = -.075, labels = labels, srt = 25, xpd = NA, cex = 1.1, font = 2)
axis(side = 2, cex.axis = 1.6)
mtext(text = "Variance Explained", side = 2, cex = 2, line = 2.6)
legend("topright", c("Observed Genotype-values", "Genotype-values expected under additivity"), col = c("lightblue", "red"), pch = 15, title = expression(bold("Simulations based on")), bty = "n")
dev.off()