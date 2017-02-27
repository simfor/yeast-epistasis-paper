#Simulate populations based on the epitatic networks sourrounding the 15 hub QTLs, as described in online methods 
#The actual simulation is performed by the function marginalEffSim_sixLoci(), and the code below is a wrapper to run th simulations in parallel
load("../data/networkInfo.RData")

networkNeighbours.hubs <- networkNeighbours[networkNeighbours$nrNeighbours > 4,]
#Skip the hubs that give the completely overlapping CopperSulfate networks
networkNeighbours.hubs <- networkNeighbours.hubs[c(1:2, 7:19), ]
library(doSNOW)
library(foreach)
cl <- makeCluster(15) 
registerDoSNOW(cl)
ptm <- proc.time()
simulations_res7 <- foreach(i = 1:15) %dopar% {
  require(GenABEL)
  require(igraph)
  source("../src/functions.R")
  load("../data/data_GenABEL.RData")
  load("../data/phenotypes.RData")
  load("../data/networkInfo.RData")
  load("../data/networks.RData")
  networkNeighbours.hubs <- networkNeighbours[networkNeighbours$nrNeighbours > 4,]
  #Skip the hubs that give the completely overlapping CopperSulfate networks
  networkNeighbours.hubs <- networkNeighbours.hubs[c(1:2, 7:19), ]
  freqs <- seq(from = .05, to = .95, length.out = 7) #simulate at these allele frequencies for all components in the networks
  
  hub <- networkNeighbours.hubs$snpIndex[i]
  trait <- networkNeighbours.hubs$trait[i]
  trait.myName <- paste("trait_", grep(pattern = trait, x = names(pheno_raw)), ".mean", sep = "") #the trait name used in data
  network <- networks.index[[trait]]
  
  #Pick the top 5 interactions around the hub (Really messy. There should be a simpler way if you understand igraph)
  hubNeighbors <- E(network)[from(as.character(hub))] #the neighboring edges
  tmp <- order(hubNeighbors$LOD, decreasing=TRUE)[1:5] #sort according to LOD
  hubNeighbors.top <- E(network)[ as.vector(hubNeighbors)[tmp] ] #sorted edge list
  tmp2 <- get.edges(graph = network, hubNeighbors.top) #indices of the nodes in hubNeighbors.top
  hubNeighbors.names <- as.numeric(V(network)[unique(c(tmp2))]$name) #take out these nodes from network
  
  simulation <- marginalEffSim_sixLoci(data = data, snps = hubNeighbors.names, freqs = freqs, trait = trait.myName)
  return(simulation)
}
stopCluster(cl)
proc.time() - ptm
save(list = c("simulations161026_res7", "networkNeighbours.hubs"), file = "../data/simulations_realGPmaps.RData")




#Simulate populations where the phenotype of every multilocus genotype is taken to be the additive expectation. See online methods for details
#The actual simulation is performed by the function marginalEffSim_sixLoci_add(), and the code below is a wrapper to run th simulations in parallel
library(doSNOW)
library(foreach)
load("../data/networkInfo.RData")
networkNeighbours.hubs <- networkNeighbours[networkNeighbours$nrNeighbours > 4,]
#Skip the hubs that give the completely overlapping CopperSulfate networks
networkNeighbours.hubs <- networkNeighbours.hubs[c(1:2, 7:19), ]
cl <- makeCluster(16)
registerDoSNOW(cl)
ptm <- proc.time()
simulationsADD <- foreach(i = 1:15) %dopar% {
  require(GenABEL)
  require(igraph)
  source("../src/functions.R")
  load("../data/data_GenABEL.RData")
  load("../data/phenotypes.RData")
  load("../data/networkInfo.RData")
  load("../data/networks.RData")
  networkNeighbours.hubs <- networkNeighbours[networkNeighbours$nrNeighbours > 4,]
  #Skip the hubs that give the completely overlapping CopperSulfate networks
  networkNeighbours.hubs <- networkNeighbours.hubs[c(1:2, 7:19), ]
  freqs <- seq(from = .05, to = .95, length.out = 7) #simulate at these allele frequencies for all components in the networks
  
  hub <- networkNeighbours.hubs$snpIndex[i]
  trait <- networkNeighbours.hubs$trait[i]
  trait.myName <- paste("trait_", grep(pattern = trait, x = names(pheno_raw)), ".mean", sep = "") #the trait name used in data
  network <- networks.index[[trait]]
  
  #Pick the top 5 interactions around the hub (Really messy. There should be a simpler way if you understand igraph)
  hubNeighbors <- E(network)[from(as.character(hub))] #the neighboring edges
  tmp <- order(hubNeighbors$LOD, decreasing=TRUE)[1:5] #sort according to LOD
  hubNeighbors.top <- E(network)[ as.vector(hubNeighbors)[tmp] ] #sorted edge list
  tmp2 <- get.edges(graph = network, hubNeighbors.top) #indices of the nodes in hubNeighbors.top
  hubNeighbors.names <- as.numeric(V(network)[unique(c(tmp2))]$name) #take out these nodes from network
  
  simulation <- marginalEffSim_sixLoci_add(data = data, snps = hubNeighbors.names, freqs = freqs, trait = trait.myName)
  return(simulation)
}
stopCluster(cl)
proc.time() - ptm
save(list = c("simulationsADD", "networkNeighbours.hubs"), file = "../data/simulationsADD.RData")



#Simulations based on the the IAA network
source("functions.R")
load("../data/data_GenABEL.RData")
snps <- c("4944074_chrVIII_114144_A_G", "9716860_chrXIV_469224_A_G", "1242017_chrIII_198615_T_G", "2358650_chrIV_998628_A_T", "8733525_chrXIII_410320_T_C", "7890567_chrXII_645539_A_G")
simIAA <- marginalEffSim_sixLoci(data = data, snps = snps, freqs = seq(from = .05, to = .95, length.out = 7), trait = "trait_8.mean")
#Simulate populations where the phenotype of every multilocus genotype is taken to be the additive expectation. See online methods for details
simIAA_add <- marginalEffSim_sixLoci_add(data = data, snps = snps, freqs = seq(from = .05, to = .95, length.out = 7), trait = "trait_8.mean")
save(list = c("simIAA", "simIAA_add"), file = "../data/simIAA.RData")


