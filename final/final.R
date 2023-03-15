 install.packages("phytools")
 if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 BiocManager::install(version = "3.12") 
 BiocManager::install(c("EBImage", "ggtree", "Biostrings"))
 install.packages("ggplot2")
 install.packages("ggimage")
 install.packages("ape")
 
# load packages
library(phytools) 
library(EBImage) # for images
library(Biostrings)
library(ggtree)
library(ggplot2)
library(ggimage)
library(ggplot2)
library(ape)
library(treeio)

# import phylogeny
tree <- read.tree("mydata.tree")

# plot MAFFT output guide tree
p1 <- ggtree(tree) +   geom_tiplab(size = 2, fontface = 4)
plot(p1)

####################################

# plot BI consensus tree
BItree <-read.nexus("mrbayes_final.tre")
is.rooted(BItree) # yes

  # better tip labels for tree
BItree$tip.label

labelKeyBI <- read.csv("labelKeyBI.csv")
labelKeyBI

genus <- labelKeyBI$Genus
species <- labelKeyBI$Species
com <- labelKeyBI$Com
d <- data.frame(label = BItree$tip.label, genus = genus,
                species = species, com = com)

cols <- labelKeyBI$Color

BIplot <- ggtree(BItree, ladderize = TRUE) %<+% d + xlim(NA, 2) +
  geom_tiplab(aes(label = paste0('italic(', genus, ')~bolditalic(', species, ')~', com)), parse=T, size=2) +
  geom_tippoint(pch=16, size=2, color=as.factor(cols)) +
  geom_nodelab(geom = "label", node = "internal", size = 2, nudge_x = -0.01, alpha = 0.5) +
  ggtitle("Flourescent Protein Bayesian Inference Phylogeny") 
  
plot(BIplot)

####################################

# plot ML output tree
MLtree <- read.nexus("maxlikelihood_final.tre")
is.rooted(MLtree) # yes

MLtree$tip.label

labelKeyML <- read.csv("labelKeyML.csv")
labelKeyML

genus <- labelKeyML$Genus
species <- labelKeyML$Species
com <- labelKeyML$Com
d <- data.frame(label = labelKeyML$tip.label, genus = genus,
                species = species, com = com)

cols <- labelKeyML$Color

MLplot <- ggtree(MLtree, ladderize = TRUE) %<+% d + xlim(NA, 2) +
  geom_tiplab(aes(label = paste0('italic(', genus, ')~bolditalic(', species, ')~', com)), parse=T, size=2) +
  geom_tippoint(pch=16, size=2, color=as.factor(cols)) +
  geom_nodelab(geom = "label", node = "internal", size = 2, nudge_x = -0.01, alpha = 0.5) +
  ggtitle("Flourescent Protein Maximum Likelihood Phylogeny") 
  
plot(MLplot)

####################################

# Comparative methods

## fluorescent protein color - 5 states
# 0 - green fluorescent protein
# 1 - cyan fluorescent protein
# 2 - yellow fluorescent protein
# 3 - red fluorescent protein
# 4 - orange fluorescent protein

MLtree$tip.label
states <- factor(c(0,0,0,0,0,0,0,
                   0,0,4,0,0,1,2,
                   3,2,0,3,0,0,0,
                   0,0,0,0,3,0,0,
                   0,1,1,0,3,3,0,
                   0))

names(states) <- MLtree$tip.label

d2 = dplyr::mutate(d, newlab = paste(genus, species, sep=' '))
d3 = dplyr::mutate(d2, newlab = paste(newlab, com, sep='|'))
tr2 = rename_taxa(MLtree, d3, label, newlab)
write.tree(tr2)

# ER reconstructions
ERreconstruction <- ace(states, MLtree, type = "d", method = "ML", CI = TRUE, model = "ER")
ERreconstruction

# ER plot
plot(tr2, font = 3, lwd = 3,
     type = "p", FALSE, label.offset = 1.6, 
     main = "Flourescent Protein Color State Equal Rate Reconstruction")
co <- c("GREEN", "CYAN", "YELLOW", "RED", "ORANGE")
tiplabels(pch = 22, bg = co[as.numeric(states)], cex = 1.5, adj = 1)
nodelabels(thermo = ERreconstruction$lik.anc, piecol = co, cex = 0.3) # log likelihood of what the prob is for the ancestor color is

# AR reconstruction
ARDreconstruction <- ace(states, MLtree, type = "d", method = "ML", CI = TRUE, model = "ARD")
ARDreconstruction

# ARD plot
plot(tr2, font = 4, lwd = 6,
     type = "p", FALSE, label.offset = 1.6, 
     main = "Flourescent Protein Color State All Rates Different Reconstruction")
co <- c("GREEN", "CYAN", "YELLOW", "RED", "ORANGE")
tiplabels(aes(label = paste0('italic(', genus, ')~bolditalic(', species, ')~', com)), parse=T, size=2)
tiplabels(pch = 22, bg = co[as.numeric(states)], cex = 1.5, adj = 1)
nodelabels(thermo = ARDreconstruction$lik.anc, piecol = co, cex = 0.3) 

    
