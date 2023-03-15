install.packages("phytools")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12") 
BiocManager::install(c("EBImage", "ggtree", "Biostrings"))
install.packages("ggplot2")
install.packages("ggimage")

# Load packages
library("phytools") 
library("EBImage") # for images
library("ggtree")
library("Biostrings")
library("ggplot2")
library("ggimage")

# Import the primate phylogeny from your working directory
tree<-read.nexus("primate_tree.nex")

# Basic plot
p1 <- ggtree(tree)
plot(p1)

# Show different layouts in multiplot
p2a <- ggtree(tree, layout="rectangular") + ggtitle("rectangular")
p2b <- ggtree(tree, layout="slanted") + ggtitle("slanted")
p2c <- ggtree(tree, layout="circular") + ggtitle("circular")
p2d <- ggtree(tree, layout="radial") + ggtitle("radial")
p2e <- ggtree(tree, layout="unrooted") + ggtitle("unrooted")

multiplot(p2a, p2b, p2c, p2d, p2e, ncol=3)

p3 <- ggtree(tree, aes(color=branch.length)) +
  xlim(0, 90) + 
  geom_tiplab(size=2, color="plum1") +
  geom_label2(aes(subset=!isTip, label=node), size=2, color="darkred", alpha=0.5) +
  theme_tree2("white") +
  scale_color_continuous(low='plum1', high='plum1', name="Branch length (my)") +
  theme(legend.position="bottom")

plot(p3)

# Highlight clades (feel free to play around with your color of interest)
p4 <- ggtree(tree) + 
  xlim(0, 125) +
  geom_tiplab(size=2, offset=0.5) +
  geom_hilight(node=124, fill="steelblue", alpha=0.5) +
  geom_hilight(node=113, fill="darkgreen", alpha=0.5) +
  geom_hilight(node=110, fill="gray", alpha=0.5) +
  geom_hilight(node=96, fill="pink", alpha=0.5) +
  geom_hilight(node=89, fill="beige", alpha=0.5) +
  geom_hilight(node=70, fill="yellow", alpha=0.5) 

plot(p4)

# label clades
p5 <- p4 + 
  geom_cladelabel(124, "Galagoidea", offset=25, barsize=2, angle=90, offset.text=1.5, hjust=0.5, fontsize=3) + 
  geom_cladelabel(113, "Lemuroidea", offset=25, barsize=2, angle=90, offset.text=1.5, hjust=0.5, fontsize=3) +
  geom_cladelabel(110, "Tarsioidea", offset=25, barsize=2, angle=75, offset.text=2.5, hjust=0.2, fontsize=2) +
  geom_cladelabel(96, "Ceboidea", offset=25, barsize=2, angle=90, offset.text=1.5, hjust=0.5, fontsize=3) +
  geom_cladelabel(89, "Hominoidea", offset=25, barsize=2, angle=90, offset.text=1.5, hjust=0.5, fontsize=3) +
  geom_cladelabel(70, "Cercopithecoidea", offset=25, barsize=2, angle=90, offset.text=1.5, hjust=0.5, fontsize=3)

plot(p5)

#Now add the images to the plot
dt <- data.frame(node = c(124,113,110,96,89,70), 
                 image = c("7fb9bea8-e758-4986-afb2-95a2c3bf983d",
                           "bac25f49-97a4-4aec-beb6-f542158ebd23",
                           "f598fb39-facf-43ea-a576-1861304b2fe4",
                           "aceb287d-84cf-46f1-868c-4797c4ac54a8",
                           "0174801d-15a6-4668-bfe0-4c421fbe51e8",
                           "72f2f854-f3cd-4666-887c-35d5c256ab0f"), 
                 name = c("a","b","c","d","e","f"))

p6 <- p5 + geom_cladelab(data = dt, 
                         mapping = aes(node = node, label = name, image = image), 
                         geom = "phylopic", imagecolor = "black", 
                         offset=24, offset.text=20)

plot(p6)

# First, let's simulate tree & data
tree<-pbtree(n=26,scale=100)
tree$tip.label<-LETTERS[26:1]

lat<-fastBM(tree,sig2=10,bounds=c(-90,90))
long<-fastBM(tree,sig2=80,bounds=c(-180,180))

# Now plot projection
xx<-phylo.to.map(tree,cbind(lat,long),plot=FALSE)

## Now try to add the tree points to map coordinates
plot(xx,type="phylogram",asp=1.3,mar=c(0.1,0.5,3.1,0.1))

# A more chaotic plot, overlapping the tree projected directly on map
plot(xx,type="direct",asp=1.3,mar=c(0.1,0.5,3.1,0.1))

library(ape)
yourdata <- chiroptera.data
#Then plot it
gzoom(yourdata, grep("Plecotus", yourdata$tip.label))

library(ggtree)

set.seed(123)
tr<-rtree(60)
p<-ggtree(tr, size=2)
print(p)

cols<-rainbow(5)

for (i in 1:5) {
  p<- p + geom_balance(identify(p), fill=cols[1])
  print(p)
}

for (i in 1:5) {
  p<- p + geom_highlight(identify(p), fill=cols[1])
  print(p)
}

for (i in 1:5) {
  p %<>% rotate(identify(p))
  print(p)
}


#### EXERCISES

#Exercise 1 - Visualizing tree with multiple sequence alignment 
#With msaplot function, you can visualize multiple sequence alignment with phylogenetic tree.
#Using the FluA_H3_AA example file from the system, plot the beast_tree using the fasta argument



