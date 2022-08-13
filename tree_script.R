rm(list = ls())
#install packages
library(tidyverse)
library(treeio)
library(ggtree)
library(ggplot2)
library(reshape2)
library(ggstance)
library(adegenet)
library(StAMPP)
library(vcfR)
library(ggplot2)
library(MASS)
library(adegraphics)

install.packages("ggstance")

#loading in tree
tree <- read.tree("###insert file path to tree in txt or nwk format")
#create first relationships - blue = genious 
plot(tree)
ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()
p <- ggtree(tree, colour="blue", size=0.4, ladderize = F )
p + geom_tiplab(size=2, colour='black', ) 
title(expression("all individuals"))
#create branch supports
ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()
p <- ggtree(tree, colour="blue", size=0.5, ladderize = F ) + geom_text(aes(label=node), hjust=1, size=2)
p + geom_tiplab(size=2.5, colour='black', hjust=0)
plot(p)
########



