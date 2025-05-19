#!/usr/bin/env Rscript

library(ggtree)
library(ape)
library(dplyr)

# load tree
tree.file <- '../data/processed_burkholderiales_tree.txt'
tree <- read.tree(tree.file)

# format leaf labels
tree$tip.label <- gsub("_", " ", tree$tip.label)

# get genus names
genus <- sapply(strsplit(tree$tip.label, " "), `[`, 1)
genus.df <- data.frame(tip.label = tree$tip.label, genus = genus)

# define color map
burk.colormap <- c(
  Massilia = rgb(0.193374, 0.018354, 0.59033),
  Burkholderia = rgb(0.299855, 0.009561, 0.631624),
  Collimonas = rgb(0.399411, 0.000859, 0.656133),
  Mycoavidus = rgb(0.494877, 0.01199, 0.657865),
  Caballerionia = rgb(0.584391, 0.068579, 0.632812),
  Mitsuaria = rgb(0.665129, 0.138566, 0.585582),
  Caballeronia = rgb(0.736019, 0.209439, 0.527908),
  Bordetella = rgb(0.798216, 0.280197, 0.469538),
  Cupriavidus = rgb(0.853319, 0.351553, 0.413734),
  Achromobacter = rgb(0.901807, 0.425087, 0.359688),
  Herbaspirillum = rgb(0.942598, 0.502639, 0.305816),
  Variovorax = rgb(0.973416, 0.585761, 0.25154),
  Paraburkholderia = rgb(0.991365, 0.675355, 0.198453),
  Acidovorax = rgb(0.993033, 0.77172, 0.154808),
  Pandoraea = rgb(0.974443, 0.874622, 0.144061)
)

p <- ggtree(tree, layout = "circular") %<+% genus.df +
  geom_tippoint(aes(color = genus), size = 2) +  
  geom_tiplab(aes(label = paste0("italic('", label, "')")), 
              parse = TRUE, 
              size = 3,
              offset = 0.015) +
  scale_color_manual(values = burk.colormap) +
  theme_tree() + 
  theme(legend.position = "none")

ggsave("../figures/parts/total_burkholderiales_tree.pdf", plot = p, device = "pdf", width = 8, height = 8.1)
ggsave("../figures/parts/total_burkholderiales_tree.png", plot = p, device = "png", width = 8, height = 8.1, dpi=600)
