library(cmapR)
library(tidyr)
library(dplyr)
my_ds <- parse_gctx("GSE92742_Broad_LINCS_Level2_GEX_delta_n49216x978.gctx")
# sample info
inst_info <- read.delim("GSE92742_Broad_LINCS_inst_info.txt")
# gene info
gene_info_delta_landmark <- read.delim("GSE92742_Broad_LINCS_gene_info_delta_landmark.txt")

# Group Clustering
library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
# Touchstone Dataset Metadata
m <- my_ds@mat
genes_info <- my_ds@rdesc

# Group Clustering
# Hierarchical clustering (samples in columns, genes in rows)
d <- dist(m)
hc <- hclust(d)
# Plot the dendrogram
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = 10)
# Visualize the dendrogram
plot(dend, main = "Hierarchical Clustering Dendrogram", xlab = "Samples", sub = "")

# Cut the tree to get clusters
clusters <- cutree(hc, h = 400000)  # Adjust the height (h) based on the dendrogram
dim(genes_info)
length(clusters)
# Combine the cluster information with your sample grouping data
result_df <- cbind(genes_info, Cluster = clusters)
unique(result_df$Cluster) # total 44 groups
