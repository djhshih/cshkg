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
clusters <- cutree(hc, h = 200000)  # Adjust the height (h) based on the dendrogram
dim(genes_info)
length(clusters)
# Combine the cluster information with your sample grouping data
result <- cbind(genes_info, Cluster = clusters)
unique(result$Cluster) # total 279 groups of genes

# Count frequency of genes in each group
result <- result %>%
  group_by(Cluster) %>%
  mutate(count_freq = n())
# Filter groups less than 10 genes
result <- subset(result, count_freq >= 10) # total 516 genes (from 978 genes)

# Convert to gene symbol
result$id <- as.numeric(result$id)
gene_info_delta_landmark$pr_gene_id <- as.numeric(gene_info_delta_landmark$pr_gene_id)
gene_filtered <- gene_info_delta_landmark[na.omit(match(result$id, gene_info_delta_landmark$pr_gene_id)),]

all(result$id == gene_filtered$pr_gene_id) # TRUE
# When FALSE, sort genes
# sort_gene <- result$id
# gene_filtered <- gene_filtered[order(match(gene_filtered$pr_gene_id, sort_gene)),]

rownames(result) <- gene_filtered$pr_gene_symbol
all(rownames(result) == gene_filtered$pr_gene_symbol) #TRUE
unique(result$Cluster) # total 20 groups of genes

save(result, file = "result.rds")

# New matrix with clustered groups of genes (98 groups of genes x 49207 samples)
m_cluster <- m_genes[intersect(rownames(m_genes), rownames(result)),]
