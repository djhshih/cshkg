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
# hc <- hclust(d, method = "complete")

# Plot the dendrogram
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = 10)
# Visualize the dendrogram
plot(dend, main = "Hierarchical Clustering Dendrogram", xlab = "Samples", sub = "")

# Try ward and average methods for clustering
heatmap(m, scale = "none")
hc1 <- agnes(m, method = "ward")
pltree(hc1, cex = 0.6, hang = -1, main = "Dendrogram of agnes: ward's method")
hc2 <- agnes(m, method = "average")
pltree(hc2, cex = 0.6, hang = -1, main = "Dendrogram of agnes: average linkage")

# Function to compute coefficient
methods_ <- c("average", "ward")
names(m) <- c("average", "ward")
ac <- function(x) {
  agnes(m, method = x)$ac
}
map_dbl(methods_, ac)

# Cut the tree to get clusters
clusters <- cutree(hc, h = 200000)  # Adjust the height (h) based on the dendrogram
clusters1 <- cutree(hc1, h = 3.0e+05)
clusters2 <- cutree(hc2, h = 200000)
dim(genes_info)
length(clusters)
# Combine the cluster information with your sample grouping data
result <- cbind(genes_info, Cluster = clusters, Cluster_ward = clusters1, Cluster_average = clusters2)
length(unique(result$Cluster)) # total 279 groups of genes
length(unique(result$Cluster_ward)) # total 306 groups of genes
length(unique(result$Cluster_average)) # total 223 groups of genes

# Count frequency of genes in each group
result <- result %>%
  group_by(Cluster) %>%
  mutate(count_freq = n())

result <- result %>%
  group_by(Cluster_ward) %>%
  mutate(count_freq1 = n())

result <- result %>%
  group_by(Cluster_average) %>%
  mutate(count_freq2 = n())

# Filter groups less than 10 genes
sub_result <- subset(result, count_freq >= 10) # total 516 genes (from 978 genes)
length(unique(sub_result$Cluster)) # total 20 groups of genes
sub_result1 <- subset(result, count_freq1 >= 10) # total 705 genes (from 978 genes)
length(unique(sub_result1$Cluster_ward)) # total 32 groups of genes
sub_result2 <- subset(result, count_freq2 >= 10) # total 674 genes (from 978 genes)
length(unique(sub_result2$Cluster_average)) # total 9 groups of genes
# Continue with the result of ward method

# Convert to gene symbol
sub_result1$id <- as.numeric(sub_result1$id)
gene_info_delta_landmark$pr_gene_id <- as.numeric(gene_info_delta_landmark$pr_gene_id)
gene_filtered <- gene_info_delta_landmark[na.omit(match(sub_result1$id, gene_info_delta_landmark$pr_gene_id)),]

all(sub_result1$id == gene_filtered$pr_gene_id) # TRUE
# When FALSE, sort genes
# sort_gene <- sub_result1$id
# gene_filtered <- gene_filtered[order(match(gene_filtered$pr_gene_id, sort_gene)),]

rownames(sub_result1) <- gene_filtered$pr_gene_symbol
all(rownames(sub_result1) == gene_filtered$pr_gene_symbol) #TRUE

# Convert result to data frame
cluster_df <- sub_result1[3]
rownames(cluster_df) <- rownames(sub_result1)

save(result, file = "result.rds")
save(sub_result1, file = "sub_result1.rds")

# Draw heatmap with gene cluster info
# pheatmap(m_genes, annotation_row = cluster_df) #vector memory exhausted
# Take a sample of samples
load("samples.rds")
library(dplyr)
library(pheatmap)
random_samples <- samples %>% sample_n(1)
m_cluster <- m_genes[intersect(rownames(m_genes), rownames(cluster_df)),
                     intersect(colnames(m_genes), random_samples$inst_id)]
# deal Error in cut.default(a, breaks = 100) : 'x' must be numeric
# m_num <- as.data.frame(lapply(m_cluster,
#                               function(x) as.numeric(as.character(x))))

pheatmap(m_cluster, annotation_row = cluster_df)
