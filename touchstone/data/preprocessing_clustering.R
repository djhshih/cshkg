library(cmapR)
library(tidyr)
library(dplyr)
my_ds <- parse_gctx("GSE92742_Broad_LINCS_Level2_GEX_delta_n49216x978.gctx")

# Group Clustering
library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
# Touchstone Dataset Metadata
m <- my_ds@mat

# group clustering
library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
# Perform hierarchical clustering on the transposed dataset (samples in rows, genes in columns)
m_rev <- na.omit(t(m))
#d <- dist(m_rev)
#hc <- hclust(d) # not working

# Subset matrix
subset_size <- 100
num_columns <- ncol(m_rev)
# Ensure that the subset size is not larger than the number of columns
if (subset_size > num_columns) {
  stop("Subset size is larger than the number of columns in the dataset.")
}

# Error: Error: vector memory exhausted (limit reached?)
# Sys.setenv("R_MAX_VSIZE" = 32000000000) # not working

subset_m <- m_rev[, sample(num_columns, size = subset_size, replace = FALSE)]
saveRDS(subset_m, file = "subset_m.rds")
hc <- hclust(dist(as.matrix(subset_m)))

# Plot the dendrogram
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = 4)
# Visualize the dendrogram
plot(dend, main = "Hierarchical Clustering Dendrogram", xlab = "Samples", sub = "")

# Cut the tree to get clusters
clusters <- cutree(hc, h = 50)  # Adjust the height (h) based on the dendrogram
dim(df_group)
length(clusters)
# Combine the cluster information with your sample grouping data
result_df <- cbind(df_group, Cluster = clusters)