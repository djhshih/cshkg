library(dplyr)
library(factoextra)
library(NbClust)
library(cluster)
library(dendextend)
library(pheatmap)
# Touchstone Matrix with gene names in rows (filtered according to inst info)
load("m_genes.rds")
# Convert log scaled matrix to z-scores
m_genes <- scale(log(m_genes + 1))

# Hierarchical Clustering
# 1. Choose Distance Metric
## Calculate Euclidean and Manhattan distances
d_euc <- dist(m_genes, method = "euclidean")
d_man <- dist(m_genes, method="manhattan")
## Hierarchically cluster dataset with both distance metric using ward method
hc_euc <- hclust(d_euc, method = "ward.D2")
hc_man <- hclust(d_man, method = "ward.D2")
## Check the correlation coefficient (distance metric vs cophenetic distance)
coph_euc <- cophenetic(hc_euc)
cor(d_euc, coph_euc) # 0.5411949
coph_man <- cophenetic(hc_man)
cor(d_man, coph_man) # 0.5461158
## Continue the analysis with Manhattan distance metric

# Determining The Optimal Number Of Clusters (k)
# Elbow method
fviz_nbclust(m_genes, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")
# Silhouette method
fviz_nbclust(m_genes, kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")
# Gap statistic
# set.seed(123)
# fviz_nbclust(m_genes, kmeans, nstart = 25,  method = "gap_stat", nboot = 50) +
#   labs(subtitle = "Gap statistic method") # vector memory exhausted
# Take a sample of samples
load("samples.rds")
random_samples <- samples %>% sample_n(1)
m_sub <- m_genes[, intersect(colnames(m_genes), random_samples$inst_id)]

set.seed(13)
clusGap(m_sub, kmeans, 10, B = 100, verbose = interactive())
## Continue the analysis with k = 2.

# Method 1. Ward's Minimum Variance Method
hc_w <- hclust(d_man, method = "ward.D2")
fviz_dend(hc_w,cex=.5)
grupward <- cutree(hc_w, k = 2)
table(grupward)
# Visualize dendogram in colors
fviz_dend(hc_w, k = 2, cex = 0.5, 
          color_labels_by_k = TRUE, rect = TRUE)
result_w <- cbind(rownames(m_genes), Cluster = grupward)
result_w <- as.data.frame(result_w[,2])
# Visualization by Heatmap
pheatmap(m_sub, annotation_row = result_w, cluster_cols = FALSE)
row_dend <- hc_w
heatmap(m_sub, cluster_cols = FALSE, 
        cluster_rows = color_branches(row_dend, k = 2))

# Method 2. Average Linkage Method
hc_avg <- hclust(d_man, method = "average")
fviz_dend(hc_avg,cex=.5)
grupavg <- cutree(hc_avg, k = 4)
table(grupavg)
head(grupavg)
# Visualize dendogram in colors
fviz_dend(hc_avg, k = 4, cex = 0.5, 
          color_labels_by_k = TRUE, rect = TRUE)
result_avg <- cbind(rownames(m_genes), Cluster = grupavg)
result_avg <- as.data.frame(result_avg[,2])
# Visualization by Heatmap
pheatmap(m_sub, annotation_row = result_avg, cluster_cols = FALSE)
row_dend2 <- hc_avg
heatmap(m_sub, cluster_cols = FALSE, 
        cluster_rows = color_branches(row_dend2, k = 4))


