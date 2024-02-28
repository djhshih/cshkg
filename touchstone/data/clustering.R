# Sample Group Clustering
library(dplyr)
library(factoextra)
library(NbClust)
library(cluster)
library(dendextend)
library(pheatmap)
# Load log-scaled Touchstone Matrix with gene names in rows (filtered according to inst info)
load("m_genes.rds")
load("m_VCAP.rds")
load("s_VCAP.rds")
load("m_MCF7.rds")
load("s_MCF7.rds")
load("m_PC3.rds")
load("s_PC3.rds")

# Prepare matrices for sample clustering
m_genes <- t(m_genes)
m_VCAP <- t(m_VCAP)
m_MCF7 <- t(m_MCF7)
m_PC3 <- t(m_PC3)

# Hierarchical Clustering
# 1. Choose Distance Metric
## Calculate Euclidean and Manhattan distances
d_euc_VCAP <- dist(m_VCAP, method = "euclidean")
d_man_VCAP <- dist(m_VCAP, method="manhattan")
## Hierarchically cluster dataset with both distance metric using ward method
hc_euc_VCAP <- hclust(d_euc_VCAP, method = "ward.D2")
hc_man_VCAP <- hclust(d_man_VCAP, method = "ward.D2")
## Check the correlation coefficient (distance metric vs cophenetic distance)
coph_euc_VCAP <- cophenetic(hc_euc_VCAP)
cor(d_euc_VCAP, coph_euc_VCAP) # 0.4887777
coph_man_VCAP <- cophenetic(hc_man_VCAP)
cor(d_man_VCAP, coph_man_VCAP) # 0.4570925
## Continue the analysis with Euclidean distance metric for VCAP cell type
# Check for MCF7 cell type
d_euc_MCF7 <- dist(m_MCF7, method = "euclidean")
d_man_MCF7 <- dist(m_MCF7, method="manhattan")
hc_euc_MCF7 <- hclust(d_euc_MCF7, method = "ward.D2")
hc_man_MCF7 <- hclust(d_man_MCF7, method = "ward.D2")
coph_euc_MCF7 <- cophenetic(hc_euc_MCF7)
coph_man_MCF7 <- cophenetic(hc_man_MCF7)
cor(d_euc_MCF7, coph_euc_MCF7)
cor(d_man_MCF7, coph_man_MCF7)
## Continue the analysis with ?? distance metric for MCF7 cell type
# Check for PC3 cell type
d_euc_PC3 <- dist(m_PC3, method = "euclidean")
d_man_PC3 <- dist(m_PC3, method="manhattan")
hc_euc_PC3 <- hclust(d_euc_PC3, method = "ward.D2")
hc_man_PC3 <- hclust(d_man_PC3, method = "ward.D2")
coph_euc_PC3 <- cophenetic(hc_euc_PC3)
coph_man_PC3 <- cophenetic(hc_man_PC3)
cor(d_euc_PC3, coph_euc_PC3)
cor(d_man_PC3, coph_man_PC3)
## Continue the analysis with ?? distance metric for PC3 cell type

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

## Number of samples in group clusters are reasonably high in average method.
## Continue with ward method where k = 4.

grupward <- cutree(hc_w, k = 4)
table(grupward)
# Visualize dendogram in colors
fviz_dend(hc_w, k = 4, cex = 0.5, 
          color_labels_by_k = TRUE, rect = TRUE)
result_w <- cbind(rownames(m_genes), Cluster = grupward)
result_w <- as.data.frame(result_w[,2])
# Visualization by Heatmap
pheatmap(m_sub, annotation_row = result_w, cluster_cols = FALSE)

# Create a clustered gene expression matrix
result_w$gene <- rownames(result_w)
result_w <- result_w[order(result_w$`result_w[, 2]`),]
result_w <- result_w[1]
m_ordered <- m_genes[rownames(result_w),]

save(m_ordered, file = "m_ordered.rds")

