# Sample Group Clustering
library(dplyr)
library(factoextra)
library(NbClust)
library(cluster)
library(dendextend)
library(pheatmap)
library(tsne)
library(Rtsne)
library(plotly)
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

# Hierarchical Clustering: `VCAP`
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

# Determining The Optimal Number Of Clusters (k)
# Elbow method
fviz_nbclust(m_VCAP, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2) +
  labs(subtitle = "Elbow method")
# Silhouette method
fviz_nbclust(m_VCAP, kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")
# Gap statistic
# ERROR: vector memory exhausted
# Take a sample of samples
random_samples <- s_VCAP %>% sample_n(1)
sub_mVCAP <- m_VCAP[intersect(rownames(m_VCAP), random_samples$inst_id), ]

clusGap(sub_mVCAP, FUN = kmeans, K.max = 10, B = 500, verbose = interactive())
fviz_nbclust(sub_mVCAP, kmeans, k.max = 10,  method = "gap_stat", nboot = 500) +
  labs(subtitle = "Gap statistic method")
## Continue the analysis with k = 2 based on PCA result.

# Method 1. Ward's Minimum Variance Method
hc_VCAP <- hclust(d_euc_VCAP, method = "ward.D2")
grupward <- cutree(hc_VCAP, k = 2) # based on PCA result k = 2
table(grupward)
# Visualize dendogram in colors
fviz_dend(hc_VCAP, k = 2, cex = 0.5, 
          color_labels_by_k = TRUE, rect = TRUE)
result_VCAP <- cbind(rownames(m_VCAP), Cluster = grupward)
result_VCAP <- as.data.frame(result_VCAP[,2])
# Visualization by Heatmap
pheatmap(m_VCAP, annotation_row = result_VCAP, cluster_cols = FALSE)
row_dend <- hc_VCAP
heatmap(m_VCAP, cluster_cols = FALSE, 
        cluster_rows = color_branches(row_dend, k = 2))

# Method 2. Average Linkage Method
hc_VCAP_avg <- hclust(d_euc_VCAP, method = "average")
grupavg <- cutree(hc_VCAP_avg, k = 2)
table(grupavg)
head(grupavg)
# Visualize dendogram in colors
fviz_dend(hc_VCAP_avg, k = 2, cex = 0.5, 
          color_labels_by_k = TRUE, rect = TRUE)
result_VCAP_avg <- cbind(rownames(m_VCAP), Cluster = grupavg)
result_VCAP_avg <- as.data.frame(result_VCAP_avg[,2])
# Visualization by Heatmap
pheatmap(m_VCAP, annotation_row = result_VCAP_avg, cluster_cols = FALSE)
row_dend2 <- hc_VCAP_avg
heatmap(m_VCAP, cluster_cols = FALSE, 
        cluster_rows = color_branches(row_dend2, k = 2))
## Continue with ward method where k = 2.

# Create a clustered expression matrix
all(s_VCAP$inst_id == rownames(result_VCAP))
result_VCAP$group <- s_VCAP$group
colnames(result_VCAP)[1] <- "cluster_group"
unique(result_VCAP$group) # Original VCAP samples had total 159 groups
VCAP_clus1 <- filter(result_VCAP, cluster_group %in% "1")
unique(VCAP_clus1$group) # total 158 groups clustered
VCAP_clus2 <- filter(result_VCAP, cluster_group %in% "2")
unique(VCAP_clus2$group) # total 82 groups clustered
# Find sample groups only in cluster 1 group but not in cluster 2 group
setdiff(VCAP_clus1[, 2], VCAP_clus2[, 2])
# Find sample groups only in cluster 2 group but not in cluster 1 group
setdiff(VCAP_clus2[, 2], VCAP_clus1[, 2])

VCAP_clusgrp <- c("VCAP_ACVR2B", "VCAP_BMPR2", "VCAP_PIK3CA", "VCAP_AXL", "VCAP_RPS6KA1", "VCAP_BMPR1B",
                  "VCAP_CHEK2", "VCAP_ARAF", "VCAP_ABL2", "VCAP_ADRBK2", "VCAP_MAP2K4", "VCAP_AK4",
                  "VCAP_ETNK2", "VCAP_AK1", "VCAP_CAMK2D", "VCAP_CDK5", "VCAP_CKS2", "VCAP_BTK",
                  "VCAP_CSNK1A1", "VCAP_CALM3", "VCAP_CDK8", "VCAP_CAMK4", "VCAP_MAP3K8", "VCAP_CDK9",
                  "VCAP_CDK11B", "VCAP_CSNK1D", "VCAP_MAPK15", "VCAP_PIK3CB", "VCAP_BMPR1A", "VCAP_BUB1B",
                  "VCAP_CSNK1G2", "VCAP_BRDT", "VCAP_CDK3", "VCAP_RFP", "VCAP_EMPTY_VECTOR", "VCAP_PANK4",
                  "VCAP_PAK7", "VCAP_MORN1", "VCAP_STK33", "VCAP_UCKL1", "VCAP_KIAA1804", "VCAP_RPS6KL1", 
                  "VCAP_MINK1", "VCAP_TTBK1", "VCAP_TEX14", "VCAP_ADPGK", "VCAP_CDKL4", "VCAP_FUK", 
                  "VCAP_LOC392265", "VCAP_LOC390877", "VCAP_LOC392226", "VCAP_LOC400301", "VCAP_MYO3B", "VCAP_MARK2P10",
                  "VCAP_DCLK3", "VCAP_ITGB1BP3", "VCAP_NLK", "VCAP_WNK1", "VCAP_AGK", "VCAP_LRRK1",
                  "VCAP_EEF2K", "VCAP_SH3BP4", "VCAP_GK5", "VCAP_NEK8", "VCAP_MASTL", "VCAP_DGKK",
                  "VCAP_NME3", "VCAP_LOC441971", "VCAP_LOC392347", "VCAP_PLXNA3", "VCAP_FGGY", "VCAP_MEX3B",
                  "VCAP_XRCC6BP1", "VCAP_FLT1", "VCAP_PNCK", "VCAP_WNK4", "VCAP_IGFN1")
s_VCAP$new_group <- s_VCAP$group
s_VCAP$new_group[s_VCAP$group %in% VCAP_clusgrp] <- "VCAP_clusgrp"
save(s_VCAP, file = "s_VCAP.rds")

# Perform t-SNE for `VCAP` cell type
library(Rtsne)
## 1. Compute t-SNE embedding for `m_VCAP`
tsne_VCAP <- Rtsne(m_VCAP, dims = 3, perplexity = 50, verbose = TRUE, max_iter = 500)
## 2. Visualize t-SNE projections
tsne_df_VCAP <- as.data.frame(tsne_VCAP$Y)
plot_ly(data = tsne_df_VCAP, x = ~V1, y = ~V2, z = ~V3, color = s_VCAP$group) |> add_markers(size = 3)
## 3. Hierarchical clustering using t-SNE embeddings
tsne_dist_VCAP <- dist(tsne_VCAP$Y)
hc_tsne_VCAP <- hclust(tsne_dist_VCAP, method = "ward.D2")
fviz_dend(hc_tsne_VCAP, k = 8, cex = 0.5, 
          color_labels_by_k = TRUE, rect = TRUE)
clusters_VCAP <- cutree(hc_tsne_VCAP, k = 8)
VCAP_tsne_grp <- cbind(rownames(m_VCAP), Cluster = clusters_VCAP)
rownames(VCAP_tsne_grp) <- VCAP_tsne_grp[,1]
VCAP_tsne_grp <- as.data.frame(VCAP_tsne_grp[,2])
# Create a clustered expression matrix
all(s_VCAP$inst_id == rownames(VCAP_tsne_grp)) # TRUE
VCAP_tsne_grp$group <- s_VCAP$group
colnames(VCAP_tsne_grp)[1] <- "cluster_group"
# m_clus_VCAP <- VCAP_tsne_grp |> group_by(cluster_group)

# Perform t-SNE for `MCF7` cell type
tsne_MCF7 <- Rtsne(m_MCF7, dims = 3, perplexity = 50, verbose = TRUE, max_iter = 1000)
tsne_df_MCF7 <- as.data.frame(tsne_MCF7$Y)
plot_ly(data = tsne_df_MCF7, x = ~V1, y = ~V2, z = ~V3, color = s_MCF7$group) |> add_markers(size = 3)

# Perform t-SNE for `PC3` cell type
tsne_PC3 <- Rtsne(m_PC3, dims = 3, perplexity = 50, verbose = TRUE, max_iter = 1000)
tsne_df_PC3 <- as.data.frame(tsne_PC3$Y)
plot_ly(data = tsne_df_PC3, x = ~V1, y = ~V2, z = ~V3, color = s_PC3$group) |> add_markers(size = 3)
