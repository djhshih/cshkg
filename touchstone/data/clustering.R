# Sample Group Clustering
library(dplyr)
library(factoextra)
library(cluster)
library(dendextend)
library(pheatmap)
library(tsne)
library(Rtsne)
library(ggplot2)
library(ggrepel)
library(plotly)
library(igraph)
library(FNN)
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

# Method: Ward's Minimum Variance Method
hc_VCAP <- hclust(d_euc_VCAP, method = "ward.D2")
grupward <- cutree(hc_VCAP, k = 2) # based on PCA result k = 2
table(grupward)
# Visualize dendogram in colors
fviz_dend(hc_VCAP, k = 2, cex = 0.5, 
          color_labels_by_k = TRUE, rect = TRUE)
result_VCAP <- cbind(rownames(m_VCAP), Cluster = grupward)
result_VCAP <- as.data.frame(result_VCAP[,2])

# Create new sample groups based on hierarchical clustering result for `VCAP` cell type
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
s_VCAP$hierarchical_group <- s_VCAP$group
s_VCAP$hierarchical_group[s_VCAP$group %in% VCAP_clusgrp] <- "VCAP_clusgrp"
### Hierarchical clustering function cannot be runned for `MCF7` and `PC3` matrixes.
### Use t-SNE and Louvain Clustering method that does not require distance matrix.

# Continue with t-SNE analysis, using initial reduction through PCA
load("VCAP_pca_prcomp_data.rds")
load("MCF7_pca_prcomp_data.rds")
load("PC3_pca_prcomp_data.rds")

# 1.`VCAP` cell type
## a. Compute t-SNE embedding for `m_VCAP`
tsne_VCAP <- Rtsne(VCAP_pca_prcomp_data, pca = FALSE)
## b. Visualize t-SNE projections
tsne_data_VCAP <- s_VCAP[, 1:2]
tsne_data_VCAP <- tsne_data_VCAP |> ungroup() |> 
  mutate(tsne1 = tsne_VCAP$Y[, 1], tsne2 = tsne_VCAP$Y[, 2])
ggplot(tsne_data_VCAP, aes(x = tsne1, y = tsne2, colour = group)) + 
  geom_point(alpha = 0.3) + theme(legend.position = "none") +
  labs(x = "tSNE dimension 1", y = "tSNE dimension 2",
        title = "tSNE result of VCAP cell type colored by perturbation type")
plot_ly(data = tsne_data_VCAP, x = ~tsne1, y = ~tsne2, color = ~group) |> add_markers(size = 3)
### k >= 12 for `VCAP` cell type
## c. Louvain clustering
knn_VCAP <- get.knn(as.matrix(tsne_VCAP$Y), k = 16)
knn_VCAP <- data.frame(from = rep(1:nrow(knn_VCAP$nn.index), 16), 
                       to = as.vector(knn_VCAP$nn.index), 
                       weight = 1/(1 + as.vector(knn_VCAP$nn.dist)))
nw_VCAP <- graph_from_data_frame(knn_VCAP, directed = FALSE)
nw_VCAP <- simplify(nw_VCAP)
lc_VCAP <- cluster_louvain(nw_VCAP)
tsne_data_VCAP$louvain_group <- as.factor(membership(lc_VCAP))
lc_data_VCAP <- tsne_data_VCAP |> group_by(louvain_group) |> 
  select(tsne1,tsne2) |> summarize_all(mean)
ggplot(tsne_data_VCAP, aes(x = tsne1, y = tsne2, colour = louvain_group)) + 
  geom_point(alpha = 0.3) + theme_bw() + 
  geom_label_repel(aes(label = louvain_group), data = lc_data_VCAP) + 
  guides(colour = FALSE) + labs(x = "tSNE dimension 1", y = "tSNE dimension 2",
                                title = "Result of Louvain Clustering of VCAP cell type")

# 2.`MCF7` cell type
## a. t-SNE
tsne_MCF7 <- Rtsne(MCF7_pca_prcomp_data, pca = FALSE)
tsne_data_MCF7 <- s_MCF7[, 1:2]
tsne_data_MCF7 <- tsne_data_MCF7 |> ungroup() |> 
  mutate(tsne1 = tsne_MCF7$Y[, 1], tsne2 = tsne_MCF7$Y[, 2])
ggplot(tsne_data_MCF7, aes(x = tsne1, y = tsne2, colour = group)) + 
  geom_point(alpha = 0.3) + theme(legend.position = "none") +
  labs(x = "tSNE dimension 1", y = "tSNE dimension 2",
       title = "tSNE result of MCF7 cell type colored by perturbation type")
plot_ly(data = tsne_data_MCF7, x = ~tsne1, y = ~tsne2, color = ~group) |> add_markers(size = 3)
### k >= 30 for `MCF7` cell type
## b. Louvain clustering
knn_MCF7 <- get.knn(as.matrix(tsne_MCF7$Y), k = 30)
knn_MCF7 <- data.frame(from = rep(1:nrow(knn_MCF7$nn.index), 30), 
                      to = as.vector(knn_MCF7$nn.index), 
                      weight = 1/(1 + as.vector(knn_MCF7$nn.dist)))
nw_MCF7 <- graph_from_data_frame(knn_MCF7, directed = FALSE)
nw_MCF7 <- simplify(nw_MCF7)
lc_MCF7 <- cluster_louvain(nw_MCF7)
tsne_data_MCF7$louvain_group <- as.factor(membership(lc_MCF7))
lc_data_MCF7 <- tsne_data_MCF7 |> group_by(louvain_group) |> 
  select(tsne1,tsne2) |> summarize_all(mean)
ggplot(tsne_data_MCF7, aes(x = tsne1, y = tsne2, colour = louvain_group)) + 
  geom_point(alpha = 0.3) + theme_bw() + 
  geom_label_repel(aes(label = louvain_group), data = lc_data_MCF7) +
  guides(colour = FALSE) + labs(x = "tSNE dimension 1", y = "tSNE dimension 2",
                                title = "Result of Louvain Clustering of MCF7 cell type")

# 3.`PC3` cell type
## a. t-SNE
tsne_PC3 <- Rtsne(PC3_pca_prcomp_data, pca = FALSE)
tsne_data_PC3 <- s_PC3[, 1:2]
tsne_data_PC3 <- tsne_data_PC3 |> ungroup() |> 
  mutate(tsne1 = tsne_PC3$Y[, 1], tsne2 = tsne_PC3$Y[, 2])
ggplot(tsne_data_PC3, aes(x = tsne1, y = tsne2, colour = group)) + 
  geom_point(alpha = 0.3) + theme(legend.position = "none") +
  labs(x = "tSNE dimension 1", y = "tSNE dimension 2",
       title = "tSNE result of PC3 cell type colored by perturbation type")
plot_ly(data = tsne_data_PC3, x = ~tsne1, y = ~tsne2, color = ~group) |> add_markers(size = 3)
### k >= 41 for `PC3` cell type
## b. Louvain clustering
knn_PC3 <- get.knn(as.matrix(tsne_PC3$Y), k = 41)
knn_PC3 <- data.frame(from = rep(1:nrow(knn_PC3$nn.index), 41), 
                       to = as.vector(knn_PC3$nn.index), 
                       weight = 1/(1 + as.vector(knn_PC3$nn.dist)))
nw_PC3 <- graph_from_data_frame(knn_PC3, directed = FALSE)
nw_PC3 <- simplify(nw_PC3)
lc_PC3 <- cluster_louvain(nw_PC3)
tsne_data_PC3$louvain_group <- as.factor(membership(lc_PC3))
lc_data_PC3 <- tsne_data_PC3 |> group_by(louvain_group) |> 
  select(tsne1,tsne2) |> summarize_all(mean)
ggplot(tsne_data_PC3, aes(x = tsne1, y = tsne2, colour = louvain_group)) + 
  geom_point(alpha = 0.3) + theme_bw() + 
  geom_label_repel(aes(label = louvain_group), data = lc_data_PC3) +
  guides(colour = FALSE) + labs(x = "tSNE dimension 1", y = "tSNE dimension 2",
                                title = "Result of Louvain Clustering of PC3 cell type")

# Save the outcomes of clustered groups
save(tsne_data_VCAP, file = "tsne_data_VCAP.rds")
save(tsne_data_MCF7, file = "tsne_data_MCF7.rds")
save(tsne_data_PC3, file = "tsne_data_PC3.rds")

# Create new sample groups
# 1. `VCAP`
tsne_data_VCAP |> select(inst_id, group, louvain_group) |> 
  group_by(louvain_group) |> summarize(number = n())

louv_groups_VCAP <- unique(tsne_data_VCAP$louvain_group)
selected_groups_VCAP <- character()
for (i in 1:length(louv_groups_VCAP)){
  louvain_groups <- louv_groups_VCAP[i]
  filtered_data <- tsne_data_VCAP[tsne_data_VCAP$louvain_group == louvain_groups, ]
  
  unique_groups <- unique(filtered_data$group)
  
  is_single_louv_group <- sapply(unique_groups, function(group){
    sum(tsne_data_VCAP$group == group & tsne_data_VCAP$louvain_group != louvain_groups) == 0
  })
  selected_groups_VCAP <- c(selected_groups_VCAP, unique_groups[is_single_louv_group])
}
## No perturbation group in single louvarian group.

# 2. `MCF7`
tsne_data_MCF7 |> select(inst_id, group, louvain_group) |> 
  group_by(louvain_group) |> summarize(number = n())

louv_groups_MCF7 <- unique(tsne_data_MCF7$louvain_group)
selected_groups_MCF7 <- character()
for (i in 1:length(louv_groups_MCF7)){
  louvain_groups <- louv_groups_MCF7[i]
  filtered_data <- tsne_data_MCF7[tsne_data_MCF7$louvain_group == louvain_groups, ]
  
  unique_groups <- unique(filtered_data$group)
  
  is_single_louv_group <- sapply(unique_groups, function(group){
    sum(tsne_data_MCF7$group == group & tsne_data_MCF7$louvain_group != louvain_groups) == 0
  })
  selected_groups_MCF7 <- c(selected_groups_MCF7, unique_groups[is_single_louv_group])
}
## No perturbation group in single louvarian group.

# 3. `PC3`
tsne_data_PC3 |> select(inst_id, group, louvain_group) |> 
  group_by(louvain_group) |> summarize(number = n())

louv_groups_PC3 <- unique(tsne_data_PC3$louvain_group)
selected_groups_PC3 <- character()
for (i in 1:length(louv_groups_PC3)){
  louvain_groups <- louv_groups_PC3[i]
  filtered_data <- tsne_data_PC3[tsne_data_PC3$louvain_group == louvain_groups, ]
  
  unique_groups <- unique(filtered_data$group)
  
  is_single_louv_group <- sapply(unique_groups, function(group){
    sum(tsne_data_PC3$group == group & tsne_data_PC3$louvain_group != louvain_groups) == 0
  })
  selected_groups_PC3 <- c(selected_groups_PC3, unique_groups[is_single_louv_group])
}
## No perturbation group in single louvarian group.
