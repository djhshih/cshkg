# Load required libraries for group clustering.
library(io)
library(dplyr)
library(factoextra)
library(cluster)
library(dendextend)
library(tsne)
library(Rtsne)
library(ggplot2)
library(ggrepel)
library(ggpubr)
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

# Prepare matrices for group clustering
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
## Continue the analysis with Euclidean distance metric for VCAP cell line

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

# Method: Ward's Minimum Variance Method
hc_VCAP <- hclust(d_euc_VCAP, method = "ward.D2")
grupward <- cutree(hc_VCAP, k = 2) # based on PCA result k = 2
table(grupward)
# Visualize dendogram in colors
fviz_dend(hc_VCAP, k = 2, cex = 0.5, 
          color_labels_by_k = TRUE, rect = TRUE)
result_VCAP <- cbind(rownames(m_VCAP), Cluster = grupward)
result_VCAP <- as.data.frame(result_VCAP[,2])
### Hierarchical clustering function cannot be runned for `MCF7` and `PC3` matrixes (memory issue).
### Group clustering with alternative method: t-SNE and Louvain clustering that does not require distance matrix.

# Continue with t-SNE analysis, using initial reduction data through PCA
load("VCAP_pca_prcomp_data.rds")
load("MCF7_pca_prcomp_data.rds")
load("PC3_pca_prcomp_data.rds")

# 1.`VCAP` cell line
## a. Compute t-SNE embedding for `m_VCAP`
tsne_VCAP <- Rtsne(VCAP_pca_prcomp_data, pca = FALSE)
## b. Visualize t-SNE projections
tsne_data_VCAP <- s_VCAP[, 1:2]
tsne_data_VCAP <- tsne_data_VCAP |> ungroup() |> 
  mutate(tsne1 = tsne_VCAP$Y[, 1], tsne2 = tsne_VCAP$Y[, 2])
p_tsne_VCAP <- ggplot(tsne_data_VCAP, aes(x = tsne1, y = tsne2, colour = group)) + 
  geom_point(alpha = 0.3) + theme(legend.position = "none") +
  labs(x = "tSNE dimension 1", y = "tSNE dimension 2",
        title = "tSNE result of VCAP cell line colored by perturbation type")
p_tsne_VCAP
### k >= 12 for `VCAP` cell line
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
p_louvain_VCAP <- ggplot(tsne_data_VCAP, aes(x = tsne1, y = tsne2, colour = louvain_group)) + 
  geom_point(alpha = 0.3) + theme_bw() + 
  geom_label_repel(aes(label = louvain_group), data = lc_data_VCAP) + 
  guides(colour = FALSE) + labs(x = "tSNE dimension 1", y = "tSNE dimension 2",
                                title = "Result of Louvain Clustering of VCAP cell line")
p_louvain_VCAP

# 2.`MCF7` cell line
## a. t-SNE
tsne_MCF7 <- Rtsne(MCF7_pca_prcomp_data, pca = FALSE)
tsne_data_MCF7 <- s_MCF7[, 1:2]
tsne_data_MCF7 <- tsne_data_MCF7 |> ungroup() |> 
  mutate(tsne1 = tsne_MCF7$Y[, 1], tsne2 = tsne_MCF7$Y[, 2])
p_tsne_MCF7 <- ggplot(tsne_data_MCF7, aes(x = tsne1, y = tsne2, colour = group)) + 
  geom_point(alpha = 0.3) + theme(legend.position = "none") +
  labs(x = "tSNE dimension 1", y = "tSNE dimension 2",
       title = "tSNE result of MCF7 cell line colored by perturbation type")
p_tsne_MCF7
### k >= 30 for `MCF7` cell line
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
p_louvain_MCF7 <- ggplot(tsne_data_MCF7, aes(x = tsne1, y = tsne2, colour = louvain_group)) + 
  geom_point(alpha = 0.3) + theme_bw() + 
  geom_label_repel(aes(label = louvain_group), data = lc_data_MCF7) +
  guides(colour = FALSE) + labs(x = "tSNE dimension 1", y = "tSNE dimension 2",
                                title = "Result of Louvain Clustering of MCF7 cell line")
p_louvain_MCF7

# 3.`PC3` cell line
## a. t-SNE
tsne_PC3 <- Rtsne(PC3_pca_prcomp_data, pca = FALSE)
tsne_data_PC3 <- s_PC3[, 1:2]
tsne_data_PC3 <- tsne_data_PC3 |> ungroup() |> 
  mutate(tsne1 = tsne_PC3$Y[, 1], tsne2 = tsne_PC3$Y[, 2])
p_tsne_PC3 <- ggplot(tsne_data_PC3, aes(x = tsne1, y = tsne2, colour = group)) + 
  geom_point(alpha = 0.3) + theme(legend.position = "none") +
  labs(x = "tSNE dimension 1", y = "tSNE dimension 2",
       title = "tSNE result of PC3 cell line colored by perturbation type")
p_tsne_PC3
### k >= 41 for `PC3` cell line
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
p_louvain_PC3 <- ggplot(tsne_data_PC3, aes(x = tsne1, y = tsne2, colour = louvain_group)) + 
  geom_point(alpha = 0.3) + theme_bw() + 
  geom_label_repel(aes(label = louvain_group), data = lc_data_PC3) +
  guides(colour = FALSE) + labs(x = "tSNE dimension 1", y = "tSNE dimension 2",
                                title = "Result of Louvain Clustering of PC3 cell line")
p_louvain_PC3

# Combine plots of t-SNE for each matrixes
qdraw(
  ggarrange(p_tsne_VCAP, p_tsne_MCF7, p_tsne_PC3, 
          labels = c("A", "B", "C"), ncol = 3, nrow = 1),
  width = 18, height = 5,
  file = "../plots/tSNE_result.png"
)

# Combine plots of Louvain clustering for each matrixes
qdraw(
  ggarrange(p_louvain_VCAP, p_louvain_MCF7, p_louvain_PC3, 
          labels = c("A", "B", "C"), ncol = 3, nrow = 1),
  width = 18, height = 5,
  file = "../plots/Louvain_result.png"
)

tsne_data_VCAP$louvain_group <- paste(tsne_data_VCAP$louvain_group, "VCAP", sep = "_")
tsne_data_MCF7$louvain_group <- paste(tsne_data_MCF7$louvain_group, "MCF7", sep = "_")
tsne_data_PC3$louvain_group <- paste(tsne_data_PC3$louvain_group, "PC3", sep = "_")
samples_tsne <- rbind(tsne_data_VCAP, tsne_data_MCF7, tsne_data_PC3)
samples_tsne <- samples_tsne |> select(c("inst_id", "group", "louvain_group"))
samples_tsne <- samples_tsne[order(match(samples_tsne$inst_id, rownames(m_genes))),]
all(rownames(m_genes) == samples_tsne$inst_id) # TRUE

# Chose representative perturbagen for each clustered groups
grouped_df <- samples_tsne |> group_by(louvain_group, group) |> count() |> ungroup()
representative <- grouped_df |> group_by(louvain_group) |> 
  ## sort by frequency in descending and perturbagen group in ascending
  arrange(desc(n), group) |> mutate(rank = row_number()) |> ungroup()
select_representative <- representative[representative$rank == 1, ]
select_representative <- select_representative[, c(1,2)]
select_representative$sep <- select_representative$louvain_group
select_representative <- select_representative |>
  separate(sep, into = c("number", "cell_line"), sep = "_", extra = "merge", fill = "right")
select_representative$sep2 <- select_representative$group
select_representative <- select_representative |>
  separate(sep2, into = c("cell_line", "perturbagen"), sep = "_", extra = "merge", fill = "right")
select_representative <- select_representative |> unite(louvain_group_name, group, number, sep = "_")
# Update on new sample dataset
samples_cluster <- merge(samples_tsne, select_representative, by = "louvain_group", all.x = TRUE)
samples_cluster <- samples_cluster[, c(2,3,1,4,5,6)]

# Save the outcomes of clustered groups
save(samples_cluster, file = "samples_cluster.rds")
