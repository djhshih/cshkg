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
## Check the correlation coefficient between distance metric and cophenetic distance
coph_euc <- cophenetic(hc_euc)
cor(d_euc, coph_euc) # 0.5411949
coph_man <- cophenetic(hc_man)
cor(d_man, coph_man) # 0.5461158
## Continue the analysis with Manhattan distance metric

# Method 1. Ward's Minimum Variance Method
library(factoextra)
hc_w <- hclust(d_euc, method = "ward.D2")
fviz_dend(hc_w,cex=.5)
grupward <- cutree(hc_w, k = 2)
table(grupward)
head(grupward)
# Visualize dendogram in colors
fviz_dend(hc_w, k = 2, cex = 0.5, 
          color_labels_by_k = TRUE, rect = TRUE)
# Visualization by Heatmap

# Method 2. Average Linkage Method