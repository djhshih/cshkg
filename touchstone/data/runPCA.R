# Run PCA to select perturbations that give similar effects.
library(corrr)
library(ggcorrplot)
library(FactoMineR)
print("loaded libraries")
# Touchstone Dataset Metadata
load("m.rds")
load("m_genes.rds")
load("touchstone_sample.rds")
load("samples.rds")
print("check dataset loading")

# Current samples are divided depends on cell id and perturbagen (total 2380 groups)
# Divide matrix by cell types "VCAP", "MCF7", "PC3" (total 3 matrixes)
# s_VCAP <- touchstone_sample[which(touchstone_sample$cell_id == "VCAP"),]
# m_VCAP <- m_genes[, intersect(colnames(m_genes), s_VCAP$inst_id)]
# dim(m_VCAP) # 978 genes x 5078 samples
# s_MCF7 <- touchstone_sample[which(touchstone_sample$cell_id == "MCF7"),]
# m_MCF7 <- m_genes[, intersect(colnames(m_genes), s_MCF7$inst_id)]
# dim(m_MCF7) # 978 genes x 21776 samples
# s_PC3 <- touchstone_sample[which(touchstone_sample$cell_id == "PC3"),]
# m_PC3 <- m_genes[, intersect(colnames(m_genes), s_PC3$inst_id)]
# dim(m_PC3) # 978 genes x 22353 samples

pca_result <- prcomp(t(m_genes), scale. = TRUE)
summary(pca_result)

p_components <- pca_result$x
pca_df <- as.data.frame(p_components)
write.table(pca_df, file="pca_results.txt", sep="\t", row.names=FALSE)