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
m_sub <- m_genes[,1:40000]
transposed_m_genes <- t(m_genes)
pca_result <- prcomp(transposed_m_genes, scale. = TRUE)
summary(pca_result)

p_components <- pca_result$x
pca_df <- as.data.frame(p_components)
write.table(pca_df, file="pca_results.txt", sep="\t", row.names=FALSE)