# Run PCA to select perturbations that give similar effects.
library(corrr)
library(ggcorrplot)
library(FactoMineR)
# Touchstone Dataset Metadata
load("m.rds")
load("m_genes.rds")
load("touchstone_sample.rds")
load("samples.rds")

# Current samples are divided depends on cell id and perturbagen (total 2380 groups)
corr_matrix <- cor(m_genes)
ggcorrplot(corr_matrix)