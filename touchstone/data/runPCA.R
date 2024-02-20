# Run PCA to select perturbations that give similar effects.
library(corrr)
library(ggplot2)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
# Touchstone Dataset Metadata
load("m.rds")
load("m_genes.rds")
load("touchstone_sample.rds")
load("samples.rds")

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

pca_prcomp <- prcomp(t(scale(m_genes)), scale = FALSE)
head(pca_prcomp$sdev)
summary(pca_prcomp)

p_components <- pca_prcomp$x
# pca_df <- as.data.frame(p_components)
# write.table(pca_df, file="pca_df.txt", sep="\t", row.names=FALSE)

# Screeplot
screeplot(pca_prcomp, bstick = TRUE, type = "l", main = NULL)

  ## Check variation of Principal Components in the dataset
  pca_var <- pca_prcomp$sdev^2 # show how much variation each PC accounts for
  pca_var_per <- round(pca_var/sum(pca_var)*100, 1) # represent as percentage %
  barplot(pca_var_per, main = "Scree plot Principal Components Variation", 
          xlab = "Principal Component (PC1~PC50)", ylab = "Percent Variation",
          xlim = c(0,50))
  ## PC1 is regarded as about 28% of variation in the dataset.

# Biplot
biplot(pca_prcomp, scale = 0, xlabs = rep(".")) # ERROR
fviz_pca_biplot(pca_prcomp, label = "var", 
                habillage = samples$group, repel = TRUE) # vector memotry exhausted

# ggplot to represent data
pca_prcomp_data <- data.frame(pca_prcomp$x)
pca_prcomp_data$plotx <- pca_prcomp_data[,1]
pca_prcomp_data$ploty <- pca_prcomp_data[,2]
all(rownames(pca_prcomp_data) == samples$inst_id)
pca_prcomp_data$group <- as.character(samples$group)

p_pca_prcomp <- ggplot(pca_prcomp_data, aes(x = plotx, y = ploty, color = group)) +
  geom_point() + labs(x = "PC1", y = "PC2", title = "PCA with sample groups") +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(color, 0.3))), 
               data = pca_prcomp_data[pca_prcomp_data$group,])
p_pca_prcomp
