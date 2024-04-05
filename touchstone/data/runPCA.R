# Load required libraries for running PCA.
library(io)
library(corrr)
library(ggplot2)
library(ggrepel)
library(ggcorrplot)
library(ggpubr)
library(FactoMineR)
library(factoextra)
# Touchstone Dataset Metadata
load("m.rds")
load("m_genes.rds")
load("touchstone_sample.rds")
load("samples.rds")

pca_prcomp <- prcomp(t(scale(m_genes)), scale = FALSE)
head(pca_prcomp$sdev)
summary(pca_prcomp)
p_components <- pca_prcomp$x

# Screeplot
screeplot(pca_prcomp, bstick = TRUE, type = "l", main = NULL)

  ## Check variation of Principal Components in the dataset
  pca_var <- pca_prcomp$sdev^2 # show how much variation each PC accounts for
  pca_var_per <- round(pca_var/sum(pca_var)*100, 1) # represent as percentage %
  barplot(pca_var_per, main = "Scree plot Principal Components Variation", 
          xlab = "Principal Component (PC1~PC50)", ylab = "Percent Variation",
          xlim = c(0,50))
  ## PC1 is regarded as about 28% of variation in the dataset.

# Biploy using ggplot to represent data of PC1 over PC2
pca_prcomp_data <- data.frame(pca_prcomp$x)
pca_prcomp_data$plotx <- pca_prcomp_data[,1]
pca_prcomp_data$ploty <- pca_prcomp_data[,2]
all(rownames(pca_prcomp_data) == samples$inst_id)
pca_prcomp_data$group <- as.character(samples$group)

VCAP_palette <- colorRampPalette(c("mediumpurple", "deepskyblue"))(length(which(grepl("^VCAP", unique(pca_prcomp_data$group)))))
MCF7_palette <- colorRampPalette(c("deeppink", "coral"))(length(which(grepl("^MCF7", unique(pca_prcomp_data$group)))))
PC3_palette <- colorRampPalette(c("orange", "limegreen"))(length(which(grepl("^PC3", unique(pca_prcomp_data$group)))))
pca_color <- c(PC3_palette, VCAP_palette, MCF7_palette)

p_pca_prcomp <- ggplot(pca_prcomp_data, aes(x = plotx, y = ploty, color = group)) +
  geom_point() + labs(x = "PC1", y = "PC2", title = "PCA with sample groups") +
  theme(legend.position = "none", legend.key.size = unit(0.01, "lines")) +
  scale_color_manual(values = pca_color)
  # stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(color, 0.3))),
  #              data = pca_prcomp_data[pca_prcomp_data$group != "VCAP_UnTrt",])
p_pca_prcomp
## RESULT: PCA showed 3 distinct groups depends on cell type.
## Continue PCA for each cell type, "VCAP", "MCF7", and "PC3"

# Current samples are divided depends on cell type and perturbagen (total 2380 groups)
# Divide matrix by cell types "VCAP", "MCF7", "PC3" (total 3 matrices)
samples$cell_type <- sub("_.*", "", samples$group)

s_VCAP <- samples[which(samples$cell_type == "VCAP"),]
m_VCAP <- m_genes[, intersect(colnames(m_genes), s_VCAP$inst_id)]
dim(m_VCAP) # 978 genes x 5078 samples

s_MCF7 <- samples[which(samples$cell_type == "MCF7"),]
m_MCF7 <- m_genes[, intersect(colnames(m_genes), s_MCF7$inst_id)]
dim(m_MCF7) # 978 genes x 21776 samples

s_PC3 <- samples[which(samples$cell_type == "PC3"),]
m_PC3 <- m_genes[, intersect(colnames(m_genes), s_PC3$inst_id)]
dim(m_PC3) # 978 genes x 22353 samples

save(m_VCAP, file = "m_VCAP.rds")
save(s_VCAP, file = "s_VCAP.rds")
save(m_MCF7, file = "m_MCF7.rds")
save(s_MCF7, file = "s_MCF7.rds")
save(m_PC3, file = "m_PC3.rds")
save(s_PC3, file = "s_PC3.rds")

# 1. Cell type == "VCAP"
VCAP_pca_prcomp <- prcomp(t(scale(m_VCAP)), scale = FALSE)
head(VCAP_pca_prcomp$sdev)
summary(VCAP_pca_prcomp)
## Screeplot
screeplot(VCAP_pca_prcomp, bstick = TRUE, type = "l", main = NULL)
### Check variation of Principal Components in the dataset
VCAP_pca_var <- VCAP_pca_prcomp$sdev^2 # show how much variation each PC accounts for
VCAP_pca_var_per <- round(VCAP_pca_var/sum(VCAP_pca_var)*100, 1) # represent as percentage %
barplot(VCAP_pca_var_per, main = "Scree plot Principal Components Variation of VCAP", 
        xlab = "Principal Component (PC1~PC50)", ylab = "Percent Variation",
        xlim = c(0,50))
### PC1 is regarded as about 23.1% of variation in the dataset.

## Biplot
biplot(VCAP_pca_prcomp, scale = 0, xlabs = rep(".", nrow(t(m_VCAP))))
## Biploy using ggplot to represent data of PC1 over PC2
VCAP_pca_prcomp_data <- data.frame(VCAP_pca_prcomp$x)
VCAP_pca_prcomp_data$plotx <- VCAP_pca_prcomp_data[,1]
VCAP_pca_prcomp_data$ploty <- VCAP_pca_prcomp_data[,2]
all(rownames(VCAP_pca_prcomp_data) == s_VCAP$inst_id)
VCAP_pca_prcomp_data$group <- as.character(s_VCAP$group)

p_pca_VCAP <- ggplot(VCAP_pca_prcomp_data, aes(x = plotx, y = ploty, color = group)) +
  geom_point() + labs(x = "PC1", y = "PC2", title = "PCA with sample groups of VCAP") +
  theme(legend.key.size = unit(0.01, "cm"), legend.text = element_text(size = 5),
        legend.position = "none") + 
  scale_colour_ordinal()
  # stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(color, 0.01))),
  #            data = VCAP_pca_prcomp_data[VCAP_pca_prcomp_data$group != "VCAP_UnTrt", ])
p_pca_VCAP

# Look for samples without >= 20 count frequency groups
VCAP_pca_prcomp_data$count_freq <- s_VCAP$count_freq
sub_VCAP_prcomp <- subset(VCAP_pca_prcomp_data, count_freq <= 20)
ggplot(sub_VCAP_prcomp, aes(x = plotx, y = ploty, color = group)) +
  geom_point() + labs(x = "PC1", y = "PC2", 
                      title = "PCA with sample groups of VCAP without high (>=20) frequency groups") +
  theme(legend.key.size = unit(0.01, "cm"), legend.text = element_text(size = 5),
        legend.position = "bottom") + 
  scale_colour_ordinal() + geom_text_repel(aes(label = group)) +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(color, 0.01))),
               data = sub_VCAP_prcomp[sub_VCAP_prcomp$group != "VCAP_UnTrt",])

unique(sub_VCAP_prcomp$group) # all in same groups except "VCAP_NNC-55-0396"
# check
sub_VCAP_prcomp <- sub_VCAP_prcomp[sub_VCAP_prcomp$group != "VCAP_NNC-55-0396", ]
ggplot(sub_VCAP_prcomp, aes(x = plotx, y = ploty, color = group)) +
  geom_point() + labs(x = "PC1", y = "PC2", 
                      title = "PCA with sample groups of VCAP without high (>=30) frequency groups and VCAP_grp1") +
  theme(legend.key.size = unit(0.01, "cm"), legend.text = element_text(size = 5),
        legend.position = "bottom") + 
  scale_colour_ordinal() + geom_text_repel(aes(label = group))
## RESULT: frequency under 23 (total 3 groups) showed none distinct effect.
## RESULT: frequency of 24 (total 51 groups) showed none distinct effect.
## RESULT: frequency above 24 (total 1 group == "BRD-K77432048") showed none distinct effect

# 2. Cell type == "MCF7"
MCF7_pca_prcomp <- prcomp(t(scale(m_MCF7)), scale = FALSE)
head(MCF7_pca_prcomp$sdev)
summary(MCF7_pca_prcomp)

# Screeplot
screeplot(MCF7_pca_prcomp, bstick = TRUE, type = "l", main = NULL)

## Check variation of Principal Components in the dataset
MCF7_pca_var <- MCF7_pca_prcomp$sdev^2 # show how much variation each PC accounts for
MCF7_pca_var_per <- round(MCF7_pca_var/sum(MCF7_pca_var)*100, 1) # represent as percentage %
barplot(MCF7_pca_var_per, main = "Scree plot Principal Components Variation of MCF7", 
        xlab = "Principal Component (PC1~PC50)", ylab = "Percent Variation",
        xlim = c(0,50))
## PC1 is regarded as about 13.4% of variation in the dataset.

# Biplot
biplot(MCF7_pca_prcomp, scale = 0, xlabs = rep(".", nrow(t(m_MCF7))))

# Biploy using ggplot to represent data of PC1 over PC2
MCF7_pca_prcomp_data <- data.frame(MCF7_pca_prcomp$x)
MCF7_pca_prcomp_data$plotx <- MCF7_pca_prcomp_data[,1]
MCF7_pca_prcomp_data$ploty <- MCF7_pca_prcomp_data[,2]
all(rownames(MCF7_pca_prcomp_data) == s_MCF7$inst_id)
MCF7_pca_prcomp_data$group <- as.character(s_MCF7$group)

p_pca_MCF7 <- ggplot(MCF7_pca_prcomp_data, aes(x = plotx, y = ploty, color = group)) +
  geom_point() + labs(x = "PC1", y = "PC2", title = "PCA with sample groups of MCF7") +
  theme(legend.key.size = unit(0.01, "cm"), legend.text = element_text(size = 5),
        legend.position = "none") + 
  scale_colour_ordinal()
  # stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(color, 0.01))),
  #              data = MCF7_pca_prcomp_data[MCF7_pca_prcomp_data$group != "MCF7_UnTrt",])
p_pca_MCF7

# 3. Cell type == "PC3"
PC3_pca_prcomp <- prcomp(t(scale(m_PC3)))
head(PC3_pca_prcomp$sdev)
summary(PC3_pca_prcomp)

# Screeplot
screeplot(PC3_pca_prcomp, bstick = TRUE, type = "l", main = NULL)

## Check variation of Principal Components in the dataset
PC3_pca_var <- PC3_pca_prcomp$sdev^2 # show how much variation each PC accounts for
PC3_pca_var_per <- round(PC3_pca_var/sum(PC3_pca_var)*100, 1) # represent as percentage %
barplot(PC3_pca_var_per, main = "Scree plot Principal Components Variation of PC3", 
        xlab = "Principal Component (PC1~PC50)", ylab = "Percent Variation",
        xlim = c(0,50))
## PC1 is regarded as about 15.5% of variation in the dataset.

# Biplot
biplot(PC3_pca_prcomp, scale = 0, xlabs = rep(".", nrow(t(m_PC3))))

# Biploy using ggplot to represent data of PC1 over PC2
PC3_pca_prcomp_data <- data.frame(PC3_pca_prcomp$x)
PC3_pca_prcomp_data$plotx <- PC3_pca_prcomp_data[,1]
PC3_pca_prcomp_data$ploty <- PC3_pca_prcomp_data[,2]
all(rownames(PC3_pca_prcomp_data) == s_PC3$inst_id)
PC3_pca_prcomp_data$group <- as.character(s_PC3$group)

p_pca_PC3 <- ggplot(PC3_pca_prcomp_data, aes(x = plotx, y = ploty, color = group)) +
  geom_point() + labs(x = "PC1", y = "PC2", title = "PCA with sample groups of PC3") +
  theme(legend.key.size = unit(0.01, "cm"), legend.text = element_text(size = 5),
        legend.position = "none") + 
  scale_colour_ordinal()
  # stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(color, 0.01))),
  #              data = PC3_pca_prcomp_data[PC3_pca_prcomp_data$group != "PC3_UnTrt",])
p_pca_PC3

# Combine plots of PCA for each matrixes
qdraw(
  ggarrange(p_pca_VCAP, p_pca_MCF7, p_pca_PC3, 
          labels = c("A", "B", "C"), ncol = 3, nrow = 1),
  width = 18, height = 5,
  file = "../plots/PCA_result.png"
)

save(VCAP_pca_prcomp_data, file = "VCAP_pca_prcomp_data.rds")
save(MCF7_pca_prcomp_data, file = "MCF7_pca_prcomp_data.rds")
save(PC3_pca_prcomp_data, file = "PC3_pca_prcomp_data.rds")
