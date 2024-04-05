# Load required libraries for pre-processing.
library(cmapR)
library(tidyr)
library(dplyr)
my_ds <- parse_gctx("GSE92742_Broad_LINCS_Level2_GEX_delta_n49216x978.gctx")
# sample info
inst_info <- read.delim("GSE92742_Broad_LINCS_inst_info.txt")
# gene info
gene_info_delta_landmark <- read.delim("GSE92742_Broad_LINCS_gene_info_delta_landmark.txt")

# Touchstone Dataset Metadata
m <- my_ds@mat
touchstone_sample <- my_ds@cdesc

# Match sample id
colnames(touchstone_sample)[1] <- "inst_id"
touchstone_sample <- merge(touchstone_sample,inst_info, by = "inst_id")
touchstone_sample <- touchstone_sample[,c(1,4:11)]
# Make sample group names by cell id, perturbagen
sample_group <- unite(touchstone_sample, group, cell_id, pert_iname, sep = "_")
sample_group <- sample_group[,c(1,3)]
# Count number of samples in each group
sample_group <- sample_group %>%
  group_by(group) %>%
  mutate(count_freq = n())
hist(sample_group$count_freq, breaks = 50, xlim = c(0,50)) # filter sample group frequency less than 10

# Filter sample group
samples <- subset(sample_group, count_freq >= 10)
length(unique(samples$group))
## [1] 2380: Total 2380 groups

# Make group information df with group names and number of samples in each group
sample_group_info <- samples[,c(2,3)]
sample_group_info <- unique(sample_group_info[c("group", "count_freq")])

# Filter samples in dataset
m_genes <- m[, intersect(colnames(m), samples$inst_id)]

# Match gene id with gene symbol
all(rownames(m_genes) %in% gene_info_delta_landmark$pr_gene_id) #TRUE
all(rownames(m_genes) == gene_info_delta_landmark$pr_gene_id) #FALSE
## sort genes by order
sort_gene <- rownames(m_genes)
gene_info_delta_landmark <- gene_info_delta_landmark[order(match(gene_info_delta_landmark$pr_gene_id, sort_gene)),]
all(rownames(m_genes) == gene_info_delta_landmark$pr_gene_id) #TRUE
## Convert gene id to gene symbol
rownames(m_genes) <- gene_info_delta_landmark$pr_gene_symbol
## Data in log scale
m_genes <- as.matrix(log(m_genes + 1))
all(rownames(m_genes) == gene_info_delta_landmark$pr_gene_symbol) #TRUE

save(m, file = "m.rds")
save(touchstone_sample, file = "touchstone_sample.rds")
save(m_genes, file = "m_genes.rds")
save(samples, file = "samples.rds")
