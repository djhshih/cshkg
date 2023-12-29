library(cmapR)
library(tidyr)
library(dplyr)
my_ds <- parse_gctx("GSE92742_Broad_LINCS_Level2_GEX_delta_n49216x978.gctx")
# sample info
GSE92742_inst_info <- read.delim("GSE92742_Broad_LINCS_inst_info.txt")
# gene info
GSE92742_gene_info_delta_landmark <- read.delim("GSE92742_Broad_LINCS_gene_info_delta_landmark.txt")

# Touchstone Dataset Metadata
touchstone_data <- my_ds@mat
touchstone_sample <- my_ds@cdesc

# Match sample id
colnames(touchstone_sample)[1] <- "inst_id"
touchstone_sample <- merge(touchstone_sample,GSE92742_inst_info, by = "inst_id")
touchstone_sample <- touchstone_sample[,c(1,4:11)]
# Make sample group names by cell id, perturbagen
sample_group <- unite(touchstone_sample, group, cell_id, pert_iname, pert_time, sep = "_")
sample_group <- sample_group[,c(1,3)]
# Count number of samples in each group
sample_group <- sample_group %>% 
  group_by(group) %>% 
  mutate(count_freq = n())
hist(sample_group$count_freq, breaks = 50, xlim = c(0,50)) # filter sample group frequency less than 10

# Filter sample group
filter_sample_group <- subset(sample_group, count_freq >= 10)
length(unique(filter_sample_group$group))
## [1] 366: Total 366 groups

# Make group information df with group names and number of samples in each group
sample_group_info <- filter_sample_group[,c(2,3)]
sample_group_info <- unique(sample_group_info[c("group", "count_freq")])

# Filter samples in dataset
touchstone_df <- touchstone_data[, intersect(colnames(touchstone_data), filter_sample_group$inst_id)]

# Match gene id with gene symbol
all(rownames(touchstone_df) %in% GSE92742_gene_info_delta_landmark$pr_gene_id) #TRUE
all(rownames(touchstone_df) == GSE92742_gene_info_delta_landmark$pr_gene_id) #FALSE
## sort genes by order
sort_gene <- rownames(touchstone_df)
GSE92742_gene_info_delta_landmark <- GSE92742_gene_info_delta_landmark[order(match(GSE92742_gene_info_delta_landmark$pr_gene_id, sort_gene)),]
all(rownames(touchstone_df) == GSE92742_gene_info_delta_landmark$pr_gene_id) #TRUE
## Convert gene id to gene symbol
rownames(touchstone_df) <- GSE92742_gene_info_delta_landmark$pr_gene_symbol
touchstone_df <- as.matrix(touchstone_df)
all(rownames(touchstone_df) == GSE92742_gene_info_delta_landmark$pr_gene_symbol) #TRUE

save(touchstone_df, file = "touchstone_df.rds")
save(filter_sample_group, file = "filter_sample_group.rds")

# Group Clustering
library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
# Touchstone Dataset Metadata
df <- my_ds@mat
df_sample <- my_ds@cdesc
# Match sample id
colnames(df_sample)[1] <- "inst_id"
df_sample <- merge(df_sample,GSE92742_inst_info, by = "inst_id")
df_sample <- df_sample[,c(1,4:11)]
# Make sample group names by cell id, perturbagen
df_group <- unite(df_sample, group, cell_id, pert_iname, sep = "_")
df_group <- df_group[,c(1,3)]
# Count number of samples in each group
df_group <- df_group %>% 
  group_by(group) %>% 
  mutate(count_freq = n())

# group clustering
library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
# Perform hierarchical clustering on the transposed dataset (samples in rows, genes in columns)
df_rev <- na.omit(t(df))
#d <- dist(df_rev)
#hc <- hclust(d)
subset_size <- 100
num_columns <- ncol(df_rev)
# Ensure that the subset size is not larger than the number of columns
if (subset_size > num_columns) {
  stop("Subset size is larger than the number of columns in the dataset.")
}
# Create a subset df
subset_df <- df_rev[, sample(num_columns, size = subset_size, replace = FALSE)]
install.packages("flashClust")
library(flashClust)
# Perform hierarchical clustering using flashClust
hc <- flashClust(dist(subset_df))
#hc <- hclust(dist(smaller_subset_df))

# Perform hierarchical clustering on the transposed dataset
hc <- agnes(df_rev, method = "euclidean")

# Plot the dendrogram
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = 4)
# Visualize the dendrogram
plot(dend, main = "Hierarchical Clustering Dendrogram", xlab = "Samples", sub = "")

# Cut the tree to get clusters
clusters <- cutree(hc, h = 50)  # Adjust the height (h) based on the dendrogram
dim(df_group)
length(clusters)
# Combine the cluster information with your sample grouping data
result_df <- cbind(df_group, Cluster = clusters)

