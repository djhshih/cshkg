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
