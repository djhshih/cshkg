# Load required libraries for the analysis.
library(io)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggsci)

# find Housekeeping Genes which have low variance among the samples

# import functions to calculate SD (at the same time, calculate mean values)
source("R/sd_function.R")

load("touchstone_df.RData")
load("filter_sample_group.RData")

# Calculate within-group and between-group SD
genes <- rownames(touchstone_df)

means <- rowMeans(touchstone_df, na.rm=TRUE)
groups <- filter_sample_group$group

within_SDs <- unlist(lapply(genes, function(gene) {
  get_within.sd(touchstone_df[gene, ], groups)
}))

between_SDs <- unlist(lapply(genes, function(gene) {
  get_btwn.sd(mat[gene, ], groups)
}))

# Make data frame of SD results
SDs_df <- data.frame(
  gene = genes, within_sd = within_SDs, between_sd = between_SDs, mean = means
  )

# Figure of within vs between SD

