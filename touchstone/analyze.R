# Load required libraries for the analysis.
library(io)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggsci)

# find Housekeeping Genes which have low variance among the samples

# import functions to calculate SD (at the same time, calculate mean values)
source("R/sd_function.R")

load("data/touchstone_df.rds")
load("data/filter_sample_group.rds")

# log scale matrix
touchstone_df <- log(touchstone_df + 1)
sum(is.na(touchstone_df))
#[1] 0 no NA

# Calculate within-group and between-group SD
genes <- rownames(touchstone_df)

hist(touchstone_df, breaks=100)

means <- rowMeans(touchstone_df, na.rm=TRUE)
groups <- filter_sample_group$group

#debug(get_within.sd)
#get_within.sd(touchstone_df[genes[1], ], groups)

within_SDs <- unlist(lapply(genes, function(gene) {
  message(gene)
  get_within.sd(touchstone_df[gene, ], groups)
}))

between_SDs <- unlist(lapply(genes, function(gene) {
  message(gene)
  get_btwn.sd(touchstone_df[gene, ], groups)
}))

# Make data frame of SD results
SDs_df <- data.frame(
  gene = genes, within_sd = within_SDs, between_sd = between_SDs, mean = means
  )

saveRDS(SDs_df, "out/sds.rds")

# Figure of within vs between SD
N <- ncol(touchstone_df)
K <- length(levels(factor(filter_sample_group$group)))

# Histogram of SDs_df to visualize
hist(SDs_df$within_sd, breaks = 50, xlim = c(0,1))
hist(SDs_df$between_sd, breaks = 50, xlim = c(0,1))
hist(SDs_df$mean, breaks = 50, xlim = c(4,10))

# Find candidates of housekeeping genes.
find_candidates <- data.frame(
  within_sd = quantile(SDs_df$within_sd, probs = c(.0,.25, .5, .75, 1.0)),
  between_sd = quantile(SDs_df$between_sd, probs = c(.0,.25, .5, .75, 1.0)),
  mean = quantile(SDs_df$mean, probs = c(.0,.25, .5, .75, 1.0))
  )
shorlist <- SDs_df[SDs_df$within_sd < 0.3282383 
                   & SDs_df$between_sd > 0.2773203 
                   & SDs_df$mean > 7.786938, ]

# Join the 'gene' column values
housekeeping <- paste0(shorlist$gene, collapse = '","')
housekeeping <- paste('"', housekeeping, '"', sep = '')
cat(housekeeping)
housekeeping <- c("CKB","KTN1","AARS","PXN","PGRMC1","TARS","LRP10","PNP",
                  "UBE2A","HSD17B10","MYC","DECR1","S100A13","POLR2K","RNH1",
                  "NCOA3","ADH5","MFSD10","DHRS7","RBM34","IARS2","HEATR1")

hkg <- dplyr::filter(SDs_df, gene %in% housekeeping)

# SDs of all genes
qdraw(
  ggplot(SDs_df, aes(x = within_sd, y = between_sd)) + 
    ggtitle("Between group SD against Within group SD of all genes") +
    geom_point(alpha = 0.3, size = 0.5) + xlim(0, 1) + ylim(0, 1.00350240) +
    geom_point(data = hkg, color = "blue", size = 0.5) +
    geom_label_repel(data = hkg, aes(label=gene), color="blue", 
                     label.padding = 0.1, max.overlaps = Inf, 
                     segment.curvature = -1e-20, segment.square = TRUE) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme_linedraw() + 
    geom_abline(intercept = 0, slope = sqrt((N - K) / (N - K - 2)),
                color="grey"),
  width = 7, height = 7,
  file = "plots/SDs_all_gene.png"
)

# Mean expression of all genes against within-group SD.
qdraw(
  ggplot(SDs_df, aes(x = within_sd, y = mean, label = gene)) +
    ggtitle("Mean expression of all genes against Within group SD") +
    geom_point(alpha = 0.3, size = 0.5) + xlim(0, 1) + ylim(4, 10) +
    geom_point(data = hkg, aes(x = within_sd, y = mean), color = "blue", size = 0.5) +
    geom_label_repel(data = hkg, color="blue", label.padding = 0.1, 
                     max.overlaps = Inf, segment.curvature = -1e-20, 
                     segment.square = TRUE) +
    theme_linedraw(),
  width = 7, height = 7,
  file = "plots/mean_expression_all_gene.png"
)
