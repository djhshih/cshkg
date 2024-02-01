# Load required libraries for the analysis.
library(io)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggsci)

# find Housekeeping Genes which have low variance among the samples

# import functions to calculate SD (at the same time, calculate mean values)
source("R/sd_function.R")

load("data/m_genes.rds")
load("data/samples.rds")
load("data/result.rds")
load("data/m_ordered.rds")

# # log scale matrix
# m_genes <- log(m_genes + 1)
# sum(is.na(m_genes))
# #[1] 0 no NA

sum(is.na(m_ordered)) # 0

# Calculate within-group and between-group SD
genes <- rownames(m_ordered)

hist(m_ordered, breaks=100)

means <- rowMeans(m_ordered, na.rm=TRUE)
groups <- samples$group

#debug(get_within.sd)
#get_within.sd(m_genes[genes[1], ], groups)

within_SDs <- unlist(lapply(genes, function(gene) {
  message(gene)
  get_within.sd(m_ordered[gene, ], groups)
}))

between_SDs <- unlist(lapply(genes, function(gene) {
  message(gene)
  get_btwn.sd(m_ordered[gene, ], groups)
}))

# Make data frame of SD results
SDs_df2 <- data.frame(
  gene = genes, within_sd = within_SDs, between_sd = between_SDs, mean = means
  )

saveRDS(SDs_df, "out/sds.rds")
saveRDS(SDs_df2, "out/sds2.rds")

# Figure of within vs between SD
N <- ncol(m_ordered)
K <- length(levels(factor(samples$group)))

# Histogram of SDs_df to visualize
hist(SDs_df2$within_sd, breaks = 50, xlim = c(0,1))
hist(SDs_df2$between_sd, breaks = 50, xlim = c(0,1))
hist(SDs_df2$mean, breaks = 50, xlim = c(4,10))

# Find candidates of housekeeping genes.
find_candidates <- data.frame(
  within_sd = quantile(SDs_df2$within_sd, probs = c(.0,.25, .5, .75, 1.0)),
  between_sd = quantile(SDs_df2$between_sd, probs = c(.0,.25, .5, .75, 1.0)),
  mean = quantile(SDs_df2$mean, probs = c(.0,.25, .5, .75, 1.0))
  )
shorlist <- SDs_df2[SDs_df2$within_sd < 0.2937889
                   & SDs_df2$between_sd > 0.33247443 
                   & SDs_df2$mean > 0.0138620, ]

# Join the 'gene' column values
housekeeping <- paste0(shorlist$gene, collapse = '","')
housekeeping <- paste('"', housekeeping, '"', sep = '')
cat(housekeeping)
housekeeping <- c("RNMT","CDH3","EPN2","STX4","CGRRF1","MAPK13","AARS",
                  "TERF2IP","NUDCD3","TMEM109","ICMT","ETFA","ARHGAP1",
                  "TCEAL4","PTPN1","TES","PIN1","DDB2","EPHA2","PPIC",
                  "NCOA3","PTK2","ADGRG1","NENF","DNAJC15","TNFRSF21",
                  "ARPP19","GNAS","PXN","TXNRD1","S100A13","FDFT1","ANXA7",
                  "DHRS7","BAG3","CD320")

hkg <- dplyr::filter(SDs_df2, gene %in% housekeeping)

# SDs of all genes
qdraw(
  ggplot(SDs_df2, aes(x = within_sd, y = between_sd)) + 
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
  ggplot(SDs_df2, aes(x = within_sd, y = mean, label = gene)) +
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
