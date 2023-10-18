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
K <- length(levels(filter_sample_group$group))

housekeeping <- c("ACTB", "UBC", "GAPDH", "TBP", "RPS18",
                  "G6PD", "HPRT1", "LDHA", "RPLP1", "RPL19",
                  "RPL18", "RPL11", "RPL32", "PGK1", "PPIA",
                  "SDHA", "ASNS", "ATP2B4", "PEX19", "RXRA",
                  "RPL13A")
hkg <- dplyr::filter(SDs_df, gene %in% housekeeping)

qdraw(
  ggplot(SDs_df, aes(x = within_sd, y = between_sd),
    ggtitle("Between group SD against Within group SD") +
    geom_point(alpha = 0.,size = 0.5) +
    geom_point(data = hkg, color = "blue"))
    #theme(plot.title = element_text(hjust = 0.5))+
    #theme_linedraw() +
    #  geom_abline(intercept = 0, slope = sqrt((N - K) / (N - K - 2)), color = "forestgreen"),
    #width = 8, height = 8)
)

