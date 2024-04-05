# Load required libraries for the analysis.
library(io)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggsci)

# find Housekeeping Genes which have low variance among the samples

# import functions to calculate SD (at the same time, calculate mean values)
source("R/sd_function.R")

load("data/m_genes.rds")
load("data/samples.rds")
load("data/samples_cluster.rds")

# check NA
sum(is.na(m_genes)) #[1] 0 no NA

# Calculate within-group and between-group SD
genes <- rownames(m_genes)

hist(m_genes, breaks=100)

means <- rowMeans(m_genes, na.rm=TRUE)
groups <- samples_cluster$louvain_group_name

#debug(get_within.sd)
#get_within.sd(m_genes[genes[1], ], groups)

within_SDs <- unlist(lapply(genes, function(gene) {
  message(gene)
  get_within.sd(m_genes[gene, ], groups)
}))

between_SDs <- unlist(lapply(genes, function(gene) {
  message(gene)
  get_btwn.sd(m_genes[gene, ], groups)
}))

# Make data frame of SD results
SDs_df <- data.frame(
  gene = genes, within_sd = within_SDs, between_sd = between_SDs, mean = means
  )

save(SDs_df, file = "out/sds.rds")

# Figure of within vs between SD
N <- ncol(m_genes)
K <- length(levels(factor(samples_cluster$louvain_group_name)))

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
shorlist <- SDs_df[SDs_df$within_sd < 0.2234330
                   & SDs_df$between_sd > 0.4394126 
                   & SDs_df$mean > 7.866634, ]

# Join the 'gene' column values
housekeeping <- paste0(shorlist$gene, collapse = '","')
housekeeping <- paste('"', housekeeping, '"', sep = '')
cat(housekeeping)
housekeeping <- c("XBP1","GRN","PXN","EBNA1BP2","LRP10","VPS26A","SEC24C",
                  "TCEAL4","S100A13","EPHB4","ERCC1","SNRPF","PPIC","CLTB",
                  "NCOA3","TSPAN6","ABAT","ANXA7","DHRS7","BAG3","DERA","CD320",
                  "CISD1","TNFRSF21","TMEM127","ZGPAT")

hkg <- dplyr::filter(SDs_df, gene %in% housekeeping)

save(hkg, file = "out/candidate_cshkg.rds")


temp <- representative[representative$louvain_group == "40_PC3",]
temp <- temp |> separate(group, into = c("cell_line", "perturbagen"), sep = "_", extra = "merge", fill = "right")
temp2 <- paste0(temp$perturbagen, collapse = ', ')
cat(temp2)

# Mark for common housekeeping genes among 978 landmarks
## Examples of well-established housekeeping genes encompass followings:
## ACTB, ATP5PB, B2M, GAPDH, HPRT1, PGK1, RER1, RPL13A, RPL27, SDHA, TBP, UBC,
## C1orf43, CHMP2A, EMC7, GPI, PSMB2, PSMB4, RAB7A, REEP5, SNRPD3, VCP, VPS29
common_hkg <- c("ACTB", "ATP5PB", "B2M", "GAPDH", "HPRT1", "PGK1", "RER1", 
                "RPL13A", "RPL27", "SDHA", "TBP", "UBC", "C1orf43", "CHMP2A", 
                "EMC7", "GPI", "PSMB2", "PSMB4", "RAB7A", "REEP5", "SNRPD3", 
                "VCP", "VPS29")
common_hkg <- dplyr::filter(SDs_df, gene %in% common_hkg)

save(common_hkg, file = "out/common_hkg.rds")

# SDs of all genes
## Candidates of cshkgs are labeled.
qdraw(
  ggplot(SDs_df, aes(x = within_sd, y = between_sd)) + 
    ggtitle("Between-group SD against Within-group SD of all genes") +
    geom_point(alpha = 0.3, size = 0.5) + xlim(0, 1) + ylim(0, 1.00350240) +
    geom_point(data = hkg, aes(x = within_sd, y = between_sd), color = "blue", size = 0.5) +
    geom_label_repel(data = hkg, aes(label=gene), color="blue", label.padding = 0.1, 
                     max.overlaps = Inf, segment.curvature = -1e-20, 
                     segment.square = TRUE) + 
    geom_point(data = common_hkg, aes(x = within_sd, y = between_sd), color = "red", size = 0.5) +
    geom_label_repel(data = common_hkg, aes(label=gene), color = "red", label.padding = 0.1, 
                     max.overlaps = Inf, segment.curvature = -1e-20, 
                     segment.square = TRUE) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme_linedraw() + 
    geom_abline(intercept = 0, slope = sqrt((N - K) / (N - K - 2)),
                color="grey"),
  width = 7, height = 7,
  file = "plots/SDs_all_gene.png"
)

# Mean expression of all genes against within-group SD.
## Candidates of cshkgs and common housekeeping genes among landmark genes are labeled.
qdraw(
  ggplot(SDs_df, aes(x = within_sd, y = mean, label = gene)) +
    ggtitle("Mean expression of all genes against Within-group SD") +
    geom_point(alpha = 0.3, size = 0.5) + xlim(0, 1) + ylim(4.5, 9.5) +
    geom_point(data = hkg, aes(x = within_sd, y = mean), color = "blue", size = 0.5) +
    geom_label_repel(data = hkg, color = "blue", label.padding = 0.1, 
                     max.overlaps = Inf, segment.curvature = -1e-20, 
                     segment.square = TRUE) +
    geom_point(data = common_hkg, aes(x = within_sd, y = mean), color = "red", size = 0.5) +
    geom_label_repel(data = common_hkg, color = "red", label.padding = 0.1, 
                     max.overlaps = Inf, segment.curvature = -1e-20, 
                     segment.square = TRUE) +
    theme_linedraw(),
  width = 7, height = 7,
  file = "plots/mean_expression_all_gene.png"
)

# ------------------------------------------------------------------------------
# Find list of potential universal hkgs that has lowest within-group SD

shorlist <- top_n(SDs_df, -10, within_sd)
shorlist <- shorlist[shorlist$between_sd < 0.3239494, ]
housekeeping <- paste0(shorlist$gene, collapse = '","')
housekeeping <- paste('"', housekeeping, '"', sep = '')
cat(housekeeping)
housekeeping <- c("RNPS1","OXSR1","FAH","C2CD5","VDAC1")
uni_hkg <- dplyr::filter(SDs_df, gene %in% housekeeping)

# Potential universal hkgs and known common hkgs are labeled
# 1. Within-group SD vs Between-group SD
qdraw(
  ggplot(SDs_df, aes(x = within_sd, y = between_sd)) + 
    ggtitle("Between-group SD against Within-group SD of all genes") +
    geom_point(alpha = 0.3, size = 0.5) + xlim(0, 1) + ylim(0, 1.00350240) +
    geom_point(data = uni_hkg, aes(x = within_sd, y = between_sd), color = "limegreen", size = 0.5) +
    geom_label_repel(data = uni_hkg, aes(label=gene), color="limegreen", label.padding = 0.1, 
                     max.overlaps = Inf, segment.curvature = -1e-20, 
                     segment.square = TRUE) + 
    geom_point(data = common_hkg, aes(x = within_sd, y = between_sd), color = "red", size = 0.5) +
    geom_label_repel(data = common_hkg, aes(label=gene), color = "red", label.padding = 0.1, 
                     max.overlaps = Inf, segment.curvature = -1e-20, 
                     segment.square = TRUE, position = position_nudge(x = 0.04, y = 0.025)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_linedraw() + 
    geom_abline(intercept = 0, slope = sqrt((N - K) / (N - K - 2)),
                color="grey"),
  width = 7, height = 7,
  file = "plots/SDs_all_gene_ver2.png"
)

# 2. Within-group SD vs Mean expression levels
qdraw(
  ggplot(SDs_df, aes(x = within_sd, y = mean, label = gene)) +
    ggtitle("Mean expression of all genes against Within-group SD") +
    geom_point(alpha = 0.3, size = 0.5) + xlim(0, 1) + ylim(4.5, 9.5) +
    geom_point(data = uni_hkg, aes(x = within_sd, y = mean), color = "limegreen", size = 0.5) +
    geom_label_repel(data = uni_hkg, color = "limegreen", label.padding = 0.1, 
                     max.overlaps = Inf, segment.curvature = -1e-20, 
                     segment.square = TRUE) +
    geom_point(data = common_hkg, aes(x = within_sd, y = mean), color = "red", size = 0.5) +
    geom_label_repel(data = common_hkg, color = "red", label.padding = 0.1, 
                     max.overlaps = Inf, segment.curvature = -1e-20, 
                     segment.square = TRUE) +
    theme_linedraw(),
  width = 7, height = 7,
  file = "plots/mean_expression_all_gene_ver2.png"
)

# ------------------------------------------------------------------------------
# Examine for each housekeeping genes

hkgplot <- function(d, mean) {
  ggplot(d, aes(x = inst_id, y = expression)) +
    facet_wrap(~ group, scales="free_x") +
    geom_point(size=0.5) + theme_classic() +
    geom_hline(yintercept = expression_mean, colour = "steelblue") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank(),
      panel.grid.major.y = element_line()
    ) +
    xlab("sample") + ylab("gene expression") +
    ggtitle(candidate)
}

# Select groups that has more than 100 samples
selected_group <- samples[samples$count_freq > 100, ]
selected_group <- paste0(unique(selected_group$group), collapse = '","')
selected_group <- paste('"', selected_group, '"', sep = '')
cat(selected_group)
selected_group <- c("VCAP_UnTrt","VCAP_GFP","VCAP_ERG","MCF7_GFP","MCF7_lacZ",
                    "MCF7_EMPTY_VECTOR","MCF7_LUCIFERASE","MCF7_UnTrt","MCF7_RFP",
                    "MCF7_pgw","PC3_GFP","PC3_lacZ","PC3_EMPTY_VECTOR",
                    "PC3_LUCIFERASE","PC3_UnTrt","PC3_RFP","PC3_pgw")
## total 17 groups

# 1. `HPRT1`: known hkg
candidate <- "HPRT1"

expression_d <- data.frame(samples_cluster, expression = m_genes[candidate, ])
expression_d_filtered <- subset(expression_d, group %in% selected_group)
expression_mean <- mean(expression_d$expression)

qdraw(
  hkgplot(expression_d_filtered, expression_mean) +
    labs(subtitle = sprintf("between-group (all groups) SD: %.3f, between-group (filtered) SD: %.3f",
                            with(expression_d, get_btwn.sd(expression, group)),
                            with(expression_d_filtered, get_btwn.sd(expression, group)))),
  width = 9, height = 5,
  file = "plots/HPRT1_known.png"
)

# 2. `REEP5`: known hkg
candidate <- "REEP5"

expression_d <- data.frame(samples_cluster, expression = m_genes[candidate, ])
expression_d_filtered <- subset(expression_d, group %in% selected_group)
expression_mean <- mean(expression_d$expression)

qdraw(
  hkgplot(expression_d_filtered, expression_mean) +
    labs(subtitle = sprintf("between-group (all groups) SD: %.3f, between-group (filtered) SD: %.3f",
                            with(expression_d, get_btwn.sd(expression, group)),
                            with(expression_d_filtered, get_btwn.sd(expression, group)))),
  width = 9, height = 5,
  file = "plots/REEP5_known.png"
)

# 3. `SNRPF`: candidate cshkg
candidate <- "SNRPF"

expression_d <- data.frame(samples_cluster, expression = m_genes[candidate, ])
expression_d_filtered <- subset(expression_d, group %in% selected_group)
expression_mean <- mean(expression_d$expression)

qdraw(
  hkgplot(expression_d_filtered, expression_mean) +
    labs(subtitle = sprintf("between-group (all groups) SD: %.3f, between-group (filtered) SD: %.3f",
                            with(expression_d, get_btwn.sd(expression, group)),
                            with(expression_d_filtered, get_btwn.sd(expression, group)))),
  width = 9, height = 5,
  file = "plots/SNRPF_candidate.png"
)

# 4. `VPS26A`: candidate cshkg
candidate <- "VPS26A"

expression_d <- data.frame(samples_cluster, expression = m_genes[candidate, ])
expression_d_filtered <- subset(expression_d, group %in% selected_group)
expression_mean <- mean(expression_d$expression)

qdraw(
  hkgplot(expression_d_filtered, expression_mean) +
    labs(subtitle = sprintf("between-group (all groups) SD: %.3f, between-group (filtered) SD: %.3f",
                            with(expression_d, get_btwn.sd(expression, group)),
                            with(expression_d_filtered, get_btwn.sd(expression, group)))),
  width = 9, height = 5,
  file = "plots/VPS26A_candidate.png"
)


# Select groups that has more than 1000 samples
samples_cluster <- samples_cluster[, c(1,4,5,6)]
colnames(samples_cluster)[2] <- "group"
samples_cluster <- samples_cluster %>%
  group_by(group) %>%
  mutate(count_freq = n())
selected_group <- samples_cluster[samples_cluster$count_freq > 500, ]
selected_group <- paste0(unique(selected_group$group), collapse = '","')
selected_group <- paste('"', selected_group, '"', sep = '')
cat(selected_group)
selected_group <- c("PC3_UnTrt_1","MCF7_UnTrt_12","PC3_UnTrt_12","MCF7_LUCIFERASE_15",
                    "PC3_UnTrt_15","MCF7_UnTrt_16","MCF7_UnTrt_18","MCF7_GFP_27",
                    "PC3_UnTrt_27","MCF7_UnTrt_3","MCF7_UnTrt_30","MCF7_LUCIFERASE_31",
                    "PC3_UnTrt_32","MCF7_pgw_33","PC3_UnTrt_34","MCF7_UnTrt_35",
                    "PC3_LUCIFERASE_36","PC3_UnTrt_38","PC3_UnTrt_39","PC3_UnTrt_40",
                    "PC3_UnTrt_41","MCF7_GFP_43","MCF7_UnTrt_45","PC3_UnTrt_47",
                    "MCF7_UnTrt_48","PC3_KRAS_50","MCF7_UnTrt_9")
## total 27 groups

# 1. `HPRT1`: known hkg
candidate <- "HPRT1"

expression_d <- data.frame(samples_cluster, expression = m_genes[candidate, ])
expression_d_filtered <- subset(expression_d, group %in% selected_group)
expression_mean <- mean(expression_d$expression)

qdraw(
  hkgplot(expression_d_filtered, expression_mean) +
    labs(subtitle = sprintf("between-group (all groups) SD: %.3f, between-group (filtered) SD: %.3f",
                            with(expression_d, get_btwn.sd(expression, group)),
                            with(expression_d_filtered, get_btwn.sd(expression, group)))),
  width = 9, height = 5,
  file = "plots/HPRT1_known_louv.png"
)

# 2. `REEP5`: known hkg
candidate <- "REEP5"

expression_d <- data.frame(samples_cluster, expression = m_genes[candidate, ])
expression_d_filtered <- subset(expression_d, group %in% selected_group)
expression_mean <- mean(expression_d$expression)

qdraw(
  hkgplot(expression_d_filtered, expression_mean) +
    labs(subtitle = sprintf("between-group (all groups) SD: %.3f, between-group (filtered) SD: %.3f",
                            with(expression_d, get_btwn.sd(expression, group)),
                            with(expression_d_filtered, get_btwn.sd(expression, group)))),
  width = 9, height = 5,
  file = "plots/REEP5_known_louv.png"
)

# 3. `SNRPF`: candidate cshkg
candidate <- "SNRPF"

expression_d <- data.frame(samples_cluster, expression = m_genes[candidate, ])
expression_d_filtered <- subset(expression_d, group %in% selected_group)
expression_mean <- mean(expression_d$expression)

qdraw(
  hkgplot(expression_d_filtered, expression_mean) +
    labs(subtitle = sprintf("between-group (all groups) SD: %.3f, between-group (filtered) SD: %.3f",
                            with(expression_d, get_btwn.sd(expression, group)),
                            with(expression_d_filtered, get_btwn.sd(expression, group)))),
  width = 9, height = 5,
  file = "plots/SNRPF_candidate_louv.png"
)

# 4. `VPS26A`: candidate cshkg
candidate <- "VPS26A"

expression_d <- data.frame(samples_cluster, expression = m_genes[candidate, ])
expression_d_filtered <- subset(expression_d, group %in% selected_group)
expression_mean <- mean(expression_d$expression)

qdraw(
  hkgplot(expression_d_filtered, expression_mean) +
    labs(subtitle = sprintf("between-group (all groups) SD: %.3f, between-group (filtered) SD: %.3f",
                            with(expression_d, get_btwn.sd(expression, group)),
                            with(expression_d_filtered, get_btwn.sd(expression, group)))),
  width = 9, height = 5,
  file = "plots/VPS26A_candidate_louv.png"
)
