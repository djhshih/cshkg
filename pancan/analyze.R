library(io)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggsci)

#import functions to calculate standard deviation
source("R/sd_function.R")

proj.info <- qread("annot/sample-project-info_pancan.tsv")
type.info <- qread("annot/sample-type-code_tcga.tsv")

x <- qread("data/expr_pancan.rds")
mat0 <- x$data;
mat <- mat0

# Grouping the tissues
idx <- match(substring(colnames(mat),1,16), proj.info$sample)
groups <- proj.info$disease_code[idx];

type_idx <- match(as.numeric(substring(colnames(mat),14,15)), type.info$code)
types <- type.info$letter_code[type_idx]

pheno <- data.frame(
  sample_id = colnames(mat),
  cancer_type = groups,
  sample_type = types
);

pheno$group <- factor(paste0(pheno$cancer_type, "-", pheno$sample_type));

levels(pheno$group)
no.groups.before <- length(levels(pheno$group)) # check the number of groups before filtering out the rare groups

#filter out the rare groups
rare.groups <- names(which(table(pheno$group) < 10));
levels(pheno$group)[levels(pheno$group) %in% rare.groups] <- NA;
idx.remove <- is.na(pheno$group);
table(idx.remove)  # 87 samples that belong to rare groups were set to NA in the mat

levels(pheno$group)
no.groups.after <- length(levels(pheno$group)) # check the number of group after filtering out the rare groups

table(pheno$group)
sum(is.na(pheno$group))

plot_1_df <- data.frame("filter" = factor(c("before", "after"), levels = c("before", "after")), "number_of_group" = c(no.groups.before, no.groups.after))

ggplot(data = plot_1_df, aes(x = filter, y = number_of_group))+ 
  ggtitle("Number of groups before and after filtering")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("filtering") + ylab("number of groups")+
  geom_bar(position="stack",stat = "identity") + 
  theme_minimal()+ 
  geom_text(aes(label= number_of_group), vjust=-1)



stopifnot(pheno$sample_id == colnames(mat))

mat <- mat[, !idx.remove]; # 87 samples that belong to rare groups were set to NA in the mat
pheno <- pheno[!idx.remove, ];

N <- ncol(mat);
K <- length(levels(pheno$group));


#Calculate standard deviation within and between the groups
genes <- rownames(mat)

n.na <- apply(mat, 1, function(z) sum(is.na(z)));

means <- rowMeans(mat, na.rm=TRUE);

within_sds <- unlist(lapply(genes, function(gene) {
  get_within.sd(mat[gene, ], pheno$group)
}));

between_sds <- unlist(lapply(genes, function(gene) {
  get_btwn.sd(mat[gene, ], pheno$group)
}));

d.sd <- data.frame(
  gene = genes, within_sd = within_sds, between_sd = between_sds, mean = means
);

qwrite(d.sd, "out/sds.rds");
#d.sd <- qread("out/sds.rds") # load d.sd 

# SSX9P has many NA and some groups have only NA
# for groups with all missing value, the mean is NaN
v.problematic <- mat["SSX9P", ];
table(is.na(v.problematic))
gmean.problematic <- tapply(v.problematic, pheno$group, mean, na.rm=TRUE);
gmean.problematic
overall_mean <- mean(gmean.problematic, na.rm=TRUE) #this makes sure overall mean is not NaN
#print(sum((group_mean-overall_mean)^2))
between_sd.problematic <- sqrt(sum((gmean.problematic-overall_mean)^2, na.rm=TRUE) / (length(gmean.problematic)-1))


# highlight common housekeeping genes and common genes
housekeeping <- c("ACTB", "UBC", "GAPDH", "TBP", "RPS18", 
                  "G6PD", "HPRT1", "LDHA", "RPLP1","RPL19",
                  "RPL18", "RPL11","RPL32", "PGK1", "PPIA", 
                  "SDHA","ASNS","ATP2B4", "PEX19","RXRA",
                  "RPL13A")  # TODO: Add more!
common_gene <- c("TP53","EGFR","AKT1")

hkg.sub <- dplyr::filter(d.sd, gene %in% housekeeping)
common.sub <- dplyr::filter(d.sd, gene %in% common_gene)

hkg.sub
common.sub

# Set ggplot options globally
options(ggrepel.max.overlaps = Inf)

# figure2: sd of all genes (btwn against within)
qdraw(
  ggplot(d.sd, aes(x=within_sd, y=between_sd))+
    ggtitle("Between group SD against Within group SD")+
    geom_point(alpha = 0., size =0.5)+
    geom_point(data= hkg.sub, color="red")+
    geom_label_repel(data = hkg.sub, aes(label=gene), color="red", label.padding = 0.1, max.overlaps = Inf, segment.curvature=-1e-20, segment.square = TRUE)+
    theme(plot.title = element_text(hjust = 0.5))+
    theme_linedraw()+
    geom_abline(intercept=0, slope=sqrt((N - K) / (N - K - 2)),
      color="grey60"),
  width = 8, height = 8,
  file = "plots/sd_all_genes_with_hkg_highlights2.png"
)  

d.sd.sub <- subset(d.sd, within_sd < 1)

# mean expression of d.sd and d.sd.sub
qdraw(
  ggplot(d.sd, aes(x=within_sd, y=mean, label=gene))+
    #geom_text_repel() +
    geom_point(alpha=0.1)+
    geom_point(data= hkg.sub, aes(x=within_sd, y=mean), color="red")+
    geom_label_repel(data = hkg.sub, color="red",label.padding = 0.1, segment.curvature=-1e-20, segment.square = TRUE)+
    theme_minimal()
  ,
  width = 8, height = 8,
  file = "plots/mean_exp_all_genes.png"
)

lapply(d.sd[, -1], mean)
lapply(hkg.sub[, -1], mean)

candidate.hkg <- filter(d.sd, mean > 12, within_sd < 0.5, between_sd < 0.5);
dim(candidate.hkg)
qwrite(candidate.hkg$gene, "out/candidate-hkg.vtr");
qwrite(candidate.hkg, "out/candidate-hkg.csv");

candidate.hkg.top <- filter(candidate.hkg, within_sd < 0.35)

hkg.both <- rbind(
  data.frame(hkg.sub, type="known"),
  data.frame(candidate.hkg.top, type="candidate")
);
hkg.both$type <- relevel(factor(hkg.both$type), "known");

qdraw(
  ggplot(d.sd.sub, aes(x=within_sd, y=mean, label=gene))+
    geom_point(alpha=0.1)+
    geom_point(
      data = hkg.both,
      aes(x=within_sd, y=mean, colour=type)
    ) +
    scale_color_npg() +
    geom_label_repel(
      data = hkg.both, aes(colour=type),
      label.padding = 0.1,
      segment.curvature = -1e-20,
      segment.square = TRUE,
      show.legend = FALSE
    ) +
    theme_minimal()
  ,
  width = 8, height = 8,
  file = "plots/mean-vs-within-sd_within.png"
)

# check how many genes are mean > 8
nrow(dplyr::filter(d.sd, mean>=8)) # 9137 genes
nrow(dplyr::filter(d.sd.sub, mean>=8)) # 6779 genes: within_sd < 0.8 AND mean > =8
# 7785 genes: within_sd <1.0 AND mean >=8

# remove genes with very low expression (remove if <8)
##min_hkg_exp <- min(hkg.sub$mean)
d.sd.hexpr <- subset(d.sd.sub, mean>=8)

# plot within_sd agains btwn for highly expressed genes
qdraw(
  ggplot(d.sd.hexpr, aes(x=within_sd, y=between_sd))+
    ggtitle("Between group SD against Within group SD")+
    geom_point(alpha = 0.1, size =0.5)+
    geom_point(data= hkg.sub, color="red")+
    scale_color_npg() +
    geom_label_repel(
      data = hkg.both, aes(colour=type, label=gene),
      label.padding = 0.1,
      segment.curvature = -1e-20,
      segment.square = TRUE,
      show.legend = FALSE
    ) +
    theme(plot.title = element_text(hjust = 0.5))+
    theme_linedraw()+
    geom_abline(intercept=0, slope=sqrt((N - K) / (N - K - 2)),
      color="grey60"),
  width = 8, height = 8,
  file = "plots/sd_high-expr-genes.png"
)  

# ---

# Identify context in which the housekeep genes have lower
# between-group SD

plot_gene <- function(d, expr.mean) {
  ggplot(d, aes(x = sample_id, y = expr)) +
    facet_wrap(~ group, scales="free_x") +
    geom_point(size=0.5) + theme_classic() +
    geom_hline(yintercept = expr.mean, colour = "royalblue") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank(),
      panel.grid.major.y = element_line()
    ) +
    xlab("sample") + ylab("gene expression") +
    ggtitle(gene)
}

make_sd_subtitle <- function(before, after) {
  sprintf("between SD: %.3f -> %.3f", 
    with(before, get_btwn.sd(expr, group)),
    with(after, get_btwn.sd(expr, group))
  )
}

hkg.improve <- filter(hkg.both, between_sd > within_sd);
hkg.improve


gene <- "GAPDH";

expr.d <- data.frame(pheno, expr = mat[gene, ]);
expr.mean <- mean(expr.d$expr);

plot_gene(expr.d, expr.mean)
# expressed more highly in tumour tissues


expr.d.no.nt <- filter(expr.d, sample_type != "NT");
plot_gene(expr.d.no.nt, expr.mean) +
  labs(subtitle = make_sd_subtitle(expr.d, expr.d.no.nt))

expr.d.nt <- filter(expr.d, sample_type == "NT");
plot_gene(expr.d.nt, mean(expr.d.nt$expr)) +
  labs(subtitle = make_sd_subtitle(expr.d, expr.d.nt))

expr.d.nt.no.lung <- filter(expr.d.nt, ! group %in% c("LUAD-NT", "LUSC-NT"));
plot_gene(expr.d.nt.no.lung, mean(expr.d.nt.no.lung$expr)) +
  labs(subtitle = make_sd_subtitle(expr.d, expr.d.nt.no.lung))

# LDHA

gene <- "LDHA";

expr.d <- data.frame(pheno, expr = mat[gene, ]);
expr.mean <- mean(expr.d$expr);

plot_gene(expr.d, expr.mean)
# upregulated in KIRC-TP, HNSC-TP
# downregulated in LGG-TP, THCA-TP

group.exclude <- c("KIRC-TP", "HNSC-TP", "LGG-TP", "THCA-TP");
expr.d.sel <- filter(expr.d, ! group %in% group.exclude);
plot_gene(expr.d.sel, mean(expr.d.sel$expr)) +
  labs(subtitle = make_sd_subtitle(expr.d, expr.d.sel))


gene <- "PTBP1";

expr.d <- data.frame(pheno, expr = mat[gene, ]);
expr.mean <- mean(expr.d$expr);

plot_gene(expr.d, expr.mean)
# upregulated in TGCT-TP, COAD-TP, HNSC-TP, UCEC-TP
# downregulated in LGG-TP, LGG-TR, BRCA-NT, GBM-TP

group.exclude <- c("TGCT-TP", "COAD-TP", "HNSC-TP", "UCEC-TP",
                   "LGG-TP", "LGG-TR", "BRCA-NT", "GBM-TP");
expr.d.sel <- filter(expr.d, ! group %in% group.exclude);
plot_gene(expr.d.sel, mean(expr.d.sel$expr)) +
  labs(subtitle = make_sd_subtitle(expr.d, expr.d.sel))


gene <- "HNRNPC";

expr.d <- data.frame(pheno, expr = mat[gene, ]);
expr.mean <- mean(expr.d$expr);

plot_gene(expr.d, expr.mean)
# upregulated in THYM-TP, OV-TP, TGCT-TP
# downregulated in ESCA-NT, SARC-TP

group.exclude <- c("THYM-TP", "OV-TP", "TGCT-TP", "ESCA-NT", "SARC-TP");
expr.d.sel <- filter(expr.d, ! group %in% group.exclude);
plot_gene(expr.d.sel, mean(expr.d.sel$expr)) +
  labs(subtitle = make_sd_subtitle(expr.d, expr.d.sel))


gene <- "UBC";

expr.d <- data.frame(pheno, expr = mat[gene, ]);
expr.mean <- mean(expr.d$expr);

plot_gene(expr.d, expr.mean)
# upregulated in KIRC-TP, HNSC-TP
# downregulated in UVM-TP, LAML-TB, READ-TP

group.exclude <- c("HNSC-TP", "KIRC-TP", "UVM-TP", "LAML-TB","READ-TP");
expr.d.sel <- filter(expr.d, ! group %in% group.exclude);
plot_gene(expr.d.sel, mean(expr.d.sel$expr)) +
  labs(subtitle = make_sd_subtitle(expr.d, expr.d.sel))

gene <- "ACTB";

expr.d <- data.frame(pheno, expr = mat[gene, ]);
expr.mean <- mean(expr.d$expr);

plot_gene(expr.d, expr.mean)
# upregulated in DLBC-TP, LUAD-NT, LUSC-NT
# downregulated in PCPG-TP, PRAD-TP, KIRC-NT, LIHC-NT

group.exclude <- c("DLBC-TP", "LUAD-NT", "LUSC-NT","PRAD-TP", "PCPG-TP", "KIRC-NT", "LIHC-NT");
expr.d.sel <- filter(expr.d, ! group %in% group.exclude);
plot_gene(expr.d.sel, mean(expr.d.sel$expr)) +
  labs(subtitle = make_sd_subtitle(expr.d, expr.d.sel))

