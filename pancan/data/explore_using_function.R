
library(io)
library(ggplot2)
library(ggrepel)


source("sd_function.R")

proj.info <- qread("../annot/sample-project-info_pancan.tsv")

type.info <- qread("../annot/sample-type-code_tcga.tsv")


x <- qread("expr_pancan.rds")
mat <- x$data;

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

rare.groups <- names(which(table(pheno$group) < 10));

levels(pheno$group)[levels(pheno$group) %in% rare.groups] <- NA;

idx.remove <- is.na(pheno$group);
table(idx.remove)


table(pheno$group)
sum(is.na(pheno$group))

stopifnot(pheno$sample_id == colnames(mat))

mat <- mat[, !idx.remove];
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

qwrite(d.sd, "sds.rds");

# SSX9P has many NA and some groups have only NA
# for groups with all missing value, the mean is NaN
v.problematic <- mat["SSX9P", ];
table(is.na(v.problematic))
gmean.problematic <- tapply(v.problematic, pheno$group, mean, na.rm=TRUE);
gmean.problematic
overall_mean <- mean(gmean.problematic, na.rm=TRUE)
#print(sum((group_mean-overall_mean)^2))
between_sd.problematic <- sqrt(sum((gmean.problematic-overall_mean)^2, na.rm=TRUE) / (length(gmean.problematic)-1))

housekeeping <- c("ACTB", "UBC", "GAPDH", "TBP", "RPS18", "G6PD", "HPRT1", "LDHA", "RPL19","RPL18","RPL11","RPL32", "PGK1", "PPIA", "RPS18", "ASNS","ATP2B4", "PEX19","RXRA");  # TODO: Add more!
non_hk <- c("TP53","EGFR","AKT1")
dplyr::filter(d.sd, gene %in% housekeeping)
dplyr::filter(d.sd, gene %in% non_hk)

d.sd.sub <- subset(d.sd, within_sd < 0.8)
hkg.sub <- subset(d.sd.sub, d.sd.sub$gene %in% housekeeping)
non_hk.sub <- subset(d.sd, d.sd$gene %in% non_hk)

dim(d.sd.sub)

ggplot(d.sd.sub, aes(x=within_sd, y=between_sd, label=gene))+
  #geom_text_repel() +
  geom_point(alpha=0.1) + coord_fixed() +
  geom_point(hkg.sub, mapping=aes(color="red"))+ geom_label_repel(hkg.sub, mapping=aes(color="red"),label.padding = 0.1)+
  geom_abline(yintercept=0, slope = sqrt((N - K) / (N - K - 2)))


ggplot(subset(d.sd.sub, within_sd < 0.25), aes(x=within_sd, y=between_sd, label=gene))+
  geom_label_repel() +
  geom_point() + coord_fixed() +
  geom_point(hkg.sub, mapping=aes(color="red"))+ geom_label_repel(hkg.sub, mapping=aes(color="red"),label.padding = 0.1)+
  geom_point(non_hk.sub, mapping=aes(color="blue"))+ geom_label(non_hk.sub, mapping=aes(color="blue"))+
  geom_abline(yintercept=0, slope = sqrt((N - K) / (N - K - 2)))


ggplot(subset(d.sd.sub, within_sd < 0.25), aes(x=within_sd, y=mean, label=gene))+
  geom_label_repel() +
  geom_point() +
  geom_point(hkg.sub, mapping=aes(color="red"))+ geom_label_repel(hkg.sub, mapping=aes(color="red"),label.padding = 0.1)+
  geom_point(non_hk.sub, mapping=aes(color="blue"))+ geom_label(non_hk.sub, mapping=aes(color="blue"))
