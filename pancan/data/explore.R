library(io)

proj.info <- qread("../annot/sample-project-info_pancan.tsv")

type.info <- qread("../annot/sample-type-code_tcga.tsv")


x <- qread("expr_pancan.rds");
mat <- x$data;

# simplify colnames in mat

#idx <- match("TCGA-OR-A5J1-01A", proj.info$sample);
#idx <- match("TCGA-OR-A5J2-01A", proj.info$sample)
#idx <- match("TCGA-OR-A5J3-01A", proj.info$sample)

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

table(pheno$group)
sum(is.na(pheno$group))


# calculate standard deviation within the group
BiocManager::install("edgeR")
library("limma")
library("edgeR")

##pheno_ACCTP <- na.omit(pheno[pheno$group == "ACC-TP",])

gene <- "UBC"

within_means = tapply(mat[gene,], pheno$group, mean)
within_sds = tapply(mat[gene,], pheno$group, sd)

##d_within <- data.frame(
##sample_id = pheno_ACCTP$sample_id,

##)

# calculate standard deviation betw3en the groups
overall_mean = mean(mat[gene,])
between_sd = sqrt(sum((within_means-overall_mean)^2)/(nrow(within_means)-1))
