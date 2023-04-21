library(io)
library(ggplot2)

source("sd_function.R")

proj.info <- qread("../annot/sample-project-info_pancan.tsv")

type.info <- qread("../annot/sample-type-code_tcga.tsv")


x <- qread("expr_pancan.rds")
mat0 <- x$data;

# NB beware of copy-on-write
mat <- mat0;

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

table(pheno$group)
sum(is.na(pheno$group))

stopifnot(pheno$sample_id == colnames(mat))


#Calculate standard deviation within and between the groups

geneofinterest <- sample(rownames(mat), 1000)
  #sample.int(nrow(mat),10)
sd_plot <- data.frame()

for (gene in geneofinterest) {
  v <- mat[gene,]
  
  within_sd = get_within.sd(v, pheno$group)
  between_sd = get_btwn.sd(v, pheno$group)
  
  # plot the graph between within the group sd and between the group sd for different genes
  sd_plot <- rbind(sd_plot, data.frame(within_sd, between_sd, row.names=gene))
  
}

ggplot(sd_plot, aes(x=within_sd, y=between_sd))+
  geom_point()+
  geom_text(
    label = rownames(sd_plot),
    hjust = 0, nudge_x = 0.005
  )
