library(io)
library(ggplot2)

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

geneofinterest <- c("UBC", "GAPDH")
dfplot <- data.frame()

for (gene in geneofinterest) {
  
  v <- mat[gene, ];
  
  within_means = tapply(v, pheno$group, mean)
  within_sds = tapply(v, pheno$group, sd)
  
  within_sd <- sqrt(
    sum((v - within_means[as.numeric(pheno$group)])^2, na.rm=TRUE) / 
      (length(v) - 1))
  
  ##within_sd_ESCANT = tapply(mat[gene,], pheno$group =="ESCA-NT", sd)
  
  # calculate standard deviation between the groups
  
  overall_mean = mean(within_means)

  between_sd = sqrt(sum((within_means-overall_mean)^2)/(length(within_means)-1))
  
  
  # plot the graph between within the group sd and between the group sd for different genes
  ##dfplot <- data.frame(within_sd, between_sd, row.names = gene)
  dfplot <- rbind(dfplot, data.frame(within_sd, between_sd, row.names = gene))
}

ggplot(dfplot, aes(x=within_sd, y=between_sd, label = gene))+
  geom_point()+
  geom_text(hjust = 0, nudge_x = 0.005)
