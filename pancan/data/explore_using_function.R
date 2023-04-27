library(io)
library(ggplot2)
library(ggrepel)

#import functions to calculate standard deviation
source("sd_function.R")

proj.info <- qread("../annot/sample-project-info_pancan.tsv")
type.info <- qread("../annot/sample-type-code_tcga.tsv")

x <- qread("expr_pancan.rds")
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
  geom_bar(stat = "identity") + theme_minimal()+ 
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

qwrite(d.sd, "sds.rds");
#d.sd <- qread("sds.rds") # load d.sd 

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
                  "RPS18","ASNS","ATP2B4", "PEX19","RXRA",
                  "RPL13A","PPIA","SDHA","PPIA")  # TODO: Add more!
common_gene <- c("TP53","EGFR","AKT1")
hkg.sub <- dplyr::filter(d.sd, gene %in% housekeeping)
dplyr::filter(d.sd, gene %in% common_gene)


# Set ggplot options globally
#options(ggrepel.max.overlaps = Inf)

# figure2: sd of all genes (btwn against within)
qdraw(
  ggplot(d.sd, aes(x=within_sd, y=between_sd))+
    ggtitle("Between group SD against Within group SD")+
    geom_point(alpha = 0.9, size =0.5)+
    geom_point(data= hkg.sub, color="red")+
    geom_label_repel(data = hkg.sub, aes(label=gene), color="red", label.padding = 0.1, max.overlaps = Inf, segment.curvature=-1e-20, segment.square = TRUE)+
    theme(plot.title = element_text(hjust = 0.5))+
    theme_linedraw()+
    geom_abline(intercept=0, slope=sqrt((N - K) / (N - K - 2))),
  width = 8, height = 8,
  file = "sd_all_genes_with_hkg_highlights2.png"
)  


# subset all genes whose within_sd < 0.8
d.sd.sub <- subset(d.sd, within_sd < 0.8)



#remove genes with very low expression
min_hkg_exp <- min(hkg.sub$mean)
d.sd.sub <- subset(d.sd.sub, mean>=min_hkg_exp)

#subset hkg and common genes
hkg.sub <- subset(d.sd.sub, d.sd.sub$gene %in% housekeeping)
common.sub <- subset(d.sd, d.sd$gene %in% common_gene)

dim(d.sd.sub)

library(io)
#within_sd <0.8
qdraw(
  ggplot(d.sd.sub, aes(x=within_sd, y=between_sd, label=gene))+
    #geom_text_repel() +
    geom_point(alpha=0.1) + coord_fixed() +
    geom_point(data= hkg.sub, aes(color=gene))+
    geom_label_repel(data = hkg.sub, aes(color=gene),label.padding = 0.1)+
    geom_abline(intercept=0, slope = sqrt((N - K) / (N - K - 2)))+
    theme_minimal()
  ,
  width = 8, height = 8,
  file = "within_sd_less_than_0.8.png"
)

#within_sd < 0.25
ggplot(subset(d.sd.sub, within_sd < 0.25), aes(x=within_sd, y=between_sd, label=gene))+
  geom_label_repel() +
  geom_point() + coord_fixed() +
  geom_point(hkg.sub, mapping=aes(color="red"))+ geom_label_repel(hkg.sub, mapping=aes(color="red"),label.padding = 0.1)+
  geom_point(common.sub, mapping=aes(color="blue"))+ geom_label(common.sub, mapping=aes(color="blue"))+
  geom_abline(yintercept=0, slope = sqrt((N - K) / (N - K - 2)))


#mean expression of genes
ggplot(subset(d.sd.sub, within_sd < 0.25), aes(x=within_sd, y=mean, label=gene))+
  geom_label_repel() +
  geom_point() +
  geom_point(hkg.sub, mapping=aes(color="red"))+ geom_label_repel(hkg.sub, mapping=aes(color="red"),label.padding = 0.1)+
  geom_point(common.sub, mapping=aes(color="blue"))+ geom_label(common.sub, mapping=aes(color="blue"))


