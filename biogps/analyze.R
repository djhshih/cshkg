library(io)
library(ggplot2)
library(dplyr)

x <- qread("gnf1h.rds");
pheno <- qread("gnf1h_pheno.rds");

# examine retina tissue
idx <- pheno$tissue == "retina";
head(x[, idx])

#gene <- "GAPDH";
#gene <- "ACTB";
gene <- "UBC";

means <- tapply(x[gene, ], pheno$group, mean);
sds <- tapply(x[gene, ], pheno$group, sd);

d.grouped <- data.frame(
	group = names(means),
	mean = means,
	lower = means - sds,
	upper = means + sds
);

ggplot(d.grouped, aes(x=group, y=mean)) +
	theme_bw() +
	geom_col() +
	geom_errorbar(aes(ymin=lower, ymax=upper), width=0.3) +
	theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
	xlab("") + ylab("expression")

d.single <- data.frame(
	sample_id = colnames(x),
	expr = x[gene, ],
	row.names = NULL
);
d.single <- left_join(d.single, pheno, by="sample_id")

qdraw(
	ggplot(d.single, aes(x=sample_id, y=expr)) +
		theme_classic() +
		geom_point() +
		facet_wrap(~ group, nrow=1, scales="free_x") +
		theme(
			axis.text.x = element_blank(),
			axis.ticks.x = element_blank(),
			panel.grid.major.y = element_line(),
			strip.background = element_blank()
		) +
		xlab("sample") + ylab("expression")
	,
	width = 15
)

