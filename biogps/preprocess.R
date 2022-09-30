library(io)
library(annotate)
library(hgu133a.db)

db <- hgu133a.db;

# load expression data
x <- qread("U133AGNF1B.gcrma.csv", type="csv");
rownames(x) <- x[,1];
x <- x[,-1];
x <- as.matrix(x);

# map probe to genes
mapping <- select(db, rownames(x), "SYMBOL");
head(mapping)
idx <- match(rownames(x), mapping$PROBEID);
genes <- mapping$SYMBOL[idx];

sum(is.na(genes))

dim(x)

# remove NA genes
valid <- !is.na(genes)
x <- x[valid, ];
genes <- genes[valid];

rownames(x) <- genes;

sum(duplicated(rownames(x)))

# remove duplicated genes
idx <- duplicated(rownames(x));
x <- x[!idx, ];

sum(duplicated(rownames(x)))
dim(x)

qwrite(x, "gnf1h.rds");

