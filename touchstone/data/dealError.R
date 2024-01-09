# Perform hierarchical clustering on the transposed dataset (samples in rows, genes in columns)
m_rev <- na.omit(t(m))
#d <- dist(m_rev)
#hc <- hclust(d) # not working
# Subset matrix
subset_size <- 10000 # worked until 10000 samples
num_rows <- nrow(m_rev)
# Ensure that the subset size is not larger than the number of columns
if (subset_size > num_rows) {
  stop("Subset size is larger than the number of columns in the dataset.")
}
# Error: Error: vector memory exhausted (limit reached?)
# Sys.setenv("R_MAX_VSIZE" = 32000000000) # not working

subset_m <- m_rev[sample(num_rows, size = subset_size, replace = FALSE),]
saveRDS(subset_m, file = "subset_m.rds")
hc <- hclust(dist(as.matrix(subset_m)))

library(io)

subset_m <- qread("./subset_m.rds");

dist_m <- dist(subset_m); # not working

sub_subset_m <- subset_m[1:10, 1:10];
dist_m <- dist(sub_subset_m);
hc <- hclust(dist_m) # working