library(io)

subset_m <- qread("./subset_m.rds");

dist_m <- dist(subset_m); # not working

sub_subset_m <- subset_m[1:10, 1:10];
dist_m <- dist(sub_subset_m); 
hc <- hclust(dist_m) # working