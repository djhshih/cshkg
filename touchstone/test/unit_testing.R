library(io)

source("../R/sd_function.R")

d <- qread("test_data.csv");

# TODO add manually calculated answers here
ans <- c(NA, NA);

values <- d$expression;
groups <- d$group;

tapply(values, groups, mean)

within_sd <- get_within.sd(values, groups);
between_sd <- get_btwn.sd(values, groups);

print(paste("within-group sd: ", within_sd))
print(paste("between-group sd: ", between_sd))

# FIXME: need ground-truth answer here!
stopifnot(abs(within_sd - ans[1]) < 1e-22);
stopifnot(abs(within_sd - ans[2]) < 1e-22);

