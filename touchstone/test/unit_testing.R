library(io)

source("R/sd_function.R")

d <- qread("test/test_data.csv");

# TODO add manually calculated answers here
ans <- c(1.324645702619340, 3.173520835823380);

values <- d$expression;
groups <- d$group;

within_sd <- get_within.sd(values, groups);
between_sd <- get_btwn.sd(values, groups);

print(paste("within-group sd: ", within_sd))
# 1.32464568298884
print(paste("between-group sd: ", between_sd))
# 3.17352085360015

# FIXME: need ground-truth answer here!
stopifnot(abs(within_sd - ans[1]) < 1e-7) # Passed!
stopifnot(abs(between_sd - ans[2]) < 1e-7) # Passed!
