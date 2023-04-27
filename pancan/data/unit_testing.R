library(io)

source("sd_function.R")

testdf <- qread("unit_testing.csv");

testidx <- c(1,1,1,2,2,2,3,3,3,4,4,4)

testmat <- testdf[,"expression"]

testwithin_sd <- get_within.sd(testmat, testidx);
  
# calculate standard deviation between the groups
testbetween_sd <- get_btwn.sd(testmat, testidx);

print(paste("within-group sd: ", testwithin_sd))
print(paste("between-group sd: ", testbetween_sd))

# FIXME: need ground-truth answer here!

