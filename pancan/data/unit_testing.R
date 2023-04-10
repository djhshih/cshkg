library(io)
library(ggplot2)

testdf <- qread("unit testing.csv");

testidx <- c(1,1,1,2,2,2,3,3,3,4,4,4)

testmat <- testdf[,"expression"]

testwithin_means = tapply(testmat, testidx, mean, na.rm=TRUE)
testwithin_sd <- sqrt(
    sum((testmat - testwithin_means[testidx])^2, na.rm=TRUE) / 
      (length(testmat) - 1))
  
  
# calculate standard deviation between the groups
testoverall_mean = mean(testwithin_means)
  
testbetween_sd = sqrt(sum((testwithin_means-testoverall_mean)^2)/(length(testwithin_means)-1))

print(paste("within-group sd: ", testwithin_sd))
print(paste("between-group sd: ", testbetween_sd))
