# 'v' is a vector containing expression of a gene across samples
# 'group' is group membership index

# Calculate within-group standard deviation
get_within.sd <- function(v, group){
  group_mean <- tapply(v, group, mean, na.rm = TRUE)
  sqrt(
    sum( (v - group_mean[group])^2, na.rm = TRUE ) / 
    (sum(!is.na(v)) - 1)
  )
}

# Calculate between-group standard deviation
get_btwn.sd <- function(v, group){
  group_mean <- tapply(v, group, mean, na.rm = TRUE)
  overall_mean <- mean(group_mean, na.rm=TRUE)
  sqrt(
    sum( (group_mean - overall_mean)^2, na.rm=TRUE ) / 
    (sum(!is.na(group_mean)) - 1)
  )
}

