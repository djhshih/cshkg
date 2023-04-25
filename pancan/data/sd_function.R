#'v' is a table containing the tissue group and expression of each sample
#'group' is index of the group

#Calculate within-group standard deviation
get_within.sd <- function(v, group){
  group_mean <- tapply(v, group, mean, na.rm = TRUE)
  within_sd = sqrt(sum((v-group_mean[group])^2, na.rm = TRUE) / (length(v)-1))
  return(within_sd)
}

#Calculate between-group standard deviation
get_btwn.sd <- function(v, group){
  group_mean <- tapply(v, group, mean, na.rm = TRUE)
  overall_mean <- mean(group_mean, na.rm=TRUE)
  #print(sum((group_mean-overall_mean)^2))
  between_sd <- sqrt(sum((group_mean-overall_mean)^2, na.rm=TRUE) / (length(group_mean)-1))
  return(between_sd)
}

