#'v' is a table containing the tissue group and expression of each sample
#'group' is index of the group

#Calculate within-group standard deviation
get_within.sd <- function(v, group){
  within_mean <- tapply(v, group, mean, na.rm = TRUE)
  within_sd = sqrt(sum((v-within_mean[group])^2, na.rm = TRUE) / (length(v)-1))
  return(within_sd)
}

#Calculate between-group standard deviation
get_btwn.sd <- function(v, group){
  within_mean <- tapply(v, group, mean, na.rm = TRUE)
  overall_mean = mean(within_mean)
  between_sd = sqrt(sum((within_mean-overall_mean)^2) / (length(within_mean)-1))
  return(between_sd)
}

