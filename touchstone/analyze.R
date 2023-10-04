# find Housekeeping Genes which have low variance among the samples

# calculate mean values
mean_df <- data.frame(matrix(ncol = length(filter_sample_group$group), nrow = nrow(touchstone_df)))
rownames(mean_df) <- rownames(touchstone_df)
colnames(mean_df) <- filter_sample_group$group
# for (group in filter_sample_group$group){
#   df_groups <- filter_sample_group[filter_sample_group$group == group, "inst_id"]
#   mean_values <- rowMeans(touchstone_df[,as.character(df_groups)])
#   df_mean[,as.character(group)] <- mean_values
# }

df_VCAP_UnTrt_24 <- touchstone_df[,filter_sample_group[filter_sample_group$group == "VCAP_UnTrt_24", "inst_id"]]
mean_values <- rowMeans(df_VCAP_UnTrt_24)
mean_values <- as.data.frame(mean_values)