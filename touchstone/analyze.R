# Load required libraries for the analysis.
library(io)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggsci)

# find Housekeeping Genes which have low variance among the samples

# import functions to calculate SD (at the same time, calculate mean values)
source("R/sd_function.R")

load("touchstone_df.RData")
