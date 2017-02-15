#!/usr/bin/R Rscript
# Purpose: quick script to compile stan models

# Main project directory
rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
.libPaths("/vega/stats/users/smm2253/rpackages")
library(rstan)

# Loop through stan models, compile, and save
stanmod_list <- c("cluster_inds_only", "knowsizes", "bb", "negbin")
for (j in 1:length(stanmod_list)) {
  stanmod_name <- stanmod_list[j]
  stanmod <- stan_model(file = paste0(rootdir, "/src/analysis/", stanmod_name, ".stan"))
  save(stanmod, file = paste0(rootdir, "/src/analysis/", stanmod_name, ".RData")) 
}

