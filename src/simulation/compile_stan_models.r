# Purpose: quick script to compile stan models

# Main project directory
rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"

# Loop through stan models, compile, and save
stanmod_list <- c("cluster_inds_only", "knowsizes", "bb", "negbin")
for (j in 1:length(stanmod_list)) {
  stanmod <- stan_model(file = paste0(rootdir, "/src/", stanmod_name, ".stan")
  save(stanmod, file = paste0(rootdir, "/src/", stanmod_name, ".RData")) 
}

