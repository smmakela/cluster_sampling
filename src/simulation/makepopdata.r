# make pop data
  args <- commandArgs(FALSE)
  use.sizes <- as.numeric(args[length(args)])
  #seed <- 1
  numclusters <- 300
  clustersize.range <- c(100, 1000)
  unitcovar.range <- c(20, 45)
  rootdir <- "/vega/stats/users/smm2253/Projects/Cluster_Sampling/"
  source("/vega/stats/users/smm2253/Projects/Cluster_Sampling/Code/Simplify/vary_K/makedata.r")
  makedata(seed = NA, numclusters, clustersize.range, unitcovar.range, rootdir, use.sizes)
