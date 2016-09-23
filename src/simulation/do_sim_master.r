# get the current task id -- we'll use it to generate new populations
args <- commandArgs(FALSE)
print(args)
len <- length(args)
sim <- -1*as.numeric(args[len])
print(sim)
seed <- NULL
numclusters <- 300
clustersize.range <- c(100, 1000)
unitcovar.range <- c(20, 45)
source("/vega/stats/users/smm2253/Projects/Cluster_Sampling/Code/Simplify/vary_K/sim_master.r")
sim_master(sim, seed, numclusters, clustersize.range, unitcovar.range)

