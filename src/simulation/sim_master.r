# Author: Susanna Makela
# Date: 13 Nov 2014
# Purpose: master file for binomial model

sim_master <- function(sim, numclusters, outcome.type) {
# Arguments:
#  sim -- simulation number; used to save and later reload sample data
#  numclusters -- number of clusters in population (needed in svy_ests.r)
#  outcome.type -- whether outcome being considered is continuous or binary

  #############################################################################
  ### Setup of directories
  #############################################################################
    libdir <- "/vega/stats/users/smm2253/rpackages"
    .libPaths(libdir)
    rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
    Sys.setenv(HOME = rootdir)

  #############################################################################
  ### Source files
  #############################################################################
    source(paste(rootdir, "/src/simulation/sampledata.r", sep = ""))
    source(paste(rootdir, "/src/simulation/runstan.r", sep = ""))
    source(paste(rootdir, "/src/simulation/lmer_compare.r", sep = ""))
    source(paste(rootdir, "/src/simulation/svy_ests.r", sep = ""))

  #############################################################################
  ### Load libraries, set up parallel stuff
  #############################################################################
    library(lme4)
    library(rstan)
    library(foreach)
    library(doParallel) 
    library(iterators)
    cl <- makeCluster(min(detectCores(), 10),
                      outfile = paste(rootdir,
                                      "/src/simulation/parallel_output.txt",
                                      sep = ""))
    registerDoParallel(cl)
    clusterEvalQ(cl, .libPaths( "/vega/stats/users/smm2253/rpackages"))
    #clusterEvalQ(cl, library(iterators, lib.loc = "/vega/stats/users/smm2253/rpackages"))
    #clusterEvalQ(cl, library(foreach, lib.loc = "/vega/stats/users/smm2253/rpackages"))
    #clusterEvalQ(cl, library(doParallel, lib.loc = "/vega/stats/users/smm2253/rpackages"))


  #############################################################################
  ### Create list of parameters to loop through for sim
  #############################################################################
    num.clusters.list <- c(15, 30, 60, 100)
    num.units.list <- c(0.05, 0.1, 0.25, 0.5, 1, 10, 50, 100)
    use.sizes.list <- c(0, 1)
    outcome.type.list <- c("continuous", "binary")
    tmplist <- list(num.clusters.list = num.clusters.list,
                    num.units.list = num.units.list,
                    use.sizes.list = use.sizes.list,
                    outcome.type.list = outcome.type.list)
    sim.params <- expand.grid(tmplist) 

  #############################################################################
  ### Loop through cluster/unit lists
  #############################################################################
    stanmod_list <- c("cluster_inds_only", "knowsizes", "bb", "negbin")
    loopres <- foreach (k = 1:nrow(sim.params),
                        .packages = c("rstan", "dplyr"),
                        .export = ls(envir = globalenv()),
                        .verbose=TRUE) %dopar% {

      # Set parameters for this simulation
      num.clusters <- sim.params[k, "num.cluster.list"]
      num.units <- sim.params[k, "num.units.list"]
      use.sizes <- sim.params[k, "use.sizes.list"]
      outcome.type <- sim.params[k, "outcome.type.list"]

      # Sample data using above parameters
      print("Sampling data")
      print(Sys.time())
      sampledata(num.clusters, num.units, use.sizes, outcome.type, rootdir, sim)

      # Run stan models
      for (p in 1:length(stanmod_list)) {
        stanmod_name <- stanmod_list[p]
        stanmod <- load(paste0(rootdir, "/src/", stanmod_name, ".RData"))
        cat("##################################################################################\n")
        cat("##################################################################################\n")
        print(paste0("Starting to run use.sizes = ", use.sizes,
                     ", stanmod = ", stanmod_name,
                     ", num.clusters = ", num.clusters,
                     ", num.units = ", num.units))
        print(Sys.time())
        stanout <- runstan(num.clusters, num.units, use.sizes,
                           rootdir, sim, stanmod, stanmod_name)
        cat("##################################################################################\n")
      } # end stanmod loop

      # Compare using lmer
      print("Running lmer_compare")
      print(Sys.time())
      lmer_compare(num.clusters, num.units, use.sizes, outcome.type,
                   rootdir, sim)

      # Estimate ybar using survey package
      print("Running svy_ests")
      print(Sys.time())
      J <- numclusters # number of clusters in the population
      svy_ests(J, num.clusters, num.units, use.size, outcome.type, rootdir, sim)
      print("##################################################################################")
      print("##################################################################################")

      # Delete the sampled data to save space
        if (num.units <= 1) {
          nunits <- paste(100*num.units, "pct", sep = "")
        } else {
          nunits <- num.units
        }
          fil1 <- paste(rootdir, "/Data/Simplify/vary_K/sampledata_usesizes_", use.sizes, "_nclusters_", num.clusters,
                        "_nunits_", nunits, "_sim_", sim, ".RData", sep = "")
          if (file.exists(fil1)) {
            file.remove(fil1)
          }
        } # end unit loop
      } # end cluster loop
    } # end use sizes loop

    #stopCluster(cl) # close workers

    # remove the pop data file with the indices
      file.remove(paste(rootdir, "/Data/Simplify/vary_K/popdata_usesizes_", use.sizes,
				 "_inds_nclusters_", num.clusters, "_nunits_", nunits,
				 "_sim_", sim, ".RData", sep = ""))
} # end function
