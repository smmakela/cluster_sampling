#!/usr/bin/R Rscript

'
Usage: sim_master.R  --simno=<simnumber> [options]

Options:
  -s --simno <simnumber>  Current simulation number
  -n --numclusters <J>    Number of clusters in population [default: 100]

' -> doc

  #############################################################################
  ### Set lib paths, source files 
  #############################################################################
    libdir <- "/vega/stats/users/smm2253/rpackages"
    .libPaths(libdir)
    rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
    Sys.setenv(HOME = rootdir)
    
    source(paste0(rootdir, "src/simulation/sampledata.r"))
    source(paste0(rootdir, "src/simulation/runstan.r"))
    source(paste0(rootdir, "src/analysis/lmer_compare.R"))
    source(paste0(rootdir, "src/analysis/svy_ests.R"))

  #############################################################################
  ### Load libraries, get options, set up parallel stuff
  #############################################################################
    require(docopt)
    require(methods)
    require(plyr)
    require(dplyr)
    require(lme4)
    require(rstan)
    require(foreach)
    require(doParallel) 
    require(iterators)

    # Store the docopt options as variables we can use in the code
    opts <- docopt(doc) 
    opts.names <- names(opts)
    opts.names <- opts.names[-grep("--", opts.names)]
    print(str(opts))
    print(opts.names)

    # The options are read in as strings, so make them numeric here
    for (j in 1:length(opts.names)) {
      opts[[j]] <- as.numeric(opts[[j]])
    }

    # Set up parallel parameters
    cl <- makeCluster(min(detectCores(), 10),
                      outfile = paste(rootdir,
                                      "/src/simulation/parallel_output.txt",
                                      sep = ""))
    registerDoParallel(cl)
    clusterEvalQ(cl, .libPaths( "/vega/stats/users/smm2253/rpackages"))
    print(getDoParWorkers()) # make sure foreach will actually run in parallel
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
                    num.units.list    = num.units.list,
                    use.sizes.list    = use.sizes.list,
                    outcome.type.list = outcome.type.list)
    sim.params <- expand.grid(tmplist)
    sim.params$outcome.type.list <- as.character(sim.params$outcome.type.list) 

  #############################################################################
  ### Loop through cluster/unit lists
  #############################################################################
    stanmod_list <- c("cluster_inds_only", "knowsizes", "bb", "negbin")
#    loopres <- foreach (k = 1:nrow(sim.params),
#                        .packages = c("rstan", "dplyr"),
#                        .export = ls(envir = globalenv()),
#                        .verbose=TRUE) %dopar% {

for (k in 1:nrow(sim.params)) {
      # Set parameters for this simulation
      num.clusters <- sim.params[k, "num.clusters.list"]
      num.units    <- sim.params[k, "num.units.list"]
      use.sizes    <- sim.params[k, "use.sizes.list"]
      outcome.type <- sim.params[k, "outcome.type.list"]

      # Print a message about which parameters we're running now
      cat("Running in parallel for", num.clusters, "clusters,", num.units,
          "units, use_sizes =", use.sizes, ", and", outcome.type, "outcomes.\n")
  
      # Sample data using above parameters
      print("Sampling data")
      print(Sys.time())
cat("current sim number is", simno)
      sampledata(num.clusters, num.units, use.sizes, outcome.type, rootdir, simno)
  
      # Run stan models
#      for (p in 1:length(stanmod_list)) {
#        stanmod_name <- stanmod_list[p]
#        stanmod <- load(paste0(rootdir, "/src/", stanmod_name, ".RData"))
#        cat("##################################################################################\n")
#        cat("##################################################################################\n")
#        print(paste0("Starting to run use.sizes = ", use.sizes,
#                     ", stanmod = ", stanmod_name,
#                     ", num.clusters = ", num.clusters,
#                     ", num.units = ", num.units))
#        print(Sys.time())
#        stanout <- runstan(num.clusters, num.units, use.sizes,
#                           rootdir, sim, stanmod, stanmod_name)
#        cat("##################################################################################\n")
#      } # end stanmod loop
#  
#      # Compare using lmer
#      print("Running lmer_compare")
#      print(Sys.time())
#      lmer_compare(num.clusters, num.units, use.sizes, outcome.type, rootdir, sim)
#  
#      # Estimate ybar using survey package
#      print("Running svy_ests")
#      print(Sys.time())
#      J <- numclusters # number of clusters in the population
#      svy_ests(J, num.clusters, num.units, use.size, outcome.type, rootdir, sim)
#      print("##################################################################################")
#      print("##################################################################################")
#  
#      # Delete the sampled data to save space
#      if (num.units <= 1) {
#        nunits <- paste(100*num.units, "pct", sep = "")
#      } else {
#        nunits <- num.units
#      }
#      fil1 <- paste(rootdir, "/Data/Simplify/vary_K/sampledata_usesizes_", use.sizes, "_nclusters_", num.clusters,
#                    "_nunits_", nunits, "_sim_", sim, ".RData", sep = "")
#      if (file.exists(fil1)) {
#        file.remove(fil1)
#      }
 } # end parallel loop

  stopCluster(cl) # close workers


