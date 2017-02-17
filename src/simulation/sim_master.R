#!/usr/bin/R Rscript

'
Usage: sim_master.R  --simno=<simnumber> --use_sizes=<us_val> --outcome_type=<ot_val> [options]

Options:
  -s --simno <simnumber>      Current simulation number
  -u --use_sizes <us_val>     Whether outcomes depend on cluster sizes
  -o --outcome_type <ot_val>  Whether outcomes are continuous or binary
  -n --numclusters <J>        Number of clusters in population [default: 100]

' -> doc

  #############################################################################
  ### Set lib paths, source files 
  #############################################################################
    libdir <- "/vega/stats/users/smm2253/rpackages"
    .libPaths(libdir)
    rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
    Sys.setenv(HOME = rootdir)
    
    source(paste0(rootdir, "src/simulation/sampledata.R"))
    source(paste0(rootdir, "src/simulation/runstan.R"))
    source(paste0(rootdir, "src/simulation/lmer_compare.R"))
    source(paste0(rootdir, "src/simulation/svy_ests.R"))

  #############################################################################
  ### Load libraries, get options, set up parallel stuff
  #############################################################################
    require(docopt)
    require(methods)
    require(plyr)
    require(dplyr)
    require(tidyr)
    require(lme4)
    require(rstan)
    require(foreach)
    require(doParallel) 
    require(iterators)
    require(survey)

    # Store the docopt options as variables we can use in the code
    opts <- docopt(doc) 
    opts.names <- names(opts)

    # The options are listed twice, once with "--" in front of the option name
    # and once without, so remove the "--" ones
    opts.names <- opts.names[-grep("--", opts.names)]
    opts <- opts[opts.names]
    print(str(opts))
    print(opts.names)

    # The options are read in as strings, so make them numeric here
    for (j in 1:length(opts.names)) {
      if (is.character(opts[[j]])) {
        print(paste0("Assigning ", opts[[j]], " to ", opts.names[j]))
        assign(opts.names[j], opts[[j]])
      } else {
        print(paste0("Assigning ", as.numeric(opts[[j]]), " to ", opts.names[j]))
        assign(opts.names[j], as.numeric(opts[[j]]))
      }
    }

    # Parameters for stan
    num.iter <- 1000
    num.chains <- 4

    # Set up parallel parameters
#    cl <- makeCluster(min(detectCores(), 10),
#                      type = "FORK",
#                      outfile = paste0(rootdir,
#                                       "/output/simulation/parallel_output_",
#                                       "usesizes_", use_sizes, "_",
#                                       outcome_type, "_simno_", simno, ".txt"))
#    registerDoParallel(cl)
#    clusterEvalQ(cl, .libPaths( "/vega/stats/users/smm2253/rpackages"))
#    print(getDoParWorkers()) # make sure foreach will actually run in parallel
#    clusterEvalQ(cl, library(iterators, lib.loc = "/vega/stats/users/smm2253/rpackages"))
#    clusterEvalQ(cl, library(foreach, lib.loc = "/vega/stats/users/smm2253/rpackages"))
#    clusterEvalQ(cl, library(doParallel, lib.loc = "/vega/stats/users/smm2253/rpackages"))

  #############################################################################
  ### Create list of parameters to loop through for sim
  #############################################################################
    num.clusters.list <- c(5)
    #num.clusters.list <- c(5, 10, 20, 50)
    num.units.list <- c(0.05)
    #num.units.list <- c(0.05, 0.1, 0.25, 0.5, 1, 10, 50, 100)
    #use.sizes.list <- c(0, 1)
    #outcome.type.list <- c("continuous", "binary")
    tmplist <- list(num.clusters.list = num.clusters.list,
                    num.units.list    = num.units.list)
                    #use.sizes.list    = use.sizes.list,
                    #outcome.type.list = outcome.type.list)
    sim.params <- expand.grid(tmplist)
    #sim.params$outcome.type.list <- as.character(sim.params$outcome.type.list) 
#print("sim params:")
#print(sim.params)
  #############################################################################
  ### Loop through cluster/unit lists
  #############################################################################
    #stanmod_list <- c("cluster_inds_only", "knowsizes", "bb", "negbin", "lognormal")
    stanmod_list <- c("knowsizes_centered", "knowsizes_noncentered")
    #stanmod_list <- c("bb", "negbin", "lognormal")
    use.sizes <- use_sizes
    outcome.type <- outcome_type
#    loopres <- foreach (k = 1:nrow(sim.params),
#                        .packages = c("rstan", "plyr", "tidyr", "lme4"),
#                        .export = ls(envir = globalenv()),
#                        .verbose=TRUE) %dopar% {

print(paste0("use sizes, outcome type: ", use.sizes, " ", outcome.type))

for (k in 1:nrow(sim.params)) {
      # Set parameters for this simulation
      num.clusters <- sim.params[k, "num.clusters.list"]
      num.units    <- sim.params[k, "num.units.list"]
      #use.sizes    <- sim.params[k, "use.sizes.list"]
      #outcome.type <- sim.params[k, "outcome.type.list"]

      # Print a message about which parameters we're running now
      cat("Running in parallel for", num.clusters, "clusters,", num.units,
          "units, use_sizes =", use.sizes, ", and", outcome.type, "outcomes.\n")
  
      # Sample data using above parameters
      cat("Sampling data\n")
      print(Sys.time())
      sampledata(num.clusters, num.units, use.sizes, outcome.type, rootdir, simno)
  
      # Run stan models
      results.list <- vector(mode = "list",
                             length = length(stanmod_list)*2 + 2)
                             #length = length(stanmod_list)*2 + 2)
      for (p in 1:length(stanmod_list)) {
        stanmod_name <- stanmod_list[p]
        load(paste0(rootdir, "/src/analysis/", stanmod_name, ".RData"))
#print(str(stanmod))
        cat("##################################################################################\n")
        cat("##################################################################################\n")
        print(paste0("Starting to run use.sizes = ", use.sizes,
                     ", stanmod = ", stanmod_name,
                     ", num.clusters = ", num.clusters,
                     ", num.units = ", num.units))
        print(Sys.time())
        stan_res <- runstan(num.clusters, num.units, use.sizes, rootdir, simno,
                            stanmod, stanmod_name, num.iter, num.chains)
        results.list[[2*p - 1]] <- stan_res[["par.ests"]]
        names(results.list)[2*p - 1] <- paste0("param_ests_", stanmod_name)
        results.list[[2*p]] <- stan_res[["ybar.ests"]]
        names(results.list)[2*p] <- paste0("ybar_ests_", stanmod_name)
        rm(stan_res)
        cat("##################################################################################\n")
      } # end stanmod loop
  
      # Compare using lmer
      print("Running lmer_compare")
      print(Sys.time())
      lmer_res <- lmer_compare(num.clusters, num.units, use.sizes, outcome.type, rootdir, simno)
      results.list[[2*length(stanmod_list) + 1]] <- lmer_res
      names(results.list)[2*length(stanmod_list) + 1] <- "param_ests_lmer"
      rm(lmer_res)
  
      # Estimate ybar using survey package
      print("Running svy_ests")
      print(Sys.time())
      J <- numclusters # number of clusters in the population
      svy_res <- svy_ests(J, num.clusters, num.units, use.sizes, outcome.type, rootdir, simno)
      #results.list[[2*length(stanmod_list) + 1]] <- svy_res
      #names(results.list)[2*length(stanmod_list) + 1] <- "ybar_ests_svy"
      results.list[[2*length(stanmod_list) + 2]] <- svy_res
      names(results.list)[2*length(stanmod_list) + 2] <- "ybar_ests_svy"
      rm(svy_res)

      print("##################################################################################")
      print("##################################################################################")
  
      # Save results
      if (num.units <= 1) {
        nunits <- paste(100*num.units, "pct", sep = "")
      } else {
        nunits <- num.units
      }
      saveRDS(results.list,
              paste0(rootdir, "output/simulation/results_usesizes_",
                     use.sizes, "_", outcome.type, "_nclusters_", num.clusters,
                     "_nunits_", nunits, "_simno_", simno, ".rds"))
      rm(results.list)
      # Delete the sampled data to save space
      fil1 <- paste0(rootdir, "output/simulation/simdata_usesizes_",
                     use.sizes, "_", outcome.type, "_nclusters_", num.clusters,
                     "_nunits_", nunits, "_simno_", simno, ".rds")
      if (file.exists(fil1)) {
        file.remove(fil1)
      }
#   return(NULL)
 } # end parallel loop

#  stopCluster(cl) # close workers

# Count the number of results files created by this file
  res.path <- paste0(rootdir, "output/simulation/")
  res.pattern <- paste0("results_usesizes_", use.sizes, "_", outcome.type,
                        ".*_simno_", simno, ".rds")
  sim.files <- list.files(path = res.path, pattern = res.pattern)
  if (length(sim.files) == nrow(sim.params)) {
    sink(paste0(rootdir, "output/simulation/sim_res_check_usesizes_",
                use.sizes, "_", outcome.type, "_simno_", simno, ".txt"))
    cat("Success! All", length(sim.files), "files for this simulation cleared successfully!\n")
    print(sim.files)
    sink()
  }

