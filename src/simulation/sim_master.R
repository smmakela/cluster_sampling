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
    require(sampling)

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

    # Set up parallel parameters
#    cl <- makeCluster(min(detectCores(), 10),
#                      type = "FORK",
#                      outfile = paste0(rootdir,
#                                       "/output/simulation/parallel_output_",
#                                       "usesizes_", use_sizes, "_",
#                                       outcome_type, ".txt"))
##print("makeCluster")
#print("registerDoParallel")
#    registerDoParallel(cl)
#print("clusterEvalQ")
#    clusterEvalQ(cl, .libPaths( "/vega/stats/users/smm2253/rpackages"))
#print("getDoParWorkers")
#    print(getDoParWorkers()) # make sure foreach will actually run in parallel
#print("clusterEvalQ specific libraries")
#    clusterEvalQ(cl, library(iterators, lib.loc = "/vega/stats/users/smm2253/rpackages"))
#    clusterEvalQ(cl, library(foreach, lib.loc = "/vega/stats/users/smm2253/rpackages"))
#    clusterEvalQ(cl, library(doParallel, lib.loc = "/vega/stats/users/smm2253/rpackages"))

  #############################################################################
  ### Create list of parameters to loop through for sim
  #############################################################################
    #num.clusters.list <- c(5)
    num.clusters.list <- c(5, 10, 20, 30)
    #num.units.list <- c(10)
    num.units.list <- c(0.05, 0.1, 0.25, 0.5, 1, 10, 30, 60)
    #use.sizes.list <- c(0, 1)
    #outcome.type.list <- c("continuous", "binary")
    modlist <- c("cluster_inds_only", "negbin", "lognormal", "knowsizes", "bb")
    if (outcome_type == "binary") {
      stanmod_list <- paste0(modlist, "_binary")
    } else {
      stanmod_list <- modlist
    }
    tmplist <- list(num.clusters.list = num.clusters.list,
                    num.units.list    = num.units.list,
                    stanmod.list      = stanmod_list)
                    #use.sizes.list    = use.sizes.list,
                    #outcome.type.list = outcome.type.list)
    sim.params <- expand.grid(tmplist)

  #############################################################################
  ### Loop through cluster/unit lists
  #############################################################################
    outcome.type <- outcome_type
    #modlist <- c("cluster_inds_only", "negbin", "lognormal", "knowsizes", "bb")
    #if (outcome.type == "binary") {
    #  stanmod_list <- paste0(modlist, "_binary")
    #} else {
    #  stanmod_list <- modlist
    #}
    #stanmod_list <- c("cluster_inds_only", "knowsizes", "lognormal", "negbin", "bb")
    #stanmod_list <- c("bb")
    use.sizes <- use_sizes
    print(paste0("use sizes, outcome type: ", use.sizes, " ", outcome.type))
    #loopres <- foreach (k = 1:nrow(sim.params),
    #                    .packages = c("rstan", "dplyr", "tidyr", "lme4", "sampling"),
    #                    .export = ls(envir = globalenv()),
    #                    .verbose=TRUE) %dopar% {


    for (k in 1:nrow(sim.params)) {
      # Set parameters for this simulation
      num.clusters <- sim.params[k, "num.clusters.list"]
      num.units    <- sim.params[k, "num.units.list"]
      stanmod_name    <- sim.params[k, "stanmod.list"]
      #use.sizes    <- sim.params[k, "use.sizes.list"]
      #outcome.type <- sim.params[k, "outcome.type.list"]

      # Print a message about which parameters we're running now
      cat("Running in parallel for", num.clusters, "clusters,", num.units,
          "units, use_sizes =", use.sizes, ", and", outcome.type, "outcomes.\n")
  
      # Sample data using above parameters
      cat("Sampling data\n")
      print(Sys.time())
      sampledata(num.clusters, num.units, use.sizes, outcome.type, rootdir, simno)
      cat("DONE sampling\n")
 
      # Run stan models
      #results.list <- vector(mode = "list", 20) # arbitrarily long to be easier
      #list_ctr <- 1
      for (p in 1:length(stanmod_list)) {
      #  stanmod_name <- stanmod_list[p]
        stanmod <- readRDS(paste0(rootdir, "/src/analysis/", stanmod_name, ".rds"))
        cat("##################################################################################\n")
        cat("##################################################################################\n")
        print(paste0("Starting to run use.sizes = ", use.sizes,
                     ", stanmod = ", stanmod_name,
                     ", num.clusters = ", num.clusters,
                     ", num.units = ", num.units))
        print(Sys.time())
        stan_res <- runstan(num.clusters, num.units, use.sizes, outcome.type,
                            rootdir, simno, stanmod, stanmod_name,
                            num.iter = 1000, num.chains = 4)
        #stan_res <- runstan(num.clusters, num.units, use.sizes, outcome.type,
        #                    rootdir, simno, stanmod, stanmod_name,
        #                    num.iter = 1000, num.chains = 4)
        #results.list[[list_ctr]] <- stan_res[["par_ests"]]
        #names(results.list)[list_ctr] <- paste0("param_ests_", stanmod_name)
        #list_ctr <- list_ctr + 1
        #results.list[[list_ctr]] <- stan_res[["draw_summ"]]
        #names(results.list)[list_ctr] <- paste0("draw_summ_", stanmod_name)
        #list_ctr <- list_ctr + 1
        #if (stanmod_name %in% c("lognormal", "negbin")) {
        #  results.list[[list_ctr]] <- stan_res[["Nj_new_means"]]
        #  names(results.list)[list_ctr] <- paste0("Nj_new_means_", stanmod_name)
        #  list_ctr <- list_ctr + 1
        #}
        #rm(stan_res)
        cat("##################################################################################\n")
      } # end stanmod loop
      print(warnings()) 
      # Compare using lmer
      print("Running lmer_compare")
      print(Sys.time())
      lmer_res <- lmer_compare(num.clusters, num.units, use.sizes, outcome.type, rootdir, simno)
      print(warnings()) 
      #results.list[[list_ctr]] <- lmer_res
      #names(results.list)[list_ctr] <- "param_ests_lmer"
      #list_ctr <- list_ctr + 1
      #rm(lmer_res)
  
      # Estimate ybar using survey package
      print("Running svy_ests")
      print(Sys.time())
      J <- numclusters # number of clusters in the population
      svy_res <- svy_ests(J, num.clusters, num.units, use.sizes, outcome.type, rootdir, simno)
      print(warnings()) 
      #results.list[[list_ctr]] <- svy_res
      #names(results.list)[list_ctr] <- "ybar_ests_svy"
      #rm(svy_res)

      print("##################################################################################")
      print("##################################################################################")
  
      # Save results
      #print("saving results")
      #saveRDS(results.list,
      #        paste0(rootdir, "output/simulation/results_usesizes_",
      #               use.sizes, "_", outcome.type, "_nclusters_", num.clusters,
      #               "_nunits_", nunits, "_sim_", simno, ".rds"))
      #rm(results.list)

      print("deleting sampled data")
      if (num.units <= 1) {
        nunits <- paste(100*num.units, "pct", sep = "")
      } else {
        nunits <- num.units
      }
      # Delete the sampled data to save space
      fil1 <- paste0(rootdir, "output/simulation/sampledata_usesizes_",
                     use.sizes, "_", outcome.type, "_nclusters_", num.clusters,
                     "_nunits_", nunits, "_sim_", simno, ".rds")
      if (file.exists(fil1)) {
        file.remove(fil1)
      }
      print(warnings()) 
    } # end parallel loop

    #stopCluster(cl) # close workers

  #############################################################################
  ### Check whether all files were successfully created
  #############################################################################
  # Count the number of results files created by this file
    res.path <- paste0(rootdir, "output/simulation/")
    res.pattern <- paste0("results_usesizes_", use.sizes, "_", outcome.type,
                          ".*_sim_", simno, ".rds")
    sim.files <- list.files(path = res.path, pattern = res.pattern)
print(sim.files)
    if (length(sim.files) == nrow(sim.params)) {
      sink(paste0(rootdir, "output/simulation/sim_res_check_usesizes_",
                  use.sizes, "_", outcome.type, "_sim_", simno, ".txt"),
           split = FALSE)
      cat("Success! All", length(sim.files), "files cleared successfully!\n")
      print(sim.files)
      sink()
    }

