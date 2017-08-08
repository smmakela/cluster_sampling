#!/usr/bin/R Rscript

'
Usage: sim_master_stan.R [options]

Options:
  -s --simno <simnumber>      Current simulation number
  -u --use_sizes <us_val>     Whether outcomes depend on cluster sizes
  -o --outcome_type <ot_val>  Whether outcomes are continuous or binary
  -z --size_model <sz_mod>    Model used in creating cluster sizes
  -n --numclusters <J>        Number of clusters in population [default: 100]
  -m --model_name <modname>   Name of stan model being run

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

  #############################################################################
  ### Load libraries, get options, set up parallel stuff
  #############################################################################
    require(docopt)
    require(methods)
    require(plyr)
    require(dplyr)
    require(tidyr)
    require(pps)
    require(rstan)
    require(sampling)

    # Store the docopt options as variables we can use in the code
    opts <- docopt(doc) 
    opts.names <- names(opts)

    # The options are listed twice, once with "--" in front of the option name
    # and once without, so remove the "--" ones
    print(opts.names)
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

  #############################################################################
  ### Create list of parameters to loop through for sim
  #############################################################################
    #num_clusters.list <- c(5, 10, 30, 50)
    #num_units.list <- c(0.05, 0.1, 0.25, 0.5, 1, 10, 30, 60)
    num_clusters.list <- 10
    num_units.list <- 0.1
    if (outcome_type == "binary") {
      stanmod_name <- paste0(model_name, "_binary")
    } else {
      stanmod_name <- model_name
    }
    tmplist <- list(num_clusters.list = num_clusters.list,
                    num_units.list    = num_units.list)
    sim.params <- expand.grid(tmplist)

  #############################################################################
  ### Loop through cluster/unit lists
  #############################################################################
    cat("##################################################################################\n")
    cat("##################################################################################\n")
    cat("THIS FILE: use_sizes =", use_sizes, ", outcome_type =", outcome_type,
        ", size_model =", size_model, ", stanmod =", stanmod_name, "\n")
    print(Sys.time())

    if (grepl("ff", size_model)) {
      num_clusters <- 16
      num_units <- 99
      #if (size_model == "ffstrat") {
      #  stanmod_name <- paste0(stanmod_name, "_strat")
      #}

      # Print a message about which parameters we're running now
      cat("CURRENTLY ON:", 
          "outcome_type =", outcome_type,
          ", use_sizes =", use_sizes,
          ", size_model =", size_model,
          ", stanmod =", stanmod_name, "\n")
  
      # If the file exists, skip it
      #nunits <- num_units
      #fil <- paste0(rootdir, "output/simulation/stan_results_usesizes_",
      #              use_sizes, "_", outcome_type, "_", size_model, "_",
      #              model_name, "_nclusters_", num_clusters,
      #              "_nunits_", nunits, "_sim_", simno, ".rds")
      #if (file.exists(fil)) {
      #  cat("----------------------------------------------------\n")
      #  cat("This file already exists, next!\n")
      #  cat(fil, "\n")
      #  cat("----------------------------------------------------\n")
      #  next
      #}

      # Sample data using above parameters
      cat("Sampling data\n")
      print(Sys.time())
      sim_data <- sampledata(J = 77, num_clusters, num_units, use_sizes,
                             outcome_type, size_model)
      cat("DONE sampling\n")
saveRDS(sim_data, file = paste0(rootdir, "/output/simulation/simdataTEMP_", use_sizes, "_", outcome_type,
                          "_", size_model, "_", model_name, "_sim_", simno, ".rds")) 
 
      # Run stan model
      cat("Running stan\n")
      print(Sys.time())
      stanmod <- readRDS(paste0(rootdir, "/src/analysis/", stanmod_name, ".rds"))
      print(Sys.time())
      stan_res <- runstan(num_clusters, num_units, use_sizes, outcome_type,
                          size_model, rootdir, simno, stanmod, stanmod_name,
                          sim_data, num_iter = 2000, num_chains = 4)
      cat("##################################################################################\n")
      print(warnings()) 
    } else {
      for (k in 1:nrow(sim.params)) {
        # Set parameters for this simulation
        num_clusters <- sim.params[k, "num_clusters.list"]
        num_units    <- sim.params[k, "num_units.list"]

# TEMP: removing to redo multinomial since we want all sims to be based on the same pop
        # If the file exists, skip it
        #if (num_units <= 1) {
        #  nunits <- paste(100*num_units, "pct", sep = "")
        #} else {
        #  nunits <- num_units
        #}
        #fil <- paste0(rootdir, "output/simulation/stan_results_usesizes_",
        #              use_sizes, "_", outcome_type, "_", size_model, "_",
        #              model_name, "_nclusters_", num_clusters,
        #              "_nunits_", nunits, "_sim_", simno, ".rds")
        #if (file.exists(fil)) {
        #  cat("----------------------------------------------------\n")
        #  cat("This file already exists, next!\n")
        #  cat(fil, "\n")
        #  cat("----------------------------------------------------\n")
        #  next
        #}

        # Print a message about which parameters we're running now
        cat("CURRENTLY ON:", 
            "num_clusters =", num_clusters,
            ", num_units =", num_units,
            ", use_sizes =", use_sizes,
            ", size_model =", size_model,
            ", stanmod =", stanmod_name, "\n")
  
        # Sample data using above parameters
        cat("Sampling data\n")
        print(Sys.time())
        sim_data <- sampledata(J = numclusters, num_clusters, num_units, use_sizes,
                               outcome_type, size_model)
        cat("DONE sampling\n")
saveRDS(sim_data, file = paste0(rootdir, "/output/simulation/simdataTEMP_", use_sizes, "_", outcome_type,
                          "_", size_model, "_", model_name, "_sim_", simno, ".rds")) 
        # Run stan model
        cat("Running stan\n")
        print(Sys.time())
        stanmod <- readRDS(paste0(rootdir, "/src/analysis/", stanmod_name, ".rds"))
        print(Sys.time())
        stan_res <- runstan(num_clusters, num_units, use_sizes, outcome_type,
                            size_model, rootdir, simno, stanmod, stanmod_name,
                            sim_data, num_iter = 2000, num_chains = 4)
        cat("##################################################################################\n")
        print(warnings()) 
      } # end sim parameters loop
    } # end if size_model = "ff" statement

  #############################################################################
  ### Check whether all files were successfully created
  #############################################################################
  # Count the number of results files created by this file
    res.path <- paste0(rootdir, "output/simulation/")
    res.pattern <- paste0("stan_results_usesizes_", use_sizes, "_", outcome_type,
                          "_", size_model, "_", model_name, ".*_sim_", simno, ".rds")
    sim.files <- list.files(path = res.path, pattern = res.pattern)
    cat("length(sim.files)=", length(sim.files), "\n")
    if ((!grepl("ff", size_model) & length(sim.files) == nrow(sim.params)) | 
        (grepl("ff", size_model) & length(sim.files) == 1)) {
      sink(paste0(rootdir, "output/simulation/stan_check_usesizes_",
                  use_sizes, "_", outcome_type, "_", size_model, "_", model_name,
                  "_sim_", simno, ".txt"), split = FALSE)
      cat("Success! All", length(sim.files), "files cleared successfully!\n")
      print(sim.files)
      sink()
    }


