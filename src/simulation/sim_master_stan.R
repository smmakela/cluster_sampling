#!/usr/bin/R Rscript

'
Usage: sim_master_stan.R  --simno=<simnumber> --use_sizes=<us_val> --outcome_type=<ot_val> --numclusters=<J> --model_name=<modname>

Options:
  -s --simno <simnumber>      Current simulation number
  -u --use_sizes <us_val>     Whether outcomes depend on cluster sizes
  -o --outcome_type <ot_val>  Whether outcomes are continuous or binary
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
    num.clusters.list <- c(5, 10, 20, 30)
    num.units.list <- c(0.05, 0.1, 0.25, 0.5, 1, 10, 30, 60)
    if (outcome_type == "binary") {
      stanmod_name <- paste0(model_name, "_binary")
    } else {
      stanmod_name <- model_name
    }
    tmplist <- list(num.clusters.list = num.clusters.list,
                    num.units.list    = num.units.list)
    sim.params <- expand.grid(tmplist)

  #############################################################################
  ### Loop through cluster/unit lists
  #############################################################################
    outcome.type <- outcome_type
    use.sizes <- use_sizes
    cat("##################################################################################\n")
    cat("##################################################################################\n")
    cat("THIS FILE: use.sizes =", use.sizes, ", stanmod =", stanmod_name, "\n")
    print(Sys.time())

    for (k in 1:nrow(sim.params)) {
      # Set parameters for this simulation
      num.clusters <- sim.params[k, "num.clusters.list"]
      num.units    <- sim.params[k, "num.units.list"]

      # If the file exists, skip it
      if (num.units <= 1) {
        nunits <- paste(100*num.units, "pct", sep = "")
      } else {
        nunits <- num.units
      }
      fil <- paste0(rootdir, "output/simulation/stan_results_usesizes_",
                    use.sizes, "_", outcome.type, "_", model_name, 
                    "_nclusters_", num.clusters,
                    "_nunits_", nunits, "_sim_", simno, ".rds")
      if (file.exists(fil)) {
        cat("----------------------------------------------------\n")
        cat("This file already exists, next!\n")
        cat(fil, "\n")
        cat("----------------------------------------------------\n")
        next
      }

      # Print a message about which parameters we're running now
      cat("CURRENTLY ON:", 
          "num.clusters =", num.clusters,
          ", num.units =", num.units,
          "use.sizes =", use.sizes,
          ", stanmod =", stanmod_name, "\n")
  
      # Sample data using above parameters
      cat("Sampling data\n")
      print(Sys.time())
      sim.data <- sampledata(num.clusters, num.units, use.sizes, outcome.type)
      cat("DONE sampling\n")
 
      # Run stan model
      cat("Running stan\n")
      print(Sys.time())
      stanmod <- readRDS(paste0(rootdir, "/src/analysis/", stanmod_name, ".rds"))
      print(Sys.time())
      stan_res <- runstan(num.clusters, num.units, use.sizes, outcome.type,
                          rootdir, simno, stanmod, stanmod_name, sim.data,
                          num.iter = 1000, num.chains = 4)
      cat("##################################################################################\n")
      print(warnings()) 
    } # end sim parameters loop


  #############################################################################
  ### Check whether all files were successfully created
  #############################################################################
  # Count the number of results files created by this file
    res.path <- paste0(rootdir, "output/simulation/")
    res.pattern <- paste0("stan_results_usesizes_", use.sizes, "_", outcome.type,
                          "_", model_name, ".*_sim_", simno, ".rds")
    sim.files <- list.files(path = res.path, pattern = res.pattern)
    if (length(sim.files) == nrow(sim.params)) {
      sink(paste0(rootdir, "output/simulation/stan_check_usesizes_",
                  use.sizes, "_", outcome.type, "_", model_name,
                  "_sim_", simno, ".txt"), split = FALSE)
      cat("Success! All", length(sim.files), "files cleared successfully!\n")
      print(sim.files)
      sink()
    }

