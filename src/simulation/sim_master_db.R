#!/usr/bin/R Rscript

'
Usage: sim_master_db.R  --simno=<simnumber> --use_sizes=<us_val> --outcome_type=<ot_val> [options]

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

  #############################################################################
  ### Create list of parameters to loop through for sim
  #############################################################################
    num.clusters.list <- c(5, 10, 20, 30)
    num.units.list <- c(0.05, 0.1, 0.25, 0.5, 1, 10, 30, 60)
    tmplist <- list(num.clusters.list = num.clusters.list,
                    num.units.list    = num.units.list)
    sim.params <- expand.grid(tmplist)

  #############################################################################
  ### Loop through cluster/unit lists
  #############################################################################
    outcome.type <- outcome_type
    use.sizes <- use_sizes
    print(paste0("use sizes, outcome type: ", use.sizes, " ", outcome.type))


    for (k in 1:nrow(sim.params)) {
      num.clusters <- sim.params[k, "num.clusters.list"]
      num.units    <- sim.params[k, "num.units.list"]

      # If the file exists, skip it
      if (num.units <= 1) {
        nunits <- paste(100*num.units, "pct", sep = "")
      } else {
        nunits <- num.units
      }
      fil1 <- paste0(rootdir, "output/simulation/lmer_ests_usesizes_",
                     use.sizes, "_", outcome.type, "_nclusters_", num.clusters,
                     "_nunits_", nunits, "_simno_", simno, ".rds")
      fil2 <- paste0(rootdir, "output/simulation/svy_ests_usesizes_",
                     use.sizes, "_", outcome.type, "_nclusters_", num.clusters,
                     "_nunits_", nunits, "_simno_", simno, ".rds")
      if (file.exists(fil1) & file.exists(fil2)) {
        cat("----------------------------------------------------\n")
        cat("Both files already exist, next!\n")
        cat(fil1, "\n")
        cat("----------------------------------------------------\n")
        next
      }

      # Print a message about which parameters we're running now
      cat("Running design-based models for", num.clusters, "clusters,",
          num.units, "units, use_sizes =", use.sizes, ", and",
          outcome.type, "outcomes.\n")
  
      # Sample data using above parameters
      cat("Sampling data\n")
      print(Sys.time())
      sim.data <- sampledata(num.clusters, num.units, use.sizes, outcome.type)
      cat("DONE sampling\n")
 
      # Compare parameter estimates from lmer
      print("Running lmer_compare")
      print(Sys.time())
      lmer_res <- lmer_compare(num.clusters, num.units, use.sizes,
                               outcome.type, sim.data, rootdir, simno)
      print(warnings()) 
  
      # Estimate ybar using survey package
      print("Running svy_ests")
      print(Sys.time())
      J <- numclusters # number of clusters in the population
      svy_res <- svy_ests(J, num.clusters, num.units, use.sizes, outcome.type,
                          sim.data, rootdir, simno)
      print(warnings()) 

      print("##################################################################################")
      print("##################################################################################")
    } # end sim parameters loop

  #############################################################################
  ### Check whether all files were successfully created
  #############################################################################
  # Count the number of results files created by this file
    res.path <- paste0(rootdir, "output/simulation/")
    res.pattern <- paste0(".*_ests_usesizes_", use.sizes, "_", outcome.type,
                          ".*_sim_", simno, ".rds")
    sim.files <- list.files(path = res.path, pattern = res.pattern)
print(sim.files)
    if (length(sim.files) == 2*nrow(sim.params)) {
      sink(paste0(rootdir, "output/simulation/db_check_usesizes_",
                  use.sizes, "_", outcome.type, "_sim_", simno, ".txt"),
           split = FALSE)
      cat("Success! All", length(sim.files), "files cleared successfully!\n")
      print(sim.files)
      sink()
    }
    print(warnings())

