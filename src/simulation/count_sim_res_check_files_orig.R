#!/usr/bin/R Rscript


'
Usage: count_sim_res_check_files.R --use_sizes=<us_val> --outcome_type=<ot_val> --num_sims=<nsims> [options]

Options:
  -u --use_sizes <us_val>     Whether outcomes depend on cluster sizes
  -o --outcome_type <ot_val>  Whether outcomes are continuous or binary
  -n --num_sims <nsims>       Total number of simulations

' -> doc

  #############################################################################
  ### Set lib paths, source files 
  #############################################################################
    libdir <- "/vega/stats/users/smm2253/rpackages"
    .libPaths(libdir)
    rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
    Sys.setenv(HOME = rootdir)
    
  #############################################################################
  ### Load libraries, get options, set up parallel stuff
  #############################################################################
    require(docopt)
    require(methods)

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
      curr.name <- opts.names[j]
      curr.name <- gsub("_", ".", curr.name)
      if (is.character(opts[[j]])) {
        print(paste0("Assigning ", opts[[j]], " to ", curr.name))
        assign(curr.name, opts[[j]])
      } else {
        print(paste0("Assigning ", as.numeric(opts[[j]]), " to ", curr.name))
        assign(curr.name, as.numeric(opts[[j]]))
      }
    }

  #############################################################################
  ### Count the number of sim check files for the current option values
  #############################################################################
    fil.list <- list.files(path = paste0(rootdir, "output/simulation"),
                           pattern = paste0("sim_res_check_usesizes_",
                                            use.sizes, "_", outcome.type,
                                            "_simno_.*.txt"))
    if (length(fil.list) == num.sims) {
      sink(paste0(rootdir, "output/simulation/total_sim_res_check_usesizes_",
                  use.sizes, "_", outcome.type, ".txt"))
      cat("All", num.sims, "simulations finished with results files.")
      sink()
    }

  #############################################################################
  ### Check whether all files were successfully created
  #############################################################################
  # Count the number of results files created by this file
    res.path <- paste0(rootdir, "output/simulation/")
    res.pattern <- paste0("stan_results_usesizes_", use.sizes, "_", outcome.type,
                          "_", model_name, ".*_sim_", simno, ".rds")
    sim.files <- list.files(path = res.path, pattern = res.pattern)
    if (length(sim.files) == nrow(sim.params)) {
      sink(paste0(rootdir, "output/simulation/sim_res_stan_check_usesizes_",
                  use.sizes, "_", outcome.type, "_", model_name,
                  "_sim_", simno, ".txt"), split = FALSE)
      cat("Success! All", length(sim.files), "files cleared successfully!\n")
      print(sim.files)
      sink()
    }
