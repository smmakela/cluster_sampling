#!/usr/bin/R Rscript


'
Usage: count_sim_res_check_files_R [options]

Options:
  -u --use_sizes <us_val>     Whether outcomes depend on cluster sizes
  -o --outcome_type <ot_val>  Whether outcomes are continuous or binary
  -z --size_model <sz_mod>    Model used to create pop cluster sizes
  -n --num_sims <nsims>       Total number of simulations

' -> doc

  #############################################################################
  ### Set lib paths, source files 
  #############################################################################
    libdir <- "/vega/stats/users/smm2253/rpackages"
    _libPaths(libdir)
    rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
    Sys_setenv(HOME = rootdir)
    
  #############################################################################
  ### Load libraries, get options, set up parallel stuff
  #############################################################################
    require(docopt)
    require(methods)

    # Store the docopt options as variables we can use in the code
    opts <- docopt(doc) 
    opts_names <- names(opts)

    # The options are listed twice, once with "--" in front of the option name
    # and once without, so remove the "--" ones
    opts_names <- opts_names[-grep("--", opts_names)]
    opts <- opts[opts_names]
    print(str(opts))
    print(opts_names)

    # The options are read in as strings, so make them numeric here
    for (j in 1:length(opts_names)) {
      curr_name <- opts_names[j]
      if (is_character(opts[[j]])) {
        print(paste0("Assigning ", opts[[j]], " to ", curr_name))
        assign(curr_name, opts[[j]])
      } else {
        print(paste0("Assigning ", as_numeric(opts[[j]]), " to ", curr_name))
        assign(curr_name, as_numeric(opts[[j]]))
      }
    }

  #############################################################################
  ### Count the number of sim check files for the current option values
  #############################################################################
    fil_list <- list_files(path = paste0(rootdir, "output/simulation"),
                           pattern = paste0("stan_check_usesizes_",
                                            use_sizes, "_", outcome_type, "_",
                                            size_model, "_sim.*.txt"))
    if (length(fil_list) == num_sims) {
      sink(paste0(rootdir, "output/simulation/total_stan_res_check_usesizes_",
                  use_sizes, "_", outcome_type, "_", size_model, ".txt"))
      cat("All", num_sims, "simulations finished with results files.")
      sink()
    }

  #############################################################################
  ### Check whether all files were successfully created
  #############################################################################
  # Count the number of results files created by this file
  #  res_path <- paste0(rootdir, "output/simulation/")
  #  res_pattern <- paste0("stan_results_usesizes_", use_sizes, "_", outcome_type,
  #                        "_", size_model, "_",  model_name, ".*sim_", simno, ".rds")
  #  sim_files <- list_files(path = res_path, pattern = res_pattern)
  #  if (length(sim_files) == nrow(sim_params)) {
  #    sink(paste0(rootdir, "output/simulation/sim_res_stan_check_usesizes_",
  #                use_sizes, "_", outcome_type, "_", size_model, "_",  model_name,
  #                "_sim_", simno, "_txt"), split = FALSE)
  #    cat("Success! All", length(sim_files), "files cleared successfully!\n")
  #    print(sim_files)
  #    sink()
  #  }
