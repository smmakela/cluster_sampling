#!/usr/bin/R Rscript


'
Usage: count_stan_check_files.R --use_sizes=<us_val> --outcome_type=<ot_val> --num_sims=<nsims> --model_name<mod_name>

Options:
  -u --use_sizes <us_val>     Whether outcomes depend on cluster sizes
  -o --outcome_type <ot_val>  Whether outcomes are continuous or binary
  -n --num_sims <nsims>       Total number of simulations
  -m --model_name <mod_name>  Name of stan model

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
                                            use_sizes, "_", outcome_type,
                                            model_name, "_sim.*.txt"))
    if (length(fil_list) == num_sims) {
      sink(paste0(rootdir, "output/simulation/total_stan_check_usesizes_",
                  use_sizes, "_", outcome_type, "_", model_name, ".txt"))
      cat("All", num_sims, "simulations finished with results files_")
      sink()
    }

