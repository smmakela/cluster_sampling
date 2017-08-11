#!/usr/bin/R Rscript

'
Usage: sim_master_db.R  [options]

Options:
  -s --simno <simnumber>      Current simulation number
  -u --use_sizes <us_val>     Whether outcomes depend on cluster sizes
  -o --outcome_type <ot_val>  Whether outcomes are continuous or binary
  -z --size_model <sz_mod>    Model used to generate pop cluster sizes
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
    source(paste0(rootdir, "src/simulation/svy_ests_strat.R"))

  #############################################################################
  ### Load libraries, get options, set up parallel stuff
  #############################################################################
    require(docopt)
    require(methods)
    require(plyr)
    require(dplyr)
    require(tidyr)
    require(lme4)
    require(data.table)
    require(pps)
    require(survey)
    require(sampling)

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
      if (is.character(opts[[j]])) {
        print(paste0("Assigning ", opts[[j]], " to ", opts_names[j]))
        assign(opts_names[j], opts[[j]])
      } else {
        print(paste0("Assigning ", as.numeric(opts[[j]]), " to ", opts_names[j]))
        assign(opts_names[j], as.numeric(opts[[j]]))
      }
    }

  #############################################################################
  ### Create list of parameters to loop through for sim
  #############################################################################
    #num_clusters_list <- c(5, 10, 20, 30)
    #num_units_list <- c(0.05, 0.1, 0.25, 0.5, 1, 10, 30, 60)
    num_clusters_list <- 10
    num_units_list <- 0.1
    tmplist <- list(num_clusters_list = num_clusters_list,
                    num_units_list    = num_units_list)
    sim_params <- expand.grid(tmplist)

  #############################################################################
  ### Loop through cluster/unit lists
  #############################################################################
    print(paste0("use sizes, outcome type: ", use_sizes, " ", outcome_type))


    if (size_model == "ff") {
      num_clusters <- 16
      num_units <- 99

      # Print a message about which parameters we're running now
      cat("CURRENTLY ON:", 
          ", use_sizes =", use_sizes,
          ", outcome_type =", outcome_type,
          ", size_model =", size_model, "\n")
  
      # Sample data using above parameters
      cat("Sampling data\n")
      print(Sys.time())
      sim_data <- sampledata(J = numclusters, num_clusters, num_units, use_sizes,
                             outcome_type, size_model)
      cat("DONE sampling\n")
 
      # Compare parameter estimates from lmer
      print("Running lmer_compare")
      print(Sys.time())
      lmer_res <- lmer_compare(num_clusters, num_units, use_sizes, outcome_type,
                               size_model, sim_data, rootdir, simno)
      print(warnings()) 
  
      # Estimate ybar using survey package
      print("Running svy_ests")
      print(Sys.time())
      J <- 77 # number of clusters in FF data
      svy_res <- svy_ests(J, num_clusters, num_units, use_sizes, 
                          outcome_type, size_model, sim_data, rootdir, simno)
      cat("##################################################################################\n")
      print(warnings()) 
    } else {
      for (k in 1:nrow(sim_params)) {
        num_clusters <- sim_params[k, "num_clusters_list"]
        num_units    <- sim_params[k, "num_units_list"]

        # If the file exists, skip it
        if (num_units <= 1) {
          nunits <- paste(100*num_units, "pct", sep = "")
        } else {
          nunits <- num_units
        }
        #fil1 <- paste0(rootdir, "output/simulation/lmer_ests_usesizes_",
        #               use_sizes, "_", outcome_type, "_", size_model,
        #               "_nclusters_", num_clusters,
        #               "_nunits_", nunits, "_simno_", simno, ".rds")
        #fil2 <- paste0(rootdir, "output/simulation/svy_ests_usesizes_",
        #               use_sizes, "_", outcome_type, "_", size_model,
        #               "_nclusters_", num_clusters,
        #               "_nunits_", nunits, "_simno_", simno, ".rds")
        #if (file.exists(fil1) & file.exists(fil2)) {
        #  cat("----------------------------------------------------\n")
        #  cat("Both files already exist, next!\n")
        #  cat(fil1, "\n")
        #  cat("----------------------------------------------------\n")
        #  next
        #}

        # Print a message about which parameters we're running now
        cat("Running design-based models for", num_clusters, "clusters,",
            num_units, "units, use_sizes =", use_sizes, ", outcome_type =",
            outcome_type, ", and size_model=", size_model, "\n")
  
        # Sample data using above parameters
        cat("Sampling data\n")
        print(Sys.time())
        sim_data <- sampledata(J = numclusters, num_clusters, num_units, use_sizes,
                               outcome_type, size_model)
        cat("DONE sampling\n")
 
        # Compare parameter estimates from lmer
        print("Running lmer_compare")
        print(Sys.time())
        lmer_res <- lmer_compare(num_clusters, num_units, use_sizes, outcome_type,
                                 size_model, sim_data, rootdir, simno)
        print(warnings()) 
  
        # Estimate ybar using survey package
        print("Running svy_ests")
        print(Sys.time())
        J <- numclusters # number of clusters in the population
        svy_res <- svy_ests(J, num_clusters, num_units, use_sizes, 
                            outcome_type, size_model, sim_data, rootdir, simno)
        print(warnings()) 

        print("##################################################################################")
        print("##################################################################################")
      } # end sim parameters loop
    } # end if size_model = "ff" statement

  #############################################################################
  ### Check whether all files were successfully created
  #############################################################################
  # Count the number of results files created by this file
    res_path <- paste0(rootdir, "output/simulation/")
    res_pattern <- paste0("_*_ests_usesizes_", use_sizes, "_", outcome_type,
                          "_", size_model, "_*_sim_", simno, ".rds")
    sim_files <- list.files(path = res_path, pattern = res_pattern)
print(sim_files)
    if (length(sim_files) == 2*nrow(sim_params)) {
      sink(paste0(rootdir, "output/simulation/db_check_usesizes_",
                  use_sizes, "_", outcome_type, "_", size_model, "_sim_",
                  simno, ".txt"),
           split = FALSE)
      cat("Success! All", length(sim_files), "files cleared successfully!\n")
      print(sim_files)
      sink()
    }
    print(warnings())

