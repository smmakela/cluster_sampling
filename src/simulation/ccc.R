#!/usr/bin/R Rscript

'
Usage: ccc.R  --use_sizes=<us_val> --outcome_type=<ot_val> --mod_name=<m_name>

Options:
  -u --use_sizes <us_val>      Whether outcomes depend on cluster sizes
  -o --outcome_type <ot_val>   Whether outcomes are continuous or binary
  -m --mod_name <m_name>       Which model to run

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
  ### Load libraries, parse docopt arguments
  #############################################################################
    require(docopt)
    require(methods)
    require(plyr)
    require(dplyr)
    require(tidyr)
    require(lme4)
    require(rstan)
    require(survey)

    # Store the docopt options as variables we can use in the code
    opts <- docopt(doc) 
    opts.names <- names(opts)

    # The options are listed twice, once with "--" in front of the option name
    # and once without, so remove the "--" ones
    opts.names <- opts.names[-grep("--", opts.names)]
    opts <- opts[opts.names]
    opts.names <- gsub("_", ".", opts.names)
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

    # add values for other important parameters
    num.clusters <- 999
    num.units <- 999
    simno <- 999
    num.iter <- 2500
    num.chains <- 4
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())

  #############################################################################
  ### Run stan/lmer
  #############################################################################
    sampledata(num.clusters, num.units, use.sizes, outcome.type, rootdir, simno)
    if (mod.name != "lmer") {
      load(paste0(rootdir, "/src/analysis/", mod.name, ".RData"))
      cat("##################################################################################\n")
      cat("##################################################################################\n")
      cat("Starting to run use.sizes =", use.sizes, ", stanmod = ", mod.name,
          "for", num.iter, "iterations and", num.chains, "chains.\n")
      print(Sys.time())
      stan_res <- runstan(num.clusters, num.units, use.sizes, rootdir, simno,
                          stanmod, mod.name, num.iter, num.chains)
      saveRDS(stan_res,
              paste0(rootdir, "output/simulation/ccc_usesizes_",
                     use.sizes, "_", outcome.type, "_", mod.name, ".rds"))
      rm(stan_res)
    } else { 
      print("Running lmer_compare")
      print(Sys.time())
      lmer_res <- lmer_compare(num.clusters, num.units, use.sizes, outcome.type,
                               rootdir, simno)
      saveRDS(lmer_res,
              paste0(rootdir, "output/simulation/ccc_usesizes_",
                     use.sizes, "_", outcome.type, "_lmer.rds"))
      rm(lmer_res)
    }
    cat("##################################################################################\n")
    cat("##################################################################################\n")


