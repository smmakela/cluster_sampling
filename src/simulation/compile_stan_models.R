#!/usr/bin/R Rscript
# Purpose: quick script to compile stan models

'
Usage: compile_stan_models.R  --model_name=<modname> 

Options:
  -m --model_name <modname>     Name of model to compile 

' -> doc

  # Main project directory
  rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
  .libPaths("/vega/stats/users/smm2253/rpackages")

  require(docopt)
  require(methods)
  require(rstan)

  # Store the docopt options as variables we can use in the code
  opts <- docopt(doc) 
  opts.names <- names(opts)

  # The options are listed twice, once with "--" in front of the option name
  # and once without, so remove the "--" ones
  opts.names <- opts.names[-grep("--", opts.names)]
  opts <- opts[opts.names]
  print(str(opts))
  print(opts.names)

  # compile, and save
  cat("*******compiling ", opts$model_name, "*************\n")
  stanmod <- stan_model(file = paste0(rootdir, "/src/analysis/",
                                      opts$model_name, ".stan"))
  saveRDS(stanmod, file = paste0(rootdir, "/src/analysis/", opts$model_name, ".rds")) 
  

