# Author: Susanna Makela
# Date: 13 Jan 2016
# Purpose: process results from simulation

################################################################################
### Get the row of sim.params that we are on
################################################################################
  args <- commandArgs(FALSE)
  cat("args:", args, "\n")
  len <- length(args)
  rownum <- -1*as.numeric(args[len])
  cat("row no:", rownum, "\n")

################################################################################
### Setup of directories and libraries
################################################################################
  libdir <- "/vega/stats/users/smm2253/rpackages"
  .libPaths(libdir)
  rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
  Sys.setenv(HOME = rootdir)

  # set working directory, figure directory
  resdir <- paste0(rootdir, "output/simulation/")
  setwd(resdir)

  # libraries
  library(dplyr)
  library(tidyr)

  # print time
  print(Sys.time())
  today <- Sys.Date()
  today <- gsub("-", "_", today)

################################################################################
# Get the values of the sim parameters for this iteration
################################################################################
  sp1 <- expand.grid(use_sizes = c(0, 1),
                     outcome_type = c("binary", "continuous"),
                     size_model = c("multinomial", "poisson"),
                     num_clusters = c(5, 10, 20, 30),
                     num_units = c(0.05, 0.1, 0.25, 0.5, 1, 10, 30, 60))
  sp2 <- expand.grid(use_sizes = c(0, 1),
                     outcome_type = c("binary", "continuous"),
                     size_model = "ff",
                     num_clusters = 16,
                     num_units = 99)

  if (rownum <= nrow(sp1)) {
    curr_params <- sp1[rownum, ]
  } else {
    curr_params <- sp2[rownum - nrow(sp1), ]
  }

  use_sizes <- curr_params$use_sizes
  size_model <- curr_params$size_model
  outcome_type <- curr_params$outcome_type
  num_clusters <- curr_params$num_clusters
  num_units <- curr_params$num_units
  if (num_units <= 1) {
    nunits <- paste(100*num_units, "pct", sep = "")
  } else {
    nunits <- num_units
  }

  # Concatenate to get current stubs
  curr_stub <- paste0("lmer_ests_usesizes_", use_sizes, "_", outcome_type,
                      "_", size_model, "_nclusters_", num_clusters, "_nunits_",
                      nunits, "_simno_.*.rds")

  cat("#####################################################################\n")
  cat("Currently on:", curr_stub, "\n")

################################################################################
# Calculate true values
################################################################################
  cat("#####################################################################\n")
  cat("Compiling true values\n")
  params_true <- data.frame()
  popdata <- readRDS(paste0(resdir, "/popdata_usesizes_", use_sizes, "_",
                            outcome_type, "_", size_model, ".rds"))

  # true parameter values
  params_true <- gather(popdata[["truepars"]], key = param_name) 
  names(params_true) <- c("param_name", "true")
  params_true$param_name <- as.character(params_true$param_name)

  print("params_true:") 
  print(params_true) 
  print(Sys.time())

################################################################################
### Preallocate space to hold all results
################################################################################
  n_sims <- 100
  stat_list <- c("est", "stderr", "lci50", "uci50", "lci95", "uci95")
  param_ests <- expand.grid(param_name  = params_true$param_name[params_true$param_name != "ybar_true"],
                            model_name  = c("with_sizes", "no_sizes"),
                            simno       = c(1:n_sims))
  param_ests$param_name <- as.character(param_ests$param_name)
  tokeep <- !(param_ests$param_name %in% c("gamma0", "gamma1") &
              param_ests$model_name == "no_sizes")
  param_ests <- param_ests[tokeep, ]
  param_ests <- left_join(param_ests, params_true, by = c("param_name"))
  param_ests[, stat_list] <- NA
  print(head(param_ests))

################################################################################
### Loop through LMER result_files
################################################################################
  cat("#####################################################################\n")
  cat("Now doing LMER_files\n")
  ci_50_mult <- qnorm(0.25, lower.tail = FALSE)
  ci_95_mult <- qnorm(0.025, lower.tail = FALSE)
  fil_list <- list.files(resdir, curr_stub)
  # loop through models, outcome_type, usesizes and compile -- use_ less
  # memory than trying to compile everything at once
  for (i in 1:length(fil_list)) {
    
    # read in current_file
    curr_name <- fil_list[i]
    curr_fil <- readRDS(curr_name)

    # parse the_filename to extract sim number
    simno <- regmatches(curr_name, regexpr("[0-9]{1,3}.rds", curr_name))
    simno <- as.numeric(gsub(".rds", "", simno))

    # add param names, model names
    curr_fil$param_name <- rownames(curr_fil)
    rownames(curr_fil) <- NULL
    curr_fil$param_name <- gsub("01", "0", curr_fil$param_name)
    curr_fil$param_name <- gsub("11", "1", curr_fil$param_name)
    curr_fil$param_name <- gsub("sigma_y1", "sigma_y", curr_fil$param_name)
    curr_fil$whichmodel[curr_fil$whichmodel == 1] <- "with_sizes"
    curr_fil$whichmodel[curr_fil$whichmodel == 2] <- "no_sizes"

    # add 50pct, 95pct CI's
    curr_fil$lci50 <- curr_fil$est - ci_50_mult * curr_fil$stderr
    curr_fil$uci50 <- curr_fil$est + ci_50_mult * curr_fil$stderr
    curr_fil$lci95 <- curr_fil$est - ci_95_mult * curr_fil$stderr
    curr_fil$uci95 <- curr_fil$est + ci_95_mult * curr_fil$stderr
    names(curr_fil) <- c("est", "stderr", "true", "model_name", "param_name",
                         "lci50", "uci50", "lci95", "uci95")
    neworder <- c("param_name", "true", "model_name", "est", "stderr", 
                  "lci50", "uci50", "lci95", "uci95")
    curr_fil <- curr_fil[, neworder]

    # get the right index
    inds <- param_ests$simno == simno
    param_ests[inds, names(curr_fil)] <- curr_fil

  } # end_file for loop
  
  # Now collapse across all sims to get means, CIs, etc
  ### PARAM ESTS
  param_ests_summ <- param_ests %>%
    dplyr::group_by(param_name, model_name) %>%
    dplyr::summarise(num.sims   = sum(!is.na(est)),
                     bias       = mean(true - est, na.rm = TRUE),
                     rel_bias   = mean((true - est)/true, na.rm = TRUE),
                     rmse       = sqrt(mean((true - est)^2, na.rm = TRUE)),
                     rel_rmse   = sqrt(mean(((true - est)/true)^2, na.rm = TRUE)),
                     covg_50    = mean(lci50 <= true & true <= uci50, na.rm = TRUE),
                     covg_95    = mean(lci95 <= true & true <= uci95, na.rm = TRUE),
                     len_50     = mean(uci50 - lci50, na.rm = TRUE),
                     len_95     = mean(uci95 - lci95, na.rm = TRUE),
                     rel_len_50 = mean((uci50 - lci50)/true, na.rm = TRUE),
                     rel_len_95 = mean((uci95 - lci95)/true, na.rm = TRUE))

  param_ests_summ$use_sizes <- use_sizes
  param_ests_summ$outcome_type <- outcome_type
  param_ests_summ$outcome_type <- as.character(param_ests_summ$outcome_type)
  param_ests_summ$size_model <- size_model
  param_ests_summ$size_model <- as.character(param_ests_summ$size_model)
  param_ests_summ$num_clusters <- num_clusters
  param_ests_summ$num_units <- num_units
  param_ests_summ$model_name <- as.character(param_ests_summ$model_name)

  print(param_ests_summ)
  print(warnings())
  print(Sys.time())


  res_fil <- paste0("compiled_lmer_results_usesizes_", use_sizes, "_",
                    outcome_type, "_", size_model, "_nclusters_", num_clusters, 
                    "_nunits_", nunits, "_", today, ".rds")
  saveRDS(param_ests_summ, file = paste0(resdir, res_fil))


