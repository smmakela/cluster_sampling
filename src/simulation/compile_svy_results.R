# Author: Susanna Makela
# Date: 13 Jan 2016
# Purpose: process results from simulation

################################################################################
### Get the row of sim_params that we are on
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
  outcome_type <- curr_params$outcome_type
  size_model <- curr_params$size_model
  num_clusters <- curr_params$num_clusters
  num_units <- curr_params$num_units
  if (num_units <= 1) {
    nunits <- paste(100*num_units, "pct", sep = "")
  } else {
    nunits <- num_units
  }

  # Concatenate to get current stubs
  curr_stub <- paste0("svy_ests_usesizes_", use_sizes, "_", outcome_type,
                      "_", size_model, "_nclusters_", num_clusters, "_nunits_",
                      nunits, "_simno_.*.rds")

  cat("#####################################################################\n")
  cat("Currently on:", curr_stub, "\n")

################################################################################
# Calculate true values
################################################################################
  cat("#####################################################################\n")
  cat("Compiling true values\n")
  ybar_true <- data.frame()
  popdata <- readRDS(paste0(resdir, "/popdata_usesizes_", use_sizes, "_",
                            outcome_type, "_", size_model, ".rds"))

  # true ybar
  ybar_true <- mean(popdata[["pop_data"]]$y)

  print("ybar_true:") 
  print(ybar_true) 
  print(Sys.time())

################################################################################
### Preallocate space to hold all results
################################################################################
  n_sims <- 100
  if (outcome_type == "continuous") {
    model_list <- c("hajek", "greg")
  } else {
    model_list <- "hajek"
  }
  ybar_ests <- expand.grid(model_name = model_list,
                           simno = c(1:n_sims))
  ybar_ests$model_name <- as.character(ybar_ests$model_name)
  stat_list <- c("ybar_hat", "ybar_se", "ybar_hat_lci50", "ybar_hat_uci50",
                 "ybar_hat_lci95", "ybar_hat_uci95")
  ybar_ests[, c(stat_list, "ybar_true")] <- NA
  ybar_ests <- ybar_ests[, c("ybar_true", "model_name", "simno", stat_list)]
  print(head(ybar_ests))

################################################################################
### Loop through SVY result files
################################################################################
  cat("#####################################################################\n")
  cat(" Now doing SVY files\n")
  ci_50_mult <- qnorm(0.25, lower.tail = FALSE)
  fil_list <- list.files(resdir, curr_stub)
  # loop through models, outcome type, usesizes and compile -- uses less
  # memory than trying to compile everything at once
  for (i in 1:length(fil_list)) {
    # read in current file
    curr_name <- fil_list[i]
    curr_fil <- readRDS(curr_name)

    # if this is a binary one, remove the greg row since it'll all be NA
    curr_fil$model_name <- as.character(curr_fil$model_name)
    curr_fil <- curr_fil[curr_fil$model_name %in% model_list, ]

    # parse the filename to extract sim number
    simno <- regmatches(curr_name, regexpr("[0-9]{1,3}.rds", curr_name))
    simno <- as.numeric(gsub(".rds", "", simno))

    # fix the 50% intervals -- they are actually 68% intervals...
    curr_fil$ybar_hat_lci50 <- curr_fil$ybar_hat - ci_50_mult * curr_fil$ybar_se
    curr_fil$ybar_hat_uci50 <- curr_fil$ybar_hat + ci_50_mult * curr_fil$ybar_se
  
    # replace the appropriate part of ybar_ests with tmp
    inds <- ybar_ests$simno == simno
    ybar_ests[inds, names(curr_fil)] <- curr_fil
 
  } # end file for loop

  ### YBAR ESTS
  ybar_ests_summ <- ybar_ests %>%
    dplyr::group_by(model_name) %>%
    dplyr::summarise(num.sims   = sum(!is.na(ybar_hat)),
                     bias       = mean(ybar_true - ybar_hat, na.rm = TRUE),
                     rel_bias   = mean((ybar_true - ybar_hat)/ybar_true, na.rm = TRUE),
                     rmse       = sqrt(mean((ybar_true - ybar_hat)^2, na.rm = TRUE)),
                     rel_rmse   = sqrt(mean(((ybar_true - ybar_hat)/ybar_true)^2, na.rm = TRUE)),
                     covg_50    = mean(ybar_hat_lci50 <= ybar_true &
                                       ybar_true <= ybar_hat_uci50, na.rm = TRUE),
                     covg_95    = mean(ybar_hat_lci95 <= ybar_true &
                                     ybar_true <= ybar_hat_uci95, na.rm = TRUE),
                     len_50     = mean(ybar_hat_uci50 - ybar_hat_lci50, na.rm = TRUE),
                     len_95     = mean(ybar_hat_uci95 - ybar_hat_lci95, na.rm = TRUE),
                     rel_len_50 = mean((ybar_hat_uci50 - ybar_hat_lci50)/ybar_true, na.rm = TRUE),
                     rel_len_95 = mean((ybar_hat_uci95 - ybar_hat_lci95)/ybar_true, na.rm = TRUE))

  ybar_ests_summ$use_sizes <- use_sizes
  ybar_ests_summ$outcome_type <- outcome_type
  ybar_ests_summ$outcome_type <- as.character(ybar_ests_summ$outcome_type)
  ybar_ests_summ$size_model <- size_model
  ybar_ests_summ$size_model <- as.character(ybar_ests_summ$size_model)
  ybar_ests_summ$num_clusters <- num_clusters
  ybar_ests_summ$num_units <- num_units
  ybar_ests_summ$model_name <- as.character(ybar_ests_summ$model_name)

  print(ybar_ests_summ)
  print(warnings())
  print(Sys.time())


  res_fil <- paste0("compiled_svy_results_usesizes_", use_sizes, "_",
                    outcome_type, "_", size_model, "_nclusters_", num_clusters, 
                    "_nunits_", nunits, "_", today, ".rds")
  saveRDS(ybar_ests_summ, file = paste0(resdir, res_fil))


