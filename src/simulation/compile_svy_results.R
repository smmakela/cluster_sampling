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
  sim.params <- expand.grid(use.sizes = c(0, 1),
                            outcome.type = c("binary", "continuous"),
                            num.clusters = c(5, 10, 20, 30),
                            num.units = c(0.05, 0.1, 0.25, 0.5, 1, 10, 30, 60))
print(nrow(sim.params))
  curr.params <- sim.params[rownum, ]
  use.sizes <- curr.params$use.sizes
  outcome.type <- curr.params$outcome.type
  num.clusters <- curr.params$num.clusters
  num.units <- curr.params$num.units
  if (num.units <= 1) {
    nunits <- paste(100*num.units, "pct", sep = "")
  } else {
    nunits <- num.units
  }

  # Concatenate to get current stubs
  curr.stub <- paste0("svy_ests_usesizes_", use.sizes, "_", outcome.type,
                      "_nclusters_", num.clusters, "_nunits_",
                      nunits, "_simno_.*.rds")

  cat("#####################################################################\n")
  cat("Currently on:", curr.stub, "\n")

################################################################################
# Calculate true values
################################################################################
  cat("#####################################################################\n")
  cat("Compiling true values\n")
  ybar.true <- data.frame()
  popdata <- readRDS(paste0(resdir, "/popdata_usesizes_", use.sizes, "_",
                            outcome.type, ".rds"))

  # true ybar
  ybar.true <- mean(popdata[["pop.data"]]$y)

  print("ybar.true:") 
  print(ybar.true) 
  print(Sys.time())

################################################################################
### Preallocate space to hold all results
################################################################################
  n.sims <- 100
  if (outcome.type == "continuous") {
    model.list <- c("hajek", "greg")
  } else {
    model.list <- "hajek"
  }
  ybar.ests <- expand.grid(model.name = model.list,
                           simno = c(1:n.sims))
  ybar.ests$model.name <- as.character(ybar.ests$model.name)
  stat.list <- c("ybar_hat", "ybar_se", "ybar_hat_lci50", "ybar_hat_uci50",
                 "ybar_hat_lci95", "ybar_hat_uci95")
  ybar.ests[, c(stat.list, "ybar.true")] <- NA
  ybar.ests <- ybar.ests[, c("ybar.true", "model.name", "simno", stat.list)]
  print(head(ybar.ests))

################################################################################
### Loop through SVY result files
################################################################################
  cat("#####################################################################\n")
  cat(" Now doing SVY files\n")
  ci_50_mult <- qnorm(0.25, lower.tail = FALSE)
  fil.list <- list.files(resdir, curr.stub)
  # loop through models, outcome type, usesizes and compile -- uses less
  # memory than trying to compile everything at once
  for (i in 1:length(fil.list)) {
    # read in current file
    curr.name <- fil.list[i]
    curr.fil <- readRDS(curr.name)

    # if this is a binary one, remove the greg row since it'll all be NA
    curr.fil$model.name <- as.character(curr.fil$model.name)
    curr.fil <- curr.fil[curr.fil$model.name %in% model.list, ]

    # parse the filename to extract sim number
    simno <- regmatches(curr.name, regexpr("[0-9]{1,3}.rds", curr.name))
    simno <- as.numeric(gsub(".rds", "", simno))

    # fix the 50% intervals -- they are actually 68% intervals...
    curr.fil$ybar_hat_lci50 <- curr.fil$ybar_hat - ci_50_mult * curr.fil$ybar_se
    curr.fil$ybar_hat_uci50 <- curr.fil$ybar_hat + ci_50_mult * curr.fil$ybar_se
  
    # replace the appropriate part of ybar.ests with tmp
    inds <- ybar.ests$simno == simno
    ybar.ests[inds, names(curr.fil)] <- curr.fil
 
  } # end file for loop

  ### YBAR ESTS
  ybar.ests.summ <- ybar.ests %>%
    dplyr::group_by(model.name) %>%
    dplyr::summarise(num.sims   = sum(!is.na(ybar_hat)),
                     bias       = mean(ybar.true - ybar_hat, na.rm = TRUE),
                     rel_bias   = mean((ybar.true - ybar_hat)/ybar.true, na.rm = TRUE),
                     rmse       = sqrt(mean((ybar.true - ybar_hat)^2, na.rm = TRUE)),
                     rel_rmse   = sqrt(mean(((ybar.true - ybar_hat)/ybar.true)^2, na.rm = TRUE)),
                     covg_50    = mean(ybar_hat_lci50 <= ybar.true &
                                       ybar.true <= ybar_hat_uci50, na.rm = TRUE),
                     covg_95    = mean(ybar_hat_lci95 <= ybar.true &
                                     ybar.true <= ybar_hat_uci95, na.rm = TRUE),
                     len_50     = mean(ybar_hat_uci50 - ybar_hat_lci50, na.rm = TRUE),
                     len_95     = mean(ybar_hat_uci95 - ybar_hat_lci95, na.rm = TRUE),
                     rel_len_50 = mean((ybar_hat_uci50 - ybar_hat_lci50)/ybar.true, na.rm = TRUE),
                     rel_len_95 = mean((ybar_hat_uci95 - ybar_hat_lci95)/ybar.true, na.rm = TRUE))

  ybar.ests.summ$use.sizes <- use.sizes
  ybar.ests.summ$outcome.type <- outcome.type
  ybar.ests.summ$outcome.type <- as.character(ybar.ests.summ$outcome.type)
  ybar.ests.summ$num.clusters <- num.clusters
  ybar.ests.summ$num.units <- num.units
  ybar.ests.summ$model.name <- as.character(ybar.ests.summ$model.name)

  print(ybar.ests.summ)
  print(warnings())
  print(Sys.time())


  res.fil <- paste0("compiled_svy_results_usesizes_", use.sizes, "_",
                    outcome.type, "_nclusters_", num.clusters, 
                    "_nunits_", nunits, "_", today, ".rds")
  saveRDS(ybar.ests.summ, file = paste0(resdir, res.fil))


