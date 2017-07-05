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
  curr.stub <- paste0("lmer_ests_usesizes_", use.sizes, "_", outcome.type,
                      "_nclusters_", num.clusters, "_nunits_",
                      nunits, "_simno_.*.rds")

  cat("#####################################################################\n")
  cat("Currently on:", curr.stub, "\n")

################################################################################
# Calculate true values
################################################################################
  cat("#####################################################################\n")
  cat("Compiling true values\n")
  params.true <- data.frame()
  popdata <- readRDS(paste0(resdir, "/popdata_usesizes_", use.sizes, "_",
                            outcome.type, ".rds"))

  # true parameter values
  params.true <- gather(popdata[["truepars"]], key = param.name) 
  names(params.true) <- c("param.name", "true")
  params.true$param.name <- as.character(params.true$param.name)

  print("params.true:") 
  print(params.true) 
  print(Sys.time())

################################################################################
### Preallocate space to hold all results
################################################################################
  n.sims <- 100
  stat.list <- c("est", "stderr", "lci50", "uci50", "lci95", "uci95")
  param.ests <- expand.grid(param.name  = params.true$param.name[params.true$param.name != "ybar_true"],
                            model.name  = c("with_sizes", "no_sizes"),
                            simno       = c(1:n.sims))
  param.ests$param.name <- as.character(param.ests$param.name)
  tokeep <- !(param.ests$param.name %in% c("gamma0", "gamma1") &
              param.ests$model.name == "no_sizes")
  param.ests <- param.ests[tokeep, ]
  param.ests <- left_join(param.ests, params.true, by = c("param.name"))
  param.ests[, stat.list] <- NA
  print(head(param.ests))

################################################################################
### Loop through LMER result files
################################################################################
  cat("#####################################################################\n")
  cat("Now doing LMER files\n")
  ci_50_mult <- qnorm(0.25, lower.tail = FALSE)
  ci_95_mult <- qnorm(0.025, lower.tail = FALSE)
  fil.list <- list.files(resdir, curr.stub)
  # loop through models, outcome type, usesizes and compile -- uses less
  # memory than trying to compile everything at once
  #for (i in 1:length(fil.list)) {
  for (i in 1:1) {
    
# read in current file
    curr.name <- fil.list[i]
    curr.fil <- readRDS(curr.name)

    # parse the filename to extract sim number
    simno <- regmatches(curr.name, regexpr("[0-9]{1,3}.rds", curr.name))
    simno <- as.numeric(gsub(".rds", "", simno))

    # add param names, model names
    curr.fil$param.name <- rownames(curr.fil)
    rownames(curr.fil) <- NULL
    curr.fil$param.name <- gsub("01", "0", curr.fil$param.name)
    curr.fil$param.name <- gsub("11", "1", curr.fil$param.name)
    curr.fil$param.name <- gsub("sigma_y1", "sigma_y", curr.fil$param.name)
    curr.fil$whichmodel[curr.fil$whichmodel == 1] <- "with_sizes"
    curr.fil$whichmodel[curr.fil$whichmodel == 2] <- "no_sizes"

    # add 50pct, 95pct CI's
    curr.fil$lci50 <- curr.fil$est - ci_50_mult * curr.fil$stderr
    curr.fil$uci50 <- curr.fil$est + ci_50_mult * curr.fil$stderr
    curr.fil$lci95 <- curr.fil$est - ci_95_mult * curr.fil$stderr
    curr.fil$uci95 <- curr.fil$est + ci_95_mult * curr.fil$stderr
    names(curr.fil) <- c("est", "stderr", "true", "model.name", "param.name",
                         "lci50", "uci50", "lci95", "uci95")
    neworder <- c("param.name", "true", "model.name", "est", "stderr", 
                  "lci50", "uci50", "lci95", "uci95")
    curr.fil <- curr.fil[, neworder]

    # get the right index
    inds <- param.ests$simno == simno
    param.ests[inds, names(curr.fil)] <- curr.fil

  } # end file for loop
  
  # Now collapse across all sims to get means, CIs, etc
  ### PARAM ESTS
  param.ests.summ <- param.ests %>%
    dplyr::group_by(param.name, model.name) %>%
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

  param.ests.summ$use.sizes <- use.sizes
  param.ests.summ$outcome.type <- outcome.type
  param.ests.summ$outcome.type <- as.character(param.ests.summ$outcome.type)
  param.ests.summ$num.clusters <- num.clusters
  param.ests.summ$num.units <- num.units
  param.ests.summ$model.name <- as.character(param.ests.summ$model.name)

  print(param.ests.summ)
  print(warnings())
  print(Sys.time())


  res.fil <- paste0("compiled_lmer_results_usesizes_", use.sizes, "_",
                    outcome.type, "_nclusters_", num.clusters, 
                    "_nunits_", nunits, "_", today, ".rds")
  saveRDS(param.ests.summ, file = paste0(resdir, res.fil))


