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
                            model.name = c("bb", "cluster_inds_only",
                                           "knowsizes", "lognormal", "negbin"),
                            num.clusters = c(5, 10, 20, 30),
                            num.units = c(0.05, 0.1, 0.25, 0.5, 1, 10, 30, 60))
  curr.params <- sim.params[rownum, ]

  cat("####################################################################\n")
  cat("current set of params is for rownum", rownum, "\n")
  print(curr.params)

  use.sizes <- curr.params$use.sizes
  outcome.type <- curr.params$outcome.type
  model.name <- curr.params$model.name
  num.clusters <- curr.params$num.clusters
  num.units <- curr.params$num.units
  if (num.units <= 1) {
    nunits <- paste(100*num.units, "pct", sep = "")
  } else {
    nunits <- num.units
  }

  # Concatenate to get current stubs
  curr.stub <- paste0("usesizes_", use.sizes, "_", outcome.type, "_",
                      model.name, "_nclusters_", num.clusters, "_nunits_",
                      nunits, "_sim_.*.rds")
  stan.stub <- paste0("stan_results_", curr.stub)

  cat("#####################################################################\n")
  cat("Currently on:", curr.stub, "\n")

################################################################################
# Calculate true values
################################################################################
  cat("#####################################################################\n")
  cat("Compiling true values\n")
  params.true <- data.frame()
  Mj.true.summ <- data.frame()
  ybar.true <- data.frame()
  popdata <- readRDS(paste0(resdir, "/popdata_usesizes_", use.sizes, "_",
                            outcome.type, ".rds"))
  Mj <- popdata[["Mj"]]

  # true parameter values
  params.true <- gather(popdata[["truepars"]], key = param.name) 
  names(params.true) <- c("param.name", "true")
  params.true$param.name <- as.character(params.true$param.name)

  # true Mj statistics
  Mj.true.summ <- data.frame(sum  = sum(Mj),
                             mean = mean(Mj),
                             sd   = sd(Mj),
                             p025 = quantile(Mj, 0.025),
                             p25  = quantile(Mj, 0.25),
                             p50  = quantile(Mj, 0.50),
                             p75  = quantile(Mj, 0.75),
                             p975 = quantile(Mj, 0.975))
  Mj.true.summ <- gather(Mj.true.summ, key = stat, value = true)
  Mj.true.summ$stat <- as.character(Mj.true.summ$stat)

  # true ybar
  ybar.true <- mean(popdata[["pop.data"]]$y)

  print("params.true:") 
  print(params.true) 
  print("Mj.true.summ:") 
  print(Mj.true.summ) 
  print("ybar.true:") 
  print(ybar.true) 
  print(Sys.time())

################################################################################
### Preallocate space to hold all results
################################################################################
  num.div.trans <- 0
  n.sims <- 100
  ybar.ests <- data.frame(true  = rep(ybar.true, n.sims),
                          mean  = rep(NA, n.sims),
                          sd    = rep(NA, n.sims),
                          p025  = rep(NA, n.sims),
                          p25   = rep(NA, n.sims),
                          p50   = rep(NA, n.sims),
                          p75   = rep(NA, n.sims),
                          p975  = rep(NA, n.sims),
                          simno = c(1:n.sims))
  Nj.ests <- expand.grid(stat  = Mj.true.summ$stat,
                         simno = c(1:n.sims))
  Nj.ests$stat <- as.character(Nj.ests$stat)
  Nj.ests <- left_join(Nj.ests, Mj.true.summ, by = "stat")
  Nj.ests$est <- NA
  stat.list <- c("mean", "sd", "p025", "p25", "p50", "p75", "p975")
  if (outcome.type == "continuous") {
    pn <- c("alpha0", "gamma0", "alpha1", "gamma1", 
            "sigma_beta0", "sigma_beta1", "sigma_y")
  } else {
    pn <- c("alpha0", "gamma0", "sigma_beta0")
  }
  if (model.name == "cluster_inds_only") { # does not use gammas
    pn <- pn[!(pn %in% c("gamma0", "gamma1"))]
  }
  param.ests <- expand.grid(param.name  = pn,
                            simno       = c(1:n.sims))
  param.ests$param.name <- as.character(param.ests$param.name)
  param.ests <- left_join(param.ests, params.true, by = c("param.name"))
  param.ests[, stat.list] <- NA
  print(head(param.ests))
  print(head(ybar.ests))
  print(head(Nj.ests))

################################################################################
### Loop through STAN result files
################################################################################
  cat("#####################################################################\n")
  cat(" Now doing STAN files\n")
  keep.param.stats <- c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%")
  fil.list <- list.files(resdir, stan.stub)
  # loop through models, outcome type, usesizes and compile -- uses less
  # memory than trying to compile everything at once
  for (i in 1:length(fil.list)) {
    # read in current file
    curr.name <- fil.list[i]
    curr.fil <- readRDS(curr.name)

    # parse the filename to extract sim number
    simno <- regmatches(curr.name, regexpr("[0-9]{1,3}.rds", curr.name))
    simno <- as.numeric(gsub(".rds", "", simno))
  
    # see if this was a case of divergent transitions we couldn't get rid of
    if (length(curr.fil) == 1) {
      num.div.trans <- num.div.trans + 1
      next
    }

    ### PARAM ESTS
    # pull out parameter estimates, make so that statistics are long
    tmp <- curr.fil[["par_ests"]]
    tmp$param.name <- rownames(tmp)

    # don't need to keep params like mu_star, etc
#print("BEFORE")
#print("unique(param.ests$param.name)")
#print(unique(param.ests$param.name))
#print(tmp$param.name)
#print("tmp$param.name")
#print(tmp$param.name)
    tmp <- tmp[tmp$param.name %in% unique(param.ests$param.name), ]
    tmp <- tmp[, c("param.name", keep.param.stats)]
    names(tmp) <- c("param.name", stat.list)

    # replace the appropriate part of param.ests with tmp
    inds <- param.ests$simno == simno
#cat("simno:", simno, "\n")
#print("tmp")
#print(tmp)
#print("str(tmp)")
#print(str(tmp))
#print("tmp$param.name")
#print(tmp$param.name)
#print("str(param.ests[inds, ])")
#print(str(param.ests[inds, ]))
#print("str(param.ests[inds, names(tmp)])")
#print(str(param.ests[inds, names(tmp)]))
#print("stat.list")
#print(stat.list)

    param.ests[inds, names(tmp)] <- tmp

    ### YBAR ESTS
    # pull out summaries of ybar_new
    draw.summ <- curr.fil[["draw_summ"]]
    if (nrow(draw.summ) == 1) {
      tmp <- draw.summ
    } else {
      tmp <- draw.summ[draw.summ$param == "ybar_new", ]
    } 

    # don't need the "param" column since everything is for ybar_new
    tmp$param <- NULL

    # rename "truth" to "true" to match what's in ybar.ests
    tmp$true <- tmp$truth
    tmp$truth <- NULL

    # replace the appropriate part of ybar.ests with tmp
    inds <- ybar.ests$simno == simno
#print("tmp")
#print(tmp)
#print("str(tmp)")
#print(str(tmp))
#print("str(ybar.ests[inds, ])")
#print(str(ybar.ests[inds, ]))
#print("str(ybar.ests[inds, names(tmp)])")
#print(str(ybar.ests[inds, names(tmp)]))
#print("str(ybar.ests[inds, ])")
#print(str(ybar.ests[inds, ]))
    ybar.ests[inds, names(tmp)] <- tmp
 
    ### NJ NEW
    # pull out summaries of Nj_new
    if (model.name %in% c("cluster_inds_only", "knowsizes")) {
      next # don't have Nj_new for these models
    }
    tmp <- curr.fil[["Nj_new_means"]]
    tmp <- tmp %>%
      dplyr::select(Nj_new) %>%
      dplyr::summarise(sum  = sum(Nj_new),
                       mean = mean(Nj_new),
                       sd   = sd(Nj_new),
                       p025 = quantile(Nj_new, 0.025),
                       p25  = quantile(Nj_new, 0.25),
                       p50  = quantile(Nj_new, 0.50),
                       p75  = quantile(Nj_new, 0.75),
                       p975 = quantile(Nj_new, 0.975))
    tmp <- tmp %>%
      tidyr::gather(key = stat, value = est)
    
    # replace the appropriate part of param.ests with tmp
    inds <- Nj.ests$simno == simno
#cat("simno:", simno, "\n")
#print("tmp:")
#print(tmp)
#print("dim(tmp):")
#print(dim(tmp))
#print("str(Nj.ests):")
#print(str(Nj.ests))
#print("Nj.ests$stat[inds]")
#print(Nj.ests$stat[inds])
#print("dim(Nj.ests[inds])")
#print(dim(Nj.ests[inds,]))
#print("tmp$stat:")
#print(tmp$stat)
    Nj.ests$stat[inds]       <- tmp$stat
#print("HERE")
    Nj.ests$est[inds]        <- tmp$est
  } # end file for loop

  # Print number of times we couln't get rid of divergent transitions
  cat("Number of times we couldn't get rid of divergent transitions:",
      num.div.trans, "\n")

  # Now collapse across all sims to get means, CIs, etc
  ### PARAM ESTS
  param.ests.summ <- param.ests %>%
    dplyr::group_by(param.name) %>%
    dplyr::summarise(num.sims   = sum(!is.na(mean)),
                     bias       = mean(true - mean, na.rm = TRUE),
                     rel_bias   = mean((true - mean)/true, na.rm = TRUE),
                     rmse       = sqrt(mean((true - mean)^2, na.rm = TRUE)),
                     rel_rmse   = sqrt(mean(((true - mean)/true)^2, na.rm = TRUE)),
                     covg_50    = mean(p25 <= true & true <= p75, na.rm = TRUE),
                     covg_95    = mean(p025 <= true & true <= p975, na.rm = TRUE),
                     len_50     = mean(p75 - p25, na.rm = TRUE),
                     len_95     = mean(p975 - p025, na.rm = TRUE),
                     rel_len_50 = mean((p75 - p25)/true, na.rm = TRUE),
                     rel_len_95 = mean((p975 - p025)/true, na.rm = TRUE))

  param.ests.summ$use.sizes <- use.sizes
  param.ests.summ$outcome.type <- outcome.type
  param.ests.summ$outcome.type <- as.character(param.ests.summ$outcome.type)
  param.ests.summ$num.clusters <- num.clusters
  param.ests.summ$num.units <- num.units
  param.ests.summ$model.name <- model.name
  param.ests.summ$model.name <- as.character(param.ests.summ$model.name)
  param.ests.summ$param.name <- as.character(param.ests.summ$param.name)

  print(param.ests.summ)
  print(warnings())
  print(Sys.time())

  ### YBAR ESTS
  ybar.ests.summ <- ybar.ests %>%
    dplyr::summarise(num.sims   = sum(!is.na(mean)),
                     bias       = mean(true - mean, na.rm = TRUE),
                     rel_bias   = mean((true - mean)/true, na.rm = TRUE),
                     rmse       = sqrt(mean((true - mean)^2, na.rm = TRUE)),
                     rel_rmse   = sqrt(mean(((true - mean)/true)^2, na.rm = TRUE)),
                     covg_50    = mean(p25 <= true & true <= p75, na.rm = TRUE),
                     covg_95    = mean(p025 <= true & true <= p975, na.rm = TRUE),
                     len_50     = mean(p75 - p25, na.rm = TRUE),
                     len_95     = mean(p975 - p025, na.rm = TRUE),
                     rel_len_50 = mean((p75 - p25)/true, na.rm = TRUE),
                     rel_len_95 = mean((p975 - p025)/true, na.rm = TRUE))

  ybar.ests.summ$use.sizes <- use.sizes
  ybar.ests.summ$outcome.type <- outcome.type
  ybar.ests.summ$outcome.type <- as.character(ybar.ests.summ$outcome.type)
  ybar.ests.summ$num.clusters <- num.clusters
  ybar.ests.summ$num.units <- num.units
  ybar.ests.summ$model.name <- model.name
  ybar.ests.summ$model.name <- as.character(ybar.ests.summ$model.name)

  print(ybar.ests.summ)
  print(warnings())
  print(Sys.time())

  ### NJ ESTS
  Nj.ests.summ <- Nj.ests %>%
    dplyr::group_by(stat) %>%
    dplyr::summarise(num.sims = sum(!is.na(est)),
                     bias     = mean(true - est, na.rm = TRUE),
                     rel_bias = mean((true - est)/true, na.rm = TRUE),
                     rmse     = sqrt(mean((true - est)^2, na.rm = TRUE)),
                     rel_rmse = sqrt(mean(((true - est)/true)^2, na.rm = TRUE)))

  Nj.ests.summ$use.sizes <- use.sizes
  Nj.ests.summ$outcome.type <- outcome.type
  Nj.ests.summ$outcome.type <- as.character(Nj.ests.summ$outcome.type)
  Nj.ests.summ$num.clusters <- num.clusters
  Nj.ests.summ$num.units <- num.units
  Nj.ests.summ$model.name <- model.name
  Nj.ests.summ$model.name <- as.character(Nj.ests.summ$model.name)
  Nj.ests.summ$stat <- as.character(Nj.ests.summ$stat)

  print(Nj.ests.summ)
  print(warnings())
  print(Sys.time())

  print(nrow(param.ests.summ))
  print(nrow(ybar.ests.summ))
  print(nrow(Nj.ests.summ))

  res <- list(param.ests.summ = param.ests.summ,
              ybar.ests.summ = ybar.ests.summ,
              Nj.ests.summ = Nj.ests.summ)
  res.fil <- paste0("compiled_stan_results_usesizes_", use.sizes, "_",
                    outcome.type, "_", model.name, "_nclusters_", num.clusters, 
                    "_nunits_", nunits, "_", today, ".rds")
  saveRDS(res, file = paste0(resdir, res.fil))
  cat("filename:", res.fil, "\n")
  rm(res)


