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
# Get the values of the sim_parameters for this iteration
################################################################################
  sp1 <- expand.grid(use_sizes = c(0, 1),
                     outcome_type = c("binary", "continuous"),
                     size_model = c("multinomial", "poisson"),
                     model_name = c("bb", "cluster_inds_only",
                                    "knowsizes", "lognormal", "negbin"),
                     num_clusters = c(5, 10, 20, 30),
                     num_units = c(0.05, 0.1, 0.25, 0.5, 1, 10, 30, 60))
  sp2 <- expand.grid(use_sizes = c(0, 1),
                     outcome_type = c("binary", "continuous"),
                     size_model = "ff",
                     model_name = c("bb", "cluster_inds_only",
                                    "knowsizes", "lognormal", "negbin"),
                     num_clusters = 16,
                     num_units = 99)
  if (rownum <= nrow(sp1)) {
    curr_params <- sp1[rownum, ]
  } else {
    curr_params <- sp2[rownum - nrow(sp1), ]
  }

  cat("####################################################################\n")
  cat("current set of params is for rownum", rownum, "\n")
  print(curr_params)

  use_sizes <- curr_params$use_sizes
  outcome_type <- curr_params$outcome_type
  size_model <- curr_params$size_model
  model_name <- curr_params$model_name
  num_clusters <- curr_params$num_clusters
  num_units <- curr_params$num_units
  if (num_units <= 1) {
    nunits <- paste(100*num_units, "pct", sep = "")
  } else {
    nunits <- num_units
  }

  # Concatenate to get curr_nt stubs
  curr_stub <- paste0("usesizes_", use_sizes, "_", outcome_type, "_",
                      size_model, "_", model_name, "_nclusters_", num_clusters,
                      "_nunits_", nunits, "_sim_.*.rds")
  stan_stub <- paste0("stan_results_", curr_stub)

  cat("#####################################################################\n")
  cat("Currently on:", curr_stub, "\n")

################################################################################
# Calculate true values
################################################################################
  cat("#####################################################################\n")
  cat("Compiling true values\n")
  params_true <- data.frame()
  Mj_true_summ <- data.frame()
  ybar_true <- data.frame()
  popdata <- readRDS(paste0(resdir, "/popdata_usesizes_", use_sizes, "_",
                            outcome_type, "_", size_model, ".rds"))
print(str(popdata))
  Mj <- popdata[["Mj"]]

  # true parameter values
  params_true <- gather(popdata[["truepars"]], key = param_name) 
  names(params_true) <- c("param_name", "true")
  params_true$param_name <- as.character(params_true$param_name)

  # true Mj statistics
  Mj_true_summ <- data.frame(sum  = sum(Mj),
                             mean = mean(Mj),
                             sd   = sd(Mj),
                             p025 = quantile(Mj, 0.025),
                             p25  = quantile(Mj, 0.25),
                             p50  = quantile(Mj, 0.50),
                             p75  = quantile(Mj, 0.75),
                             p975 = quantile(Mj, 0.975))
  Mj_true_summ <- gather(Mj_true_summ, key = stat, value = true)
  Mj_true_summ$stat <- as.character(Mj_true_summ$stat)

  # true ybar
  ybar_true <- mean(popdata[["pop_data"]]$y)

  print("params_true:") 
  print(params_true) 
  print("Mj_true_summ:") 
  print(Mj_true_summ) 
  print("ybar_true:") 
  print(ybar_true) 
  print(Sys.time())

################################################################################
### Preallocate space to hold all results
################################################################################
  num_div_trans <- 0
  n_sims <- 100
  ybar_ests <- data.frame(true  = rep(ybar_true, n_sims),
                          mean  = rep(NA, n_sims),
                          sd    = rep(NA, n_sims),
                          p025  = rep(NA, n_sims),
                          p25   = rep(NA, n_sims),
                          p50   = rep(NA, n_sims),
                          p75   = rep(NA, n_sims),
                          p975  = rep(NA, n_sims),
                          simno = c(1:n_sims))
  Nj_ests <- expand.grid(stat  = Mj_true_summ$stat,
                         simno = c(1:n_sims))
  Nj_ests$stat <- as.character(Nj_ests$stat)
  Nj_ests <- left_join(Nj_ests, Mj_true_summ, by = "stat")
  Nj_ests$est <- NA
  stat_list <- c("mean", "sd", "p025", "p25", "p50", "p75", "p975")
  if (outcome_type == "continuous") {
    pn <- c("alpha0", "gamma0", "alpha1", "gamma1", 
            "sigma_beta0", "sigma_beta1", "sigma_y")
  } else {
    pn <- c("alpha0", "gamma0", "sigma_beta0")
  }
  if (model_name == "cluster_inds_only") { # does not use gammas
    pn <- pn[!(pn %in% c("gamma0", "gamma1"))]
  }
  param_ests <- expand.grid(param_name  = pn,
                            simno       = c(1:n_sims))
  param_ests$param_name <- as.character(param_ests$param_name)
  param_ests <- left_join(param_ests, params_true, by = c("param_name"))
  param_ests[, stat_list] <- NA
  print(head(param_ests))
  print(head(ybar_ests))
  print(head(Nj_ests))

################################################################################
### Loop through STAN result files
################################################################################
  cat("#####################################################################\n")
  cat(" Now doing STAN files\n")
  keep_param_stats <- c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%")
  fil_list <- list.files(resdir, stan_stub)
  # loop through models, outcome_type, usesizes and compile -- uses less
  # memory than trying to compile everything at once
  for (i in 1:length(fil_list)) {
    # read in curr_nt file
    curr_name <- fil_list[i]
    curr_fil <- readRDS(curr_name)

    # parse the filename to extract sim_number
    simno <- regmatches(curr_name, regexpr("[0-9]{1,3}.rds", curr_name))
    simno <- as.numeric(gsub(".rds", "", simno))
  
    # see if this was a case of divergent transitions we couldn't get rid of
    if (length(curr_fil) == 1) {
      num_div_trans <- num_div_trans + 1
      next
    }

    ### PARAM ESTS
    # pull out parameter estimates, make so that statistics are long
    tmp <- curr_fil[["par_ests"]]
    tmp$param_name <- rownames(tmp)

    # don't need to keep params like mu_star, etc
#print("BEFORE")
#print("unique(param_ests$param_name)")
#print(unique(param_ests$param_name))
#print(tmp$param_name)
#print("tmp$param_name")
#print(tmp$param_name)
    tmp <- tmp[tmp$param_name %in% unique(param_ests$param_name), ]
    tmp <- tmp[, c("param_name", keep_param_stats)]
    names(tmp) <- c("param_name", stat_list)

    # replace the appropriate part of param_ests with tmp
    inds <- param_ests$simno == simno
#cat("simno:", simno, "\n")
#print("tmp")
#print(tmp)
#print("str(tmp)")
#print(str(tmp))
#print("tmp$param_name")
#print(tmp$param_name)
#print("str(param_ests[inds, ])")
#print(str(param_ests[inds, ]))
#print("str(param_ests[inds, names(tmp)])")
#print(str(param_ests[inds, names(tmp)]))
#print("stat_list")
#print(stat_list)

    param_ests[inds, names(tmp)] <- tmp

    ### YBAR ESTS
    # pull out summaries of ybar_new
    draw.summ <- curr_fil[["draw_summ"]]
    if (nrow(draw.summ) == 1) {
      tmp <- draw.summ
    } else {
      tmp <- draw.summ[draw.summ$param == "ybar_new", ]
    } 

    # don't need the "param" column since everything is for ybar_new
    tmp$param <- NULL

    # rename "truth" to "true" to match what's in ybar_ests
    tmp$true <- tmp$truth
    tmp$truth <- NULL

    # replace the appropriate part of ybar_ests with tmp
    inds <- ybar_ests$simno == simno
#print("tmp")
#print(tmp)
#print("str(tmp)")
#print(str(tmp))
#print("str(ybar_ests[inds, ])")
#print(str(ybar_ests[inds, ]))
#print("str(ybar_ests[inds, names(tmp)])")
#print(str(ybar_ests[inds, names(tmp)]))
#print("str(ybar_ests[inds, ])")
#print(str(ybar_ests[inds, ]))
    ybar_ests[inds, names(tmp)] <- tmp
 
    ### NJ NEW
    # pull out summaries of Nj_new
    if (model_name %in% c("cluster_inds_only", "knowsizes")) {
      next # don't have Nj_new for these models
    }
    tmp <- curr_fil[["Nj_new_means"]]
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
print("tmp:")
print(tmp)
print(warnings())
    tmp <- tmp %>%
      tidyr::gather(key = stat, value = est)
print("tmp:")
print(tmp)
print(warnings())
    
    # replace the appropriate part of param_ests with tmp
    inds <- Nj_ests$simno == simno
#cat("simno:", simno, "\n")
#print("tmp:")
#print(tmp)
#print("dim(tmp):")
#print(dim(tmp))
#print("str(Nj_ests):")
#print(str(Nj_ests))
#print("Nj_ests$stat[inds]")
#print(Nj_ests$stat[inds])
#print("dim(Nj_ests[inds])")
#print(dim(Nj_ests[inds,]))
#print("tmp$stat:")
#print(tmp$stat)
    Nj_ests$stat[inds]       <- tmp$stat
#print("HERE")
    Nj_ests$est[inds]        <- tmp$est
  } # end file for loop

  # Print number of times we couln't get rid of divergent transitions
  cat("Number of times we couldn't get rid of divergent transitions:",
      num_div_trans, "\n")

  # Now collapse across all sim_ to get means, CIs, etc
  ### PARAM ESTS
  param_ests_summ <- param_ests %>%
    dplyr::group_by(param_name) %>%
    dplyr::summarise(num_sims   = sum(!is.na(mean)),
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

  param_ests_summ$use_sizes <- use_sizes
  param_ests_summ$outcome_type <- outcome_type
  param_ests_summ$outcome_type <- as.character(param_ests_summ$outcome_type)
  param_ests_summ$size_model <- size_model
  param_ests_summ$size_model <- as.character(param_ests_summ$size_model)
  param_ests_summ$num_clusters <- num_clusters
  param_ests_summ$num_units <- num_units
  param_ests_summ$model_name <- model_name
  param_ests_summ$model_name <- as.character(param_ests_summ$model_name)
  param_ests_summ$param_name <- as.character(param_ests_summ$param_name)

  print(param_ests_summ)
  print(warnings())
  print(Sys.time())

  ### YBAR ESTS
  ybar_ests_summ <- ybar_ests %>%
    dplyr::summarise(num_sims   = sum(!is.na(mean)),
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

  ybar_ests_summ$use_sizes <- use_sizes
  ybar_ests_summ$outcome_type <- outcome_type
  ybar_ests_summ$outcome_type <- as.character(ybar_ests_summ$outcome_type)
  ybar_ests_summ$size_model <- size_model
  ybar_ests_summ$size_model <- as.character(ybar_ests_summ$size_model)
  ybar_ests_summ$num_clusters <- num_clusters
  ybar_ests_summ$num_units <- num_units
  ybar_ests_summ$model_name <- model_name
  ybar_ests_summ$model_name <- as.character(ybar_ests_summ$model_name)

  print(ybar_ests_summ)
  print(warnings())
  print(Sys.time())

  ### NJ ESTS
  Nj_ests_summ <- Nj_ests %>%
    dplyr::group_by(stat) %>%
    dplyr::summarise(num_sims = sum(!is.na(est)),
                     bias     = mean(true - est, na.rm = TRUE),
                     rel_bias = mean((true - est)/true, na.rm = TRUE),
                     rmse     = sqrt(mean((true - est)^2, na.rm = TRUE)),
                     rel_rmse = sqrt(mean(((true - est)/true)^2, na.rm = TRUE)))

  Nj_ests_summ$use_sizes <- use_sizes
  Nj_ests_summ$outcome_type <- outcome_type
  Nj_ests_summ$outcome_type <- as.character(Nj_ests_summ$outcome_type)
  Nj_ests_summ$size_model <- size_model
  Nj_ests_summ$size_model <- as.character(Nj_ests_summ$size_model)
  Nj_ests_summ$num_clusters <- num_clusters
  Nj_ests_summ$num_units <- num_units
  Nj_ests_summ$model_name <- model_name
  Nj_ests_summ$model_name <- as.character(Nj_ests_summ$model_name)
  Nj_ests_summ$stat <- as.character(Nj_ests_summ$stat)

  print(Nj_ests_summ)
  print(warnings())
  print(Sys.time())

  print(nrow(param_ests_summ))
  print(nrow(ybar_ests_summ))
  print(nrow(Nj_ests_summ))

  res <- list(param_ests_summ = param_ests_summ,
              ybar_ests_summ = ybar_ests_summ,
              Nj_ests_summ = Nj_ests_summ)
  res.fil <- paste0("compiled_stan_results_usesizes_", use_sizes, "_",
                    outcome_type, "_", size_model, "_", model_name,
                    "_nclusters_", num_clusters, 
                    "_nunits_", nunits, "_", today, ".rds")
  saveRDS(res, file = paste0(resdir, res.fil))
  cat("filename:", res.fil, "\n")
  rm(res)


