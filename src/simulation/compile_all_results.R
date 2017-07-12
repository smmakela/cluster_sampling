# Author: Susanna Makela
# Date: 13 Jan 2016
# Purpose: process results from simulation


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
# Loop through sim params
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
  tot_rows <- nrow(sp1) + nrow(sp2)
 
  # model names for stan
  model_list = c("bb", "cluster_inds_only", "knowsizes", "lognormal", "negbin")

  # allocate space for all results (columns calculated from outputs of the
  # individual compile files) 
  lmer_all <- data.frame(matrix(NA, nrow = 12*tot_rows, ncol = 18))
  svy_all <- data.frame(matrix(NA, nrow = 2*tot_rows, ncol = 17))
  stan_pars_all <- data.frame(matrix(NA, nrow = 5*7*tot_rows, ncol = 18))
  stan_ybar_all <- data.frame(matrix(NA, nrow = 5*tot_rows, ncol = 17))
  stan_Nj_all <- data.frame(matrix(NA, nrow = 5*8*tot_rows, ncol = 12))

  # counters for rows
  start_lmer <- 1
  start_svy <- 1
  start_pars <- 1
  start_ybar <- 1
  start_Nj <- 1

  for (i in 1:tot_rows) {
    if (i <= nrow(sp1)) {
      curr_params <- sp1[i, ]
    } else {
      curr_params <- sp2[i - nrow(sp1), ]
    }
cat("curr_params for i =", i, "\n")
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

    # concatenate to get current stubs
    curr_stub <- paste0("usesizes_", use_sizes, "_", outcome_type, "_",
                        size_model, "_nclusters_", num_clusters, "_nunits_",
                        nunits, "_", today, ".rds")
    cat("curr_stub:", curr_stub, "\n")
    svy_stub <- paste0("compiled_svy_results_", curr_stub)
    lmer_stub <- paste0("compiled_lmer_results_", curr_stub)

    # load files
    if (!file.exists(paste0(resdir, svy_stub))) {
      cat("This file does not exist!\n")
      cat(svy_stub, "\n")
      next
    }
    curr_svy <- readRDS(paste0(resdir, svy_stub))
    if (!file.exists(paste0(resdir, lmer_stub))) {
      cat("This file does not exist!\n")
      cat(lmer_stub, "\n")
      next
    }
    curr_lmer <- readRDS(paste0(resdir, lmer_stub))

    # rename *_all variables
    names(svy_all) <- names(curr_svy)
    names(lmer_all) <- names(curr_lmer)

#print("str(curr_svy)")
#print(str(curr_svy))
#print("str(curr_lmer)")
#print(str(curr_lmer))
#print("nrow(curr_svy)")
#print(nrow(curr_svy))
#print("length svy")
#print(length(start_svy:(start_svy+nrow(curr_svy)-1)))
#print("length lmer")
#print(length(start_lmer:(start_lmer+nrow(curr_lmer)-1)))
#print("nrow(curr_lmer)")
#print(nrow(curr_lmer))
#print("start_svy")
#print(start_svy)
#print("start_lmer")
#print(start_lmer)
#print("str(svy_all)")
#print(str(svy_all))
#print("str(lmer_all)")
#print(str(lmer_all))

    # add to output
    svy_all[start_svy:(start_svy+nrow(curr_svy)-1), ] <- curr_svy
    lmer_all[start_lmer:(start_lmer+nrow(curr_lmer)-1), ] <- curr_lmer

    # update start values
    start_svy <- start_svy + nrow(curr_svy) + 1
    start_lmer <- start_lmer + nrow(curr_lmer) + 1

    # loop through models to do stan files
    for (m in model_list) {
      stan_stub <- paste0("compiled_stan_results_usesizes_", use_sizes, "_",
                          outcome_type, "_", size_model, "_", m, "_nclusters_",
                          num_clusters, "_nunits_", nunits, "_", today, ".rds")
      if (!file.exists(paste0(resdir, stan_stub))) {
        cat("This file does not exist!\n")
        cat(stan_stub, "\n")
        next
      }
      stan_res <- readRDS(paste0(resdir, stan_stub))

      # pull summaries out of the list
      curr_pars <- stan_res[["param_ests_summ"]]
      curr_ybar <- stan_res[["ybar_ests_summ"]]
      curr_Nj <- stan_res[["Nj_ests_summ"]]

#print("str(curr_pars)")
#print(str(curr_pars))
#print("str(curr_ybar)")
#print(str(curr_ybar))
#print("str(curr_Nj)")
#print(str(curr_Nj))
#print("start_pars")
#print(start_pars)
#print("start_ybar")
#print(start_ybar)
#print("start_Nj")
#print(start_Nj)

      # rename *_all variables
      names(stan_pars_all) <- names(curr_pars)
      names(stan_ybar_all) <- names(curr_ybar)
      names(stan_Nj_all) <- names(curr_Nj)

      # add to output
      stan_pars_all[start_pars:(start_pars+nrow(curr_pars)-1), ] <- curr_pars
      stan_ybar_all[start_ybar:(start_ybar+nrow(curr_ybar)-1), ] <- curr_ybar
      stan_Nj_all[start_Nj:(start_Nj+nrow(curr_Nj)-1), ] <- curr_Nj

      # update start values
      start_pars <- start_pars + nrow(curr_pars) + 1
      start_ybar <- start_ybar + nrow(curr_ybar) + 1
      start_Nj <- start_Nj + nrow(curr_Nj) + 1
    } # end stan model loop
  } # end sim_params loop
print("str(svy_all)")
print(str(svy_all))
print("str(lmer_all)")
print(str(lmer_all))
print("str(stan_pars_all)")
print(str(stan_pars_all))
print("str(stan_ybar_all)")
print(str(stan_ybar_all))
print("str(stan_Nj_all)")
print(str(stan_Nj_all))


################################################################################
# Get rid of remaining all-NA rows
################################################################################
  count_na <- function(x) sum(is.na(x))

  ncol_svy <- ncol(svy_all)
  svy_all <- svy_all %>%
    dplyr::mutate(num_na = apply(., 1, count_na)) %>%
    dplyr::filter(!(num_na == ncol_svy)) %>%
    dplyr::select(-num_na)
  
  ncol_lmer <- ncol(lmer_all)
  lmer_all <- lmer_all %>%
    dplyr::mutate(num_na = apply(., 1, count_na)) %>%
    dplyr::filter(!(num_na == ncol_lmer)) %>%
    dplyr::select(-num_na)

  ncol_pars <- ncol(stan_pars_all)
  stan_pars_all <- stan_pars_all %>%
    dplyr::mutate(num_na = apply(., 1, count_na)) %>%
    dplyr::filter(!(num_na == ncol_pars)) %>%
    dplyr::select(-num_na)

  ncol_ybar <- ncol(stan_ybar_all)
  stan_ybar_all <- stan_ybar_all %>%
    dplyr::mutate(num_na = apply(., 1, count_na)) %>%
    dplyr::filter(!(num_na == ncol_ybar)) %>%
    dplyr::select(-num_na)

  ncol_Nj <- ncol(stan_Nj_all)
  stan_Nj_all <- stan_Nj_all %>%
    dplyr::mutate(num_na = apply(., 1, count_na)) %>%
    dplyr::filter(!(num_na == ncol_Nj)) %>%
    dplyr::select(-num_na)

################################################################################
# SAVE
################################################################################

  saveRDS(svy_all, paste0(resdir, "all_svy_results_", today, ".rds"))
  saveRDS(lmer_all, paste0(resdir, "all_lmer_results_", today, ".rds"))
  saveRDS(stan_pars_all, paste0(resdir, "all_pars_results_", today, ".rds"))
  saveRDS(stan_ybar_all, paste0(resdir, "all_ybar_results_", today, ".rds"))
  saveRDS(stan_Nj_all, paste0(resdir, "all_Nj_results_", today, ".rds"))


