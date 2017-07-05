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
  #today <- Sys.Date()
  #today <- gsub("-", "_", today)
  today <- "2017_05_24"

################################################################################
# Loop through sim params
################################################################################
  sim.params <- expand.grid(use.sizes = c(0, 1),
                            outcome.type = c("binary", "continuous"),
                            num.clusters = c(5, 10, 20, 30),
                            num.units = c(0.05, 0.1, 0.25, 0.5, 1, 10, 30, 60))
 
  # model names for stan
  model.list = c("bb", "cluster_inds_only", "knowsizes", "lognormal", "negbin")

  # allocate space for all results (columns calculated from outputs of the
  # individual compile files) 
  lmer.all <- data.frame(matrix(NA, nrow = 12*nrow(sim.params), ncol = 17))
  svy.all <- data.frame(matrix(NA, nrow = 2*nrow(sim.params), ncol = 16))
  stan.pars.all <- data.frame(matrix(NA, nrow = 5*7*nrow(sim.params), ncol = 17))
  stan.ybar.all <- data.frame(matrix(NA, nrow = 5*nrow(sim.params), ncol = 16))
  stan.Nj.all <- data.frame(matrix(NA, nrow = 5*8*nrow(sim.params), ncol = 11))

  # counters for rows
  start.lmer <- 1
  start.svy <- 1
  start.pars <- 1
  start.ybar <- 1
  start.Nj <- 1

  for (i in 1:nrow(sim.params)) {
    curr.params <- sim.params[i, ]
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

    # concatenate to get current stubs
    curr.stub <- paste0("usesizes_", use.sizes, "_", outcome.type,
                        "_nclusters_", num.clusters, "_nunits_",
                        nunits, "_", today, ".rds")
    cat("curr.stub:", curr.stub, "\n")
    svy.stub <- paste0("compiled_svy_results_", curr.stub)
    lmer.stub <- paste0("compiled_lmer_results_", curr.stub)

    # load files
    curr.svy <- readRDS(paste0(resdir, svy.stub))
    curr.lmer <- readRDS(paste0(resdir, lmer.stub))

    # rename *.all variables
    names(svy.all) <- names(curr.svy)
    names(lmer.all) <- names(curr.lmer)

#print("str(curr.svy)")
#print(str(curr.svy))
#print("str(curr.lmer)")
#print(str(curr.lmer))
#print("nrow(curr.svy)")
#print(nrow(curr.svy))
#print("length svy")
#print(length(start.svy:(start.svy+nrow(curr.svy)-1)))
#print("length lmer")
#print(length(start.lmer:(start.lmer+nrow(curr.lmer)-1)))
#print("nrow(curr.lmer)")
#print(nrow(curr.lmer))
#print("start.svy")
#print(start.svy)
#print("start.lmer")
#print(start.lmer)
#print("str(svy.all)")
#print(str(svy.all))
#print("str(lmer.all)")
#print(str(lmer.all))

    # add to output
    svy.all[start.svy:(start.svy+nrow(curr.svy)-1), ] <- curr.svy
    lmer.all[start.lmer:(start.lmer+nrow(curr.lmer)-1), ] <- curr.lmer

    # update start values
    start.svy <- start.svy + nrow(curr.svy) + 1
    start.lmer <- start.lmer + nrow(curr.lmer) + 1

    # loop through models to do stan files
    for (m in model.list) {
      stan.stub <- paste0("compiled_stan_results_usesizes_", use.sizes, "_",
                          outcome.type, "_", m, "_nclusters_",
                          num.clusters, "_nunits_", nunits, "_", today, ".rds")
      stan.res <- readRDS(paste0(resdir, stan.stub))

      # pull summaries out of the list
      curr.pars <- stan.res[["param.ests.summ"]]
      curr.ybar <- stan.res[["ybar.ests.summ"]]
      curr.Nj <- stan.res[["Nj.ests.summ"]]

#print("str(curr.pars)")
#print(str(curr.pars))
#print("str(curr.ybar)")
#print(str(curr.ybar))
#print("str(curr.Nj)")
#print(str(curr.Nj))
#print("start.pars")
#print(start.pars)
#print("start.ybar")
#print(start.ybar)
#print("start.Nj")
#print(start.Nj)

      # rename *.all variables
      names(stan.pars.all) <- names(curr.pars)
      names(stan.ybar.all) <- names(curr.ybar)
      names(stan.Nj.all) <- names(curr.Nj)

      # add to output
      stan.pars.all[start.pars:(start.pars+nrow(curr.pars)-1), ] <- curr.pars
      stan.ybar.all[start.ybar:(start.ybar+nrow(curr.ybar)-1), ] <- curr.ybar
      stan.Nj.all[start.Nj:(start.Nj+nrow(curr.Nj)-1), ] <- curr.Nj

      # update start values
      start.pars <- start.pars + nrow(curr.pars) + 1
      start.ybar <- start.ybar + nrow(curr.ybar) + 1
      start.Nj <- start.Nj + nrow(curr.Nj) + 1
    } # end stan model loop
  } # end sim.params loop
print("str(svy.all)")
print(str(svy.all))
print("str(lmer.all)")
print(str(lmer.all))
print("str(stan.pars.all)")
print(str(stan.pars.all))
print("str(stan.ybar.all)")
print(str(stan.ybar.all))
print("str(stan.Nj.all)")
print(str(stan.Nj.all))


################################################################################
# Get rid of remaining all-NA rows
################################################################################
  count_na <- function(x) sum(is.na(x))

  ncol.svy <- ncol(svy.all)
  svy.all <- svy.all %>%
    dplyr::mutate(num_na = apply(., 1, count_na)) %>%
    dplyr::filter(!(num_na == ncol.svy)) %>%
    dplyr::select(-num_na)
  
  ncol.lmer <- ncol(lmer.all)
  lmer.all <- lmer.all %>%
    dplyr::mutate(num_na = apply(., 1, count_na)) %>%
    dplyr::filter(!(num_na == ncol.lmer)) %>%
    dplyr::select(-num_na)

  ncol.pars <- ncol(stan.pars.all)
  stan.pars.all <- stan.pars.all %>%
    dplyr::mutate(num_na = apply(., 1, count_na)) %>%
    dplyr::filter(!(num_na == ncol.pars)) %>%
    dplyr::select(-num_na)

  ncol.ybar <- ncol(stan.ybar.all)
  stan.ybar.all <- stan.ybar.all %>%
    dplyr::mutate(num_na = apply(., 1, count_na)) %>%
    dplyr::filter(!(num_na == ncol.ybar)) %>%
    dplyr::select(-num_na)

  ncol.Nj <- ncol(stan.Nj.all)
  stan.Nj.all <- stan.Nj.all %>%
    dplyr::mutate(num_na = apply(., 1, count_na)) %>%
    dplyr::filter(!(num_na == ncol.Nj)) %>%
    dplyr::select(-num_na)

################################################################################
# SAVE
################################################################################

  saveRDS(svy.all, paste0(resdir, "all_svy_results_", today, ".rds"))
  saveRDS(lmer.all, paste0(resdir, "all_lmer_results_", today, ".rds"))
  saveRDS(stan.pars.all, paste0(resdir, "all_pars_results_", today, ".rds"))
  saveRDS(stan.ybar.all, paste0(resdir, "all_ybar_results_", today, ".rds"))
  saveRDS(stan.Nj.all, paste0(resdir, "all_Nj_results_", today, ".rds"))


