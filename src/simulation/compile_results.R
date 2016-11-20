# Author: Susanna Makela
# Date: 13 Jan 2016
# Purpose: process results from simulation

##########################################
### Setup of directories and libraries
##########################################
  libdir <- "/vega/stats/users/smm2253/rpackages"
  .libPaths(libdir)
  rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
  Sys.setenv(HOME = rootdir)

  # set working directory, figure directory
  resdir <- paste0(rootdir, "output/simulation/", sep = "")
  setwd(resdir)

  # libraries
  library(plyr)

  # print time
  print(Sys.time())

  # list of unit and cluster numbers used
    num.clusters.list <- c(5, 10, 20, 50)
    #num.units.list <- c(0.05, 0.1)
    num.units.list <- c(0.05, 0.1, 0.25, 0.5, 1, 10, 50, 100)
    num.units.list2 <- c("5pct", "10pct", "25pct", "50pct", "100pct", 10, 50, 100)
    u.list <- c(0, 1)
    outcome.types <- c("continuous")

  # calculate true Mj total and true ybar
    true.values <- data.frame()
    for (usz in u.list) {
      for (ot in outcome.types) {
        popdata <- readRDS(paste0(resdir, "/popdata_usesizes_",
                                  usz, "_", ot, ".rds"))
        truepars.long <- gather(popdata[["truepars"]], key = param.name)
        tmp <- data.frame(use.sizes = usz,
                          outcome.type = ot,
                          Mj = sum(popdata[["Mj"]]),
                          ybar_true = mean(popdata[["pop.data"]]$yi),
                          truepars.long)
        true.values <- rbind(true.values, tmp)
      }
    }

  # which model to pull results for
    modelname.list <- c("knowsizes", "cluster_inds_only", "negbin", "bb")
    #modelname.list <- c("negbin", "bb")

##########################################
### Loop through results files
##########################################
  fil.list <- list.files(resdir, "results_.*_continuous_.*.rds")
  ybar.ests <- data.frame()
  param.ests <- data.frame()
  for (i in 1:length(fil.list)) {
    curr.name <- fil.list[i]
  
    # parse the filename to extract number of clusters, units, etc
    name.parts <- unlist(strsplit(curr.name, split = "_"))
    usz <- as.numeric(name.parts[3])
    out.tp <- name.parts[4]
    num.clusters <- as.numeric(name.parts[6])
    num.units <- name.parts[8]
    simno <- name.parts[10]
 
    # read in current file
    curr.fil <- readRDS(curr.name)

    # separate results on ybar and parameter estimates
    param.inds <- which(substr(names(curr.fil), 1, 5) == "param")
    param.res <- curr.fil[param.inds]
    param.res.names <- names(param.res)
    param.res.models <- gsub("param_ests_", "", param.res.names)
    for (j in 1:length(param.res.models)) {
      tmp <- param.res[[j]]
      if (param.res.models[j] != "lmer") {
        tmp$model.name <- param.res.models[j]
        tmp <- left_join(tmp, truepars.long, by = "param.name")
      } else { # add missing columns to lmer output so rbind() doesn't complain
        tmp[["se_mean"]] <- NA
        tmp[["50%"]] <- NA
        tmp[["n_eff"]] <- NA
        tmp[["Rhat"]] <- NA
      } # end if
      param.ests <- rbind(param.ests, tmp)
    } # end param.res.models loop
    param.ests$use.sizes <- usz
    param.ests$outcome.type <- out.tp
    param.ests$num.clusters <- num.clusters
    param.ests$num.units <- num.units
    param.ests$simno <- simno

    ybar.inds <- which(substr(names(curr.fil), 1, 4) == "ybar")
    ybar.res <- curr.fil[ybar.inds]
    ybar.res.names <- names(ybar.res)
    ybar.res.models <- gsub("ybar_ests_", "", ybar.res.names)
    ybar.ests <- data.frame()
    for (j in 1:length(ybar.res.models)) {
      tmp <- ybar.res[[j]]
      if (ybar.res.models[j] != "svy") {
        tmp$model.name <- ybar.res.models[j]
      }
      ybar.ests <- rbind(ybar.ests, tmp)
    } # end ybar.res.models loop
  } # end file for loop

print("str(param.ests):")
print(str(param.ests))
print("str(ybar.ests)")
print(str(ybar.ests))
##########################################
### Meta
##########################################
#  parests_all$metadata <- "J = 300, Nj ~ (100, 1000)"

##########################################
### save
##########################################
  today <- Sys.Date()
  today <- gsub("-", "_", today)
  res <- list(param.ests = param.ests, ybar.ests = ybar.ests)
  saveRDS(res, file = paste0(resdir, "/simulation_results_", today, ".rds"))


