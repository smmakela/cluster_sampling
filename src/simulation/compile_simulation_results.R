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
  library(dplyr)
  library(tidyr)

  # print time
  print(Sys.time())

  # calculate true Mj total and true ybar
    true.values <- data.frame()
    for (usz in c(0, 1)) {
      for (ot in c("continuous")) {
        popdata <- readRDS(paste0(resdir, "/popdata_usesizes_",
                                  usz, "_", ot, ".rds"))
        ybar.true <- mean(popdata[["pop.data"]]$y)
        truepars.long <- gather(popdata[["truepars"]], key = param.name)
        names(truepars.long) <- c("param.name", "truth")
        truepars.long <- rbind(truepars.long, c("ybar_new", ybar.true))
        tmp <- data.frame(use.sizes = usz,
                          outcome.type = ot,
                          M_tot = sum(popdata[["Mj"]]),
                          truepars.long)
        true.values <- rbind(true.values, tmp)
      }
    }
true.values$param.name <- as.character(true.values$param.name)
print(true.values)  

##########################################
### Loop through results files
##########################################
  fil.list <- list.files(resdir, "results_.*_continuous_.*.rds")
  ybar.ests <- data.frame()
  param.ests <- data.frame()
  for (i in 1:length(fil.list)) {
    # read in current file
    curr.name <- fil.list[i]
    curr.fil <- readRDS(curr.name)
    print(paste0("Currently on: ", i, " of ", length(fil.list), " files." ))

    # parse the filename to extract number of clusters, units, etc
    curr.name <- gsub(".rds", "", curr.name) 
    name.parts <- unlist(strsplit(curr.name, split = "_"))
    usz <- as.numeric(name.parts[3])
    out.tp <- name.parts[4]
    num.clusters <- as.numeric(name.parts[6])
    num.units <- name.parts[8]
    simno <- as.numeric(name.parts[10])
 
    # separate results on ybar and parameter estimates
    param.inds <- which(substr(names(curr.fil), 1, 5) == "param")
    param.res <- curr.fil[param.inds]
    param.res.names <- names(param.res)
    param.res.models <- gsub("param_ests_", "", param.res.names)
    for (j in 1:length(param.res.models)) {
      tmp <- param.res[[j]]
      if (grepl("lmer", param.res.models[j]) == FALSE) {
        tmp$model.name <- param.res.models[j]
      } else { # add missing columns to lmer output so rbind() doesn't complain
        tmp[["se_mean"]] <- NA
        tmp[["50%"]] <- NA
        tmp[["n_eff"]] <- NA
        tmp[["Rhat"]] <- NA
        tmp$truth <- NULL
      } # end if
      true.values.tmp <- dplyr::filter(true.values, use.sizes == usz)
      tmp <- left_join(tmp, true.values.tmp, by = "param.name")
      tmp$use.sizes <- usz
      tmp$outcome.type <- out.tp
      tmp$num.clusters <- num.clusters
      tmp$num.units <- num.units
      tmp$simno <- simno
      param.ests <- rbind(param.ests, tmp)
    } # end param.res.models loop

    ybar.inds <- which(substr(names(curr.fil), 1, 4) == "ybar")
    ybar.res <- curr.fil[ybar.inds]
    ybar.res.names <- names(ybar.res)
    ybar.res.models <- gsub("ybar_ests_", "", ybar.res.names)
    for (j in 1:length(ybar.res.models)) {
      tmp <- ybar.res[[j]]
      if (grepl("svy", ybar.res.models[j]) == FALSE) {
        tmp$model.name <- ybar.res.models[j]
      }
      tmp$use.sizes <- usz
      tmp$outcome.type <- out.tp
      tmp$num.clusters <- num.clusters
      tmp$num.units <- num.units
      tmp$simno <- simno
      ybar.ests <- rbind(ybar.ests, tmp)
    } # end ybar.res.models loop
  } # end file for loop
print(warnings())
tmp <- dplyr::summarise(group_by(param.ests, model.name, use.sizes, num.clusters,
                                 num.units),
                        num.obs = n())
print("str(param.ests):")
print(str(param.ests))
print(tmp)
                       
tmp <- dplyr::summarise(group_by(ybar.ests, use.sizes, num.clusters, num.units),
                        num.obs = n())
print("str(ybar.ests)")
print(str(ybar.ests))
print(tmp)
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
  saveRDS(res, file = paste0(resdir, "/compiled_simulation_results.rds"))
  saveRDS(res, file = paste0(resdir, "/compiled_simulation_results_",
                             today, ".rds"))
