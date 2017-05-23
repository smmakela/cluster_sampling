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
# Calculate true values
################################################################################
  cat("#####################################################################\n")
  cat("#####################################################################\n")
  cat("Compiling true values\n")
  params.true <- data.frame()
  Mj.true.summ <- data.frame()
  ybar.true <- data.frame()
  for (usz in c(0, 1)) {
    for (ot in c("continuous")) {
      popdata <- readRDS(paste0(resdir, "/popdata_usesizes_",
                                usz, "_", ot, ".rds"))
      Mj <- popdata[["Mj"]]
      truepars.long <- gather(popdata[["truepars"]], key = param.name)
      names(truepars.long) <- c("param.name", "true")
      tmp <- data.frame(use.sizes = usz,
                        outcome.type = ot,
                        truepars.long)
      params.true <- rbind(params.true, tmp)
      tmp2 <- data.frame(use.sizes = usz,
                         outcome.type = ot,
                         sum  = sum(Mj),
                         mean = mean(Mj),
                         sd   = sd(Mj),
                         p025 = quantile(Mj, 0.025),
                         p25  = quantile(Mj, 0.25),
                         p50  = quantile(Mj, 0.50),
                         p75  = quantile(Mj, 0.75),
                         p975 = quantile(Mj, 0.975))
      Mj.true.summ <- rbind(Mj.true.summ, tmp2)
      tmp3 <- data.frame(use.sizes = usz,
                         outcome.type = ot,
                         ybar_true = mean(popdata[["pop.data"]]$y))
      ybar.true <- rbind(ybar.true, tmp3)
    }
  }
  params.true$param.name <- as.character(params.true$param.name)
  print("params.true:") 
  print(params.true) 
  Mj.true.summ <- Mj.true.summ %>%
    tidyr::gather(key = stat, value = true, -c(use.sizes, outcome.type))
  print("Mj.true.summ:") 
  print(Mj.true.summ) 
  print("ybar.true:") 
  print(ybar.true) 
  print(Sys.time())

################################################################################
### Loop through STAN result files
################################################################################
  cat("#####################################################################\n")
  cat("#####################################################################\n")
  cat(" Now doing STAN files\n")
  fil.list <- list.files(resdir, "stan_results_.*.rds")
  ybar.ests <- data.frame()
  param.ests <- data.frame()
  Nj.ests <- data.frame()
  div.trans <- expand.grid(use.sizes = c(0, 1),
                           outcome.type = c("binary", "continuous"),
                           model.name = c("bb", "cluster_inds_only",
                                          "knowsizes", "lognormal", "negbin"),
                           num.clusters = c(5, 10, 20, 30),
                           num.units = c(0.05, 0.1, 0.25, 0.5, 1, 10, 30, 60))
  for (i in 1:length(fil.list)) {
    # read in current file
    curr.name <- fil.list[i]
    curr.fil <- readRDS(curr.name)
     if ((i %% 500) == 0) {
       cat("Currently on: ", i, " of ", length(fil.list), " stan files.\n")
       cat(curr.name, "\n")
    }
    # parse the filename to extract number of clusters, units, etc
    curr.name <- gsub(".rds", "", curr.name) 
    curr.name <- gsub("cluster_inds_only", "cluster.inds.only", curr.name) 
    name.parts <- unlist(strsplit(curr.name, split = "_"))
    usz <- as.numeric(name.parts[4])
    out.tp <- name.parts[5]
    model.name <- name.parts[6]
    num.clusters <- as.numeric(name.parts[8])
    num.units <- name.parts[10]
    simno <- as.numeric(name.parts[12])

    # see if this was a case of divergent transitions we couldn't get rid of
    if (length(curr.fil) == 1) {
      ind <- div.trans$use.sizes    == usz &
             div.trans$outcome.type == out.tp &
             div.trans$model.name   == model.name &
             div.trans$num.clusters == num.clusters &
             div.trans$num.units    == num.units
      div.trans$counter[ind] <- div.trans$counter[ind] + 1
      next
    }

    # pull out parameter estimates, make so that statistics are long
    tmp <- curr.fil[["par_ests"]]
    tmp$param.name <- rownames(tmp)
    tmp <- tmp %>%
      tidyr::gather(key = stat, value = value, -param.name)
    params.true.tmp <- dplyr::filter(params.true,
                                     use.sizes == usz & outcome.type == out.tp)
    tmp <- left_join(tmp, params.true.tmp, by = "param.name")
    tmp$use.sizes <- usz
    tmp$outcome.type <- out.tp
    tmp$model.name <- model.name
    tmp$num.clusters <- num.clusters
    tmp$num.units <- num.units
    tmp$simno <- simno
    param.ests <- bind_rows(param.ests, tmp)

    # pull out summaries of ybar_new
    draw.summ <- curr.fil[["draw_summ"]]
    if (nrow(draw.summ) == 1) {
      tmp <- draw.summ
    } else {
      tmp <- draw.summ[draw.summ$param == "ybar_new", ]
    }
    tmp$use.sizes <- usz
    tmp$outcome.type <- out.tp
    tmp$model.name <- model.name
    tmp$num.clusters <- num.clusters
    tmp$num.units <- num.units
    tmp$simno <- simno
    tmp$ybar_true <- ybar.true$ybar_true[ybar.true$use.sizes == usz &
                                         ybar.true$outcome.type == out.tp]
    ybar.ests <- bind_rows(ybar.ests, tmp)

    # pull out summaries of Nj_new
    if (model.name %in% c("bb", "lognormal", "negbin")) {
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
      tmp$use.sizes <- usz
      tmp$outcome.type <- out.tp
      tmp$model.name <- model.name
      tmp$num.clusters <- num.clusters
      tmp$num.units <- num.units
      tmp$simno <- simno
      tmp <- left_join(tmp, Mj.true.summ, by = c("use.sizes", "outcome.type"))
      Nj.ests <- bind_rows(Nj.ests, tmp)
    }

  } # end file for loop
  print(warnings())
  print(Sys.time())

  res <- list(param.ests = param.ests, ybar.ests = ybar.ests,
              Nj.ests = Nj.ests, div.trans = div.trans)
  saveRDS(res, file = paste0(resdir, "/compiled_stan_results_", today, ".rds"))
  rm(res)

  # summarise how many sims we got for everything
  param.ests.nobs <- param.ests %>%
    dplyr::(param.name, model.name, use.sizes, outcome.type,
            num.clusters, num.units) %>%
    dplyr::summarise(num.obs = n())
  print("str(param.ests):")
  print(str(param.ests))
  print(param.ests.nobs)
  print("num obs:")
  print(param.ests.nobs$num.obs)
  print("div trans:")
  print(div.trans$count)
                       
################################################################################
### Loop through LMER result files
################################################################################
  cat("#####################################################################\n")
  cat("#####################################################################\n")
  cat(" Now doing LMER files\n")
  fil.list <- list.files(resdir, "lmer_ests_.*.rds")
  param.ests.lmer <- data.frame()
  for (i in 1:length(fil.list)) {
    # read in current file
    curr.name <- fil.list[i]
    curr.fil <- readRDS(curr.name)
    if ((i %% 500) == 0) {
      print(paste0("Currently on: ", i, " of ", length(fil.list), " lmer files." ))
      cat(curr.name, "\n")
    }

    # parse the filename to extract number of clusters, units, etc
    curr.name <- gsub(".rds", "", curr.name) 
    curr.name <- gsub("cluster_inds_only", "cluster.inds.only", curr.name) 
    name.parts <- unlist(strsplit(curr.name, split = "_"))
    usz <- as.numeric(name.parts[4])
    out.tp <- name.parts[5]
    num.clusters <- as.numeric(name.parts[7])
    num.units <- name.parts[9]
    simno <- as.numeric(name.parts[11])

    # pull out parameter estimates
    curr.fil$whichmodel[curr.fil$whichmodel == 1] <- "lmer_with_sizes"
    curr.fil$whichmodel[curr.fil$whichmodel == 2] <- "lmer_no_sizes"
    names(curr.fil) <- c("mean", "sd", "true", "model.name")
    curr.fil$param.name <- rownames(curr.fil)
    curr.fil$param.name <- gsub("01", "0", curr.fil$param.name)
    curr.fil$param.name <- gsub("11", "1", curr.fil$param.name)
    params.true.tmp <- dplyr::filter(params.true,
                                     use.sizes == usz & outcome.type == out.tp)
    curr.fil <- left_join(curr.fil, params.true.tmp, by = "param.name")
    curr.fil$use.sizes <- usz
    curr.fil$outcome.type <- out.tp
    curr.fil$num.clusters <- num.clusters
    curr.fil$num.units <- num.units
    curr.fil$simno <- simno
    param.ests.lmer <- rbind(param.ests.lmer, curr.fil)
  }
  print(warnings())
  print(Sys.time())

  res <- list(param.ests.lmer = param.ests.lmer)
  saveRDS(res, file = paste0(resdir, "/compiled_lmer_results_", today, ".rds"))
  rm(res)
  
  # summarise how many sims we got for everything
  param.ests.nobs <- param.ests.lmer %>%
    dplyr::filter(param.name == "sigma_beta0" & model.name == "lmer_with_sizes") %>%
    dplyr::group_by(use.sizes, outcome.type, num.clusters, num.units) %>%
    dplyr::summarise(num.obs = n())
  print("str(param.ests.lmer):")
  print(str(param.ests.lmer))
  print(param.ests.nobs)
  print("num obs")
  print(param.ests.nobs$num.obs)

################################################################################
### Loop through SVY result files
################################################################################
  cat("#####################################################################\n")
  cat("#####################################################################\n")
  cat(" Now doing SVY files\n")
  fil.list <- list.files(resdir, "svy_ests_.*.rds")
  ybar.ests.svy <- data.frame()
  for (i in 1:length(fil.list)) {
    # read in current file
    curr.name <- fil.list[i]
    curr.fil <- readRDS(curr.name)
    if ((i %% 500) == 0) {
      print(paste0("Currently on: ", i, " of ", length(fil.list), " svy files." ))
      cat(curr.name, "\n")
    }

    # parse the filename to extract number of clusters, units, etc
    curr.name <- gsub(".rds", "", curr.name) 
    curr.name <- gsub("cluster_inds_only", "cluster.inds.only", curr.name) 
    name.parts <- unlist(strsplit(curr.name, split = "_"))
    usz <- as.numeric(name.parts[4])
    out.tp <- name.parts[5]
    num.clusters <- as.numeric(name.parts[7])
    num.units <- name.parts[9]
    simno <- as.numeric(name.parts[11])

    # pull out parameter estimates
    curr.fil$use.sizes <- usz
    curr.fil$outcome.type <- out.tp
    curr.fil$num.clusters <- num.clusters
    curr.fil$num.units <- num.units
    curr.fil$simno <- simno
    ybar.ests.svy <- rbind(ybar.ests.svy, curr.fil)
  }
  print(warnings())
  print(Sys.time())

  res <- list(ybar.ests.svy = ybar.ests.svy)
  saveRDS(res, file = paste0(resdir, "/compiled_svy_results_", today, ".rds"))
  rm(res)
  
  # summarise how many sims we got for everything
  ybar.ests.nobs <- ybar.ests.svy %>%
    dplyr::filter(model.name == "hajek") %>%
    dplyr::group_by(use.sizes, outcome.type, num.clusters, num.units) %>%
    dplyr::summarise(num.obs = n())
  print("str(ybar.ests.svy):")
  print(str(ybar.ests.svy))
  print(ybar.ests.nobs)
  print("num obs")
  print(ybar.ests.nobs$num.obs)


