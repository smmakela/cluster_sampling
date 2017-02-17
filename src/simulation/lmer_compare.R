lmer_compare <- function(num.clusters, num.units, use.sizes, outcome.type, rootdir, sim) {
  # num.clusters -- number of clusters to sample
  # num.units -- number of units to sample
  # use.sizes -- 0/1 for whether cluster sizes used in pop data
  # outcome.type -- whether outcomes are continuous or binary
  # rootdir -- root directory for cluster_sampling
  # sim -- current iteration; used so that multiple instances aren't trying to write to the same file

  ##########################################
  ### Setup of libraries, load data
  ##########################################
    if (num.units <= 1) {
      nunits <- paste(num.units*100, "pct", sep = "")
    } else {
      nunits <- num.units
    }
    simdata <- readRDS(paste0(rootdir, "output/simulation/simdata_usesizes_",
                              use.sizes, "_", outcome.type, "_nclusters_",
                              num.clusters, "_nunits_", nunits, "_simno_", sim,
                              ".rds"))
print(str(simdata))
print(names(simdata))
    for (j in names(simdata)) {
      assign(j, simdata[[j]])
    }
    rm(simdata)
    popdata <- readRDS(paste0(rootdir, "/output/simulation/popdata_usesizes_",
                       use.sizes, "_", outcome.type, ".rds"))
    truepars <- popdata[["truepars"]]
    truepars.long <- gather(popdata[["truepars"]],
                            key = param.name, value = truth)
    print(truepars)
    rm(popdata)

  ##########################################
  ### Run lmer
  ##########################################
    #sample.data$xi_c <- sample.data$xi - mean(sample.data$xi)
    #sample.data$logMj_c <- log(sample.data$Mj) - mean(log(sample.data$Mj))
    # model 1: uses cluster sizes
    lmer_with_cluster_sizes <- lmer(y ~ x + logMj_c + x:logMj_c + (1 + x | cluster.id), data = sample.data)
    print(summary(lmer_with_cluster_sizes))
    # model 2: only uses cluster indicators
    lmer_cluster_indicators_only <- lmer(y ~ x + (1 + x | cluster.id), data = sample.data)
    print(summary(lmer_cluster_indicators_only))
    modlist <- c("lmer_with_cluster_sizes", "lmer_cluster_indicators_only")

  ##########################################
  ### Compare results
  ##########################################
    allres <- data.frame()
    for (j in 1:2) {
      assign("currmod", get(modlist[j]))
      lmer_summ <- summary(currmod)

      # pull out coefficient estimates from current model
      coefmat <- data.frame(lmer_summ$coefficients)
      if (j == 1) {
        coefmat$param.name <- c("gamma0", "alpha0", "gamma1", "alpha1")
      } else {
        coefmat$param.name <- c("gamma0", "alpha0")
      }
      coefmat <- left_join(coefmat, truepars.long, by = "param.name")
print(str(coefmat))
      coefmat <- coefmat[, c("param.name", "truth", "Estimate", "Std..Error")]
      # make names correspond to what's in the stan parameter estimate results
      names(coefmat) <- c("param.name", "truth", "mean", "sd")

      # pull out variance estimates from current model
      sdmat <- as.data.frame(VarCorr(currmod))[c(1,2,4), ]
      sdmat$param.name <- c("sigma_beta0", "sigma_beta1", "sigma_y")
      sdmat <- left_join(sdmat, truepars.long, by = "param.name")
      sdmat <- sdmat[, c("param.name", "sdcor", "truth")]
      names(sdmat) <- c("param.name", "mean", "truth")
      sdmat$sd <- NA
      both <- rbind(coefmat, sdmat)
      both$model.name <- modlist[j]
      both[["25%"]] <- both[["mean"]] - both[["sd"]]
      both[["75%"]] <- both[["mean"]] + both[["sd"]]
      both[["2.5%"]] <- both[["mean"]] - 1.96*both[["sd"]]
      both[["97.5%"]] <- both[["mean"]] + 1.96*both[["sd"]]
      allres <- rbind(allres, both)
    }

  ##########################################
  ### Print and save
  ##########################################
    print("######################################################################")
    print("compare to lmer results:")
    print(allres)
    print("######################################################################")
    
    #write.table(allres, file = paste(rootdir, "/Results/Simplify/vary_K/parests_lmer_usesizes_", use.sizes,
    #                                 "_nclusters_", num.clusters, "_nunits_", nunits,
    #                                 "_sim_", sim, ".txt", sep = ""), sep = ",")
    
    return(allres)                               

}
