lmer_compare <- function(num.clusters, num.units, use.sizes, outcome.type,
                         rootdir, sim) {
  # num.clusters -- number of clusters to sample
  # num.units -- number of units to sample
  # use.sizes -- whether y depends on cluster sizes or not
  # rootdir -- root directory where Code, Data folders are
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
    for (j in names(simdata)) {
      assign(j, simdata[[j]])
    }
    rm(simdata)
    popdata <- readRDS(paste0(rootdir, "/output/simulation/popdata_usesizes_",
                       use.sizes, "_", outcome.type, ".rds"))
    truepars <- popdata[["truepars"]]
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
      coefmat <- data.frame(lmer_summ$coefficients)
      if (j == 1) {
        rownames(coefmat) <- c("gamma0", "alpha0", "gamma1", "alpha1")
      } else {
        rownames(coefmat) <- c("gamma0", "alpha0")
      }
      truth <- truepars[,rownames(coefmat)]
      coefmat <- cbind(coefmat, truth = t(truth))
      sdmat <- as.data.frame(VarCorr(currmod))[c(1,2,4),]
      rownames(sdmat) <- c("sigma_beta0", "sigma_beta1", "sigma_y")
      truth <- truepars[,rownames(sdmat)]
      sdmat <- cbind(sdmat, truth = t(truth))
      sdmat <- sdmat[, c("sdcor", "truth")]
      names(sdmat) <- c("est", "truth")
      coefmat <- coefmat[, c("Estimate", "Std..Error", "truth")]
      names(coefmat) <- c("est", "stderr", "truth")
      sdmat$stderr <- NA
      sdmat <- sdmat[, c("est", "stderr", "truth")]
      both <- rbind(coefmat, sdmat)
      both$whichmodel <- j
      allres <- rbind(allres, both)
    }

  ##########################################
  ### Print and save
  ##########################################
    print("######################################################################")
    print("compare to lmer results:")
    print(allres)
    print("######################################################################")
    
    return(allres)
    #write.table(allres, file = paste(rootdir, "/Results/Simplify/vary_K/parests_lmer_usesizes_", use.sizes,
    #                                 "_nclusters_", num.clusters, "_nunits_", nunits,
    #                                 "_sim_", sim, ".txt", sep = ""), sep = ",")
                                   

}
