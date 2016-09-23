lmer_compare <- function(num.clusters, num.units, rootdir, sim, popseed, use.sizes) {
  # num.clusters -- number of clusters to sample
  # num.units -- number of units to sample
  # rootdir -- root directory where Code, Data folders are
  # sim -- current iteration; used so that multiple instances aren't trying to write to the same file
  # popseed -- which seed was used to make the pop data

  ##########################################
  ### Setup of libraries, load data
  ##########################################
    library(lme4)
    if (num.units <= 1) {
      nunits <- paste(num.units*100, "pct", sep = "")
    } else {
      nunits <- num.units
    }
    load(paste(rootdir, "/Data/Simplify/vary_K/sampledata_usesizes_", use.sizes, "_nclusters_", num.clusters,
               "_nunits_", nunits, "_sim_", sim, ".RData", sep = ""))
    load(file = paste(rootdir, "/Data/Simplify/vary_K/popdata_usesizes_", use.sizes, "_seed_", popseed, ".RData", sep = ""))
    truepars <- popdata[["truepars"]]
    print(truepars)
    rm(popdata)

  ##########################################
  ### Run lmer
  ##########################################
    #sample.data$xi_c <- sample.data$xi - mean(sample.data$xi)
    #sample.data$logMj_c <- log(sample.data$Mj) - mean(log(sample.data$Mj))
    # model 1: uses cluster sizes
    lmer_mod1 <- lmer(yi ~ xi + logMj_c + xi:logMj_c + (1 + xi | cluster.id), data = sample.data)
    print(summary(lmer_mod1))
    # model 2: only uses cluster indicators
    lmer_mod2 <- lmer(yi ~ xi + (1 + xi | cluster.id), data = sample.data)
    print(summary(lmer_mod2))
    modlist <- c("lmer_mod1", "lmer_mod2")

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
    
    write.table(allres, file = paste(rootdir, "/Results/Simplify/vary_K/parests_lmer_usesizes_", use.sizes,
                                     "_nclusters_", num.clusters, "_nunits_", nunits,
                                     "_sim_", sim, ".txt", sep = ""), sep = ",")
                                   

}
