lmer_compare <- function(num_clusters, num_units, use_sizes, outcome_type,
                         size_model, sim_data, rootdir, simno) {
  # num_clusters -- number of clusters to sample
  # num_units -- number of units to sample
  # use_sizes -- whether y depends on cluster sizes or not
  # outcome_type -- whether y is continuous or binary
  # size_model -- model used to generate pop cluster sizes
  # sim_data -- list of data for simulation - pop data, sampled data, etc
  # rootdir -- root directory where Code, Data folders are
  # simno -- current iteration; used so that multiple instances aren't trying to write to the same file

  ##########################################
  ### Setup of libraries, load data
  ##########################################
    if (num_units <= 1) {
      nunits <- paste(num_units*100, "pct", sep = "")
    } else {
      nunits <- num_units
    }
    for (j in names(sim_data)) {
      assign(j, sim_data[[j]])
    }
    popdata <- readRDS(paste0(rootdir, "/output/simulation/popdata_usesizes_",
                              use_sizes, "_", outcome_type, "_", size_model, ".rds"))
    truepars <- popdata[["truepars"]]
    print(truepars)
    rm(popdata)

  ##########################################
  ### Run lmer/glmer
  ##########################################
  if (outcome_type == "continuous") {
    # model 1: uses cluster sizes
    lmer_with_cluster_sizes <- lmer(y ~ x + logMj_c + x:logMj_c +
                                    (1 + x | cluster_id), data = sample_data)
    print(summary(lmer_with_cluster_sizes))
    # model 2: only uses cluster indicators
    lmer_cluster_indicators_only <- lmer(y ~ x + (1 + x | cluster_id),
                                         data = sample_data)
    print(summary(lmer_cluster_indicators_only))
    modlist <- c("lmer_with_cluster_sizes", "lmer_cluster_indicators_only")
  } else {
    # model 1: uses cluster sizes
    glmer_with_cluster_sizes <- glmer(y ~ logMj_c + (1 | cluster_id),
                                      family = binomial(link = "logit"),
                                      data = sample_data)
    print(summary(glmer_with_cluster_sizes))
    # model 2: only uses cluster indicators
    glmer_cluster_indicators_only <- glmer(y ~ (1 | cluster_id),
                                           family = binomial(link = "logit"),
                                           data = sample_data)
    print(summary(glmer_cluster_indicators_only))
    modlist <- c("glmer_with_cluster_sizes", "glmer_cluster_indicators_only")
  }

  ##########################################
  ### Compare results
  ##########################################
    allres <- data.frame()
    for (j in 1:2) {
      assign("currmod", get(modlist[j]))
      lmer_summ <- summary(currmod)
      coefmat <- data.frame(lmer_summ$coefficients)
      if (outcome_type == "continuous") {
        if (j == 1) {
          rownames(coefmat) <- c("alpha0", "alpha1", "gamma0", "gamma1")
        } else {
          rownames(coefmat) <- c("alpha0", "alpha1")
        }
print("---------------------------")
print(coefmat)
print(truepars)
print(as.data.frame(VarCorr(currmod)))
        truth <- truepars[,rownames(coefmat)]
        coefmat <- cbind(coefmat, truth = t(truth))
        sdmat <- as.data.frame(VarCorr(currmod))[c(1,2,4),]
        rownames(sdmat) <- c("sigma_beta0", "sigma_beta1", "sigma_y")
      } else {
        if (j == 1) {
          rownames(coefmat) <- c("alpha0", "gamma0")
        } else {
          rownames(coefmat) <- c("alpha0")
        }
print("---------------------------")
print(coefmat)
print(truepars)
print(as.data.frame(VarCorr(currmod)))
        truth <- truepars[,rownames(coefmat)]
        coefmat <- cbind(coefmat, truth = t(truth))
        sdmat <- as.data.frame(VarCorr(currmod))[1,]
        rownames(sdmat) <- c("sigma_beta0")
      }
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
    print("compare to (g)lmer results:")
    print(allres)
    print("######################################################################")
    
    saveRDS(allres,
            paste0(rootdir, "output/simulation/lmer_ests_usesizes_",
                   use_sizes, "_", outcome_type, "_", size_model, "_nclusters_",
                   num_clusters, "_nunits_", nunits, "_simno_", simno, ".rds"))
    return(NULL)
}

