# Author: Susanna Makela
# Date: 21 Apr 2014
# Purpose: generate population data

makedata <- function(seed = NA, numclusters, clustersize.range, unitcovar.range, rootdir, use.sizes) {
  # seed -- integer for setting seed so we can reproduce the results (optional)
  # numclusters -- number of clusters in population
  # clustersize.range -- min and max of possible clustersizes
  # unitcovar.range -- min and max of possible unit-level covariate
  # rootdir -- root directory where Code, Data folders are
  # use.sizes -- 0/1 for whether to use cluster sizes in pop model

  ##########################################
  ### Checks
  ##########################################
  if (length(clustersize.range) != 2) {
    stop("clustersize.range must be of length 2")
  }
  if (clustersize.range[1] >= clustersize.range[2]) {
    stop("clustersize.range must be of the form (min.clustersize, max.clustersize)")
  }
  if (length(unitcovar.range) != 2) {
    stop("unitcovar.range must be of length 2")
  }
  if (unitcovar.range[1] >= unitcovar.range[2]) {
    stop("unitcovar.range must be of the form (min.unitcovar, max.unitcovar)")
  }
  
  ##########################################
  ### Setup of libraries and directories
  ##########################################
    libdir <- "/vega/stats/users/smm2253/rpackages"
    .libPaths(libdir)
    library(plyr)
  
  ##########################################
  ### Generate data
  ##########################################
    # set the seed for the random number generator so we can reproduce the results
      if (!is.na(seed)) {
        set.seed(seed) 
      }
    
    # number of clusters, cluster size, unit-level covariate
      J <- numclusters
      Mj <- sample(c(clustersize.range[1]:clustersize.range[2]), J, replace = TRUE)
      xi <- sample(c(unitcovar.range[1]:unitcovar.range[2]), sum(Mj), replace = TRUE)
      xi <- xi - mean(xi)
 
    # load hyperparameters
      gamma0 <- rnorm(1)
      alpha0 <- rnorm(1)
      if (use.sizes == 1) {
        gamma1 <- rnorm(1)
        alpha1 <- rnorm(1)
      }  else {
        gamma1 <- 0
        alpha1 <- 0
      }
      sigma_beta0 <- abs(rnorm(1, 0, 0.5))
      sigma_beta1 <- abs(rnorm(1, 0, 0.5))
      sigma_y <- abs(rnorm(1, 0, 0.75))
      #sigma_beta0 <- abs(rcauchy(1, 0, 2.5))
      #sigma_beta1 <- abs(rcauchy(1, 0, 2.5))
      #sigma_y <- abs(rcauchy(1, 0, 2.5))
      
    # draw population parameters and yij's
      logMj_c <- log(Mj) - mean(log(Mj))
      #logMj_c <- log(Mj)
      beta0j <- rnorm(n = J, mean = gamma0 + gamma1*logMj_c, sd = sigma_beta0)
      b0j <- rep(beta0j, Mj)
      beta1j <- rnorm(n = J, mean = alpha0 + alpha1*logMj_c, sd = sigma_beta1)
      b1j <- rep(beta1j, Mj)
      ymean <- b0j + b1j*xi 
      yi <- rnorm(ymean, mean = ymean, sd = sigma_y)

    # make data frame of pop data
      pop.data <- data.frame(yi, xi, Mj = rep(Mj, Mj), logMj_c = rep(logMj_c, Mj))
      pop.data$cluster.id <- rep(c(1:J), times = Mj)
      pop.data$unit.id <- unlist(lapply(Mj, seq_len))
print("HEREEEEE") 
  # make list of things to save
    truepars <- data.frame(gamma0, gamma1, alpha0, alpha1, sigma_beta0, sigma_beta1, sigma_y)
    popdata <- list(pop.data = pop.data, J = J, Mj = Mj, logMj_c = logMj_c, beta0j = beta0j, beta1j = beta1j, truepars = truepars)
    if (is.na(seed)) {
      popseed <- 1
    } else {
      popseed <- seed
    }
    save(popdata, file = paste(rootdir, "/Data/Simplify/vary_K/popdata_usesizes_", use.sizes, "_seed_", popseed, ".RData", sep = ""))
}
