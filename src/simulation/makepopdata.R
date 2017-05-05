#/usr/bin/R Rscript
# Author: Susanna Makela
# Date: 21 Apr 2014
# Purpose: generate population data

'Usage: makepopdata.r [options]

  -s --seed <seedval>                      Integer for setting seed so we can reproduce the results [default: NULL]
  -n --numclusters <J>                     Number of clusters in population [default: 100]
  -u --unitcovar_range <x_range>           Min, max values for unit covariate [default: 20,45]
  -z --use_sizes <usz>                     Whether outcomes are related to cluster sizes [default: 0]
  -o --outcome_type <y_type>               Whether outcomes are continuous or binary [default: continuous]
  -r --rootdir <dirname>                   Root project directory [default: /vega/stats/users/smm2253/cluster_sampling]
' -> doc

  print("before")
  print(doc)
  print("after")
  
  #############################################################################
  ### Setup of libraries and directories
  #############################################################################
    libdir <- "/vega/stats/users/smm2253/rpackages"
    .libPaths(libdir)
    library(dplyr)
    require(docopt)
    require(methods)

    codedir <- "/vega/stats/users/smm2253/cluster_sampling/src/simulation"
    source(paste0(codedir, "/draw_pop_cluster_sizes_for_sim.R"))
  
    # Store the docopt options as variables we can use in the code
    opts <- docopt(doc) 
    opts.names <- names(opts)
    opts.names <- opts.names[-grep("--", opts.names)]
    print(str(opts))
    print(opts.names)
    for (j in 1:length(opts.names)) {
      curr.name <- opts.names[j]
      curr.name2 <- gsub("_", ".", curr.name)
      print(curr.name2)
      # All arguments get read in as characters, so convert the numeric ones to numeric
      # First check if the argument contains a comma. If it does, the contents
      #   are of the form "i, j", where i and j are numbers, and we need to
      #   convert them into a numeric vector
      # Then check if the argument can be coerced to numeric -- if not, it's a
      #   string and we don't need to do anything
      # If the argument *can* be converted to numeric, do that
      if (as.numeric(regexpr(",", opts[[curr.name]])) > 0) {
        assign(curr.name2, as.numeric(unlist(strsplit(opts[[curr.name]], ","))))
      } else if (is.na(as.numeric(opts[[curr.name]]))) {
        assign(curr.name2, opts[[curr.name]])
        if (curr.name2 == "seed" & opts[[curr.name]] == "NULL") {
          seed <- NULL
        }
      } else {
        assign(curr.name2, as.numeric(opts[[curr.name]]))
      }
      print(str(get(curr.name2)))
    }
    J <- numclusters # since J is used in the rest of the code

  #############################################################################
  ### Useful function 
  #############################################################################
    inv_logit <- function(x) {
      y <- exp(x) / (1 + exp(x))
      return(y)
    }

  #############################################################################
  ### Checks
  #############################################################################
    #if (length(clustersize.range) != 2) {
    #  stop("clustersize.range must be of length 2")
    #}
    #if (clustersize.range[1] >= clustersize.range[2]) {
    #  stop(paste0("clustersize.range must be of the form
    #               (min.clustersize, max.clustersize)"))
    #}
    if (length(unitcovar.range) != 2) {
      stop("unitcovar.range must be of length 2")
    }
    if (unitcovar.range[1] >= unitcovar.range[2]) {
      stop("unitcovar.range must be of the form (min.unitcovar, max.unitcovar)")
    }
  
  #############################################################################
  ### Generate data
  #############################################################################
    # Set the seed for the random number generator to reproduce the results
    if (!is.null(seed)) {
      set.seed(seed) 
    }
    
    # Draw number of clusters, cluster size, unit-level covariate
    # Here cluster sizes are taken from the observed pop cluster sizes by
    # fitting a gamma distribution to the observed sizes, drawing from the
    # fitted gamma, and then exponentiating; note we have to ROUND here so we
    # get integer sizes
    print("ENTERING DRAW_POP_CLUSTER_SIZES_FOR_SIM.R")
    Mj <- round(draw_pop_cluster_sizes_for_sim(J))
    print("str(Mj):")
    print(str(Mj))
    #Mj <- sample(c(clustersize.range[1]:clustersize.range[2]), J,
    #             replace = TRUE)
    logMj_c <- log(Mj) - mean(log(Mj))
    x <- base::sample(c(unitcovar.range[1]:unitcovar.range[2]), sum(Mj),
                      replace = TRUE)
    x <- x - mean(x)
print(length(Mj))
print(sum(Mj))
print(length(x))   
print(summary(x))
print(summary(Mj))

    # Draw hyperparameters, varying slopes and coefficients, and outcomes
    if (outcome.type == "continuous") {
      alpha0 <- rnorm(1)
      alpha1 <- rnorm(1)
      if (use.sizes == 1) {
        gamma0 <- rnorm(1)
        gamma1 <- rnorm(1)
      }  else {
        gamma0 <- 0
        gamma1 <- 0
      }
      sigma_beta0 <- abs(rnorm(1, 0, 0.5))
      sigma_beta1 <- abs(rnorm(1, 0, 0.5))
      sigma_y <- abs(rnorm(1, 0, 0.75))
        
      beta0 <- rnorm(n = J, mean = alpha0 + gamma0 * logMj_c, sd = sigma_beta0)
      beta0_rep <- rep(beta0, Mj)
      beta1 <- rnorm(n = J, mean = alpha1 + gamma1 * logMj_c, sd = sigma_beta1)
      beta1_rep <- rep(beta1, Mj)
      ymean <- beta0_rep + beta1_rep * x 
      y <- rnorm(ymean, mean = ymean, sd = sigma_y)
    } else {
      alpha0 <- rnorm(1)
      if (use.sizes == 1) {
        gamma0 <- rnorm(1)
      }  else {
        gamma0 <- 0
      }
      sigma_beta0 <- abs(rnorm(1, 0, 0.5))
      sigma_y <- abs(rnorm(1, 0, 0.75))
        
      beta0 <- rnorm(n = J, mean = alpha0 + gamma0*logMj_c, sd = sigma_beta0)
      beta0_rep <- rep(beta0, Mj)
      y_prob <- inv_logit(beta0_rep) 
      y <- rbinom(y_prob, size = 1, prob = y_prob)
      # make x = 0 everywhere so we can still use svy_ests.R
      x_new <- rep(0, length(x))
      x <- x_new
    }

    # Make data frame of pop data
    pop.data <- data.frame(y, x, Mj = rep(Mj, Mj), logMj_c = rep(logMj_c, Mj))
    pop.data$cluster.id <- rep(c(1:J), times = Mj)
    pop.data$unit.id <- unlist(lapply(Mj, seq_len))

    # Make list of things to save
    ybar_true <- mean(pop.data$y)
    if (outcome.type == "continuous") {
      truepars <- data.frame(alpha0, gamma0, alpha1, gamma1, sigma_beta0,
                             sigma_beta1, sigma_y, ybar_true)
    } else {
      truepars <- data.frame(alpha0, gamma0,sigma_beta0,
                             sigma_y, ybar_true)
      beta1 = NA
    }
    popdata <- list(pop.data = pop.data, J = J, Mj = Mj, logMj_c = logMj_c,
                    beta0 = beta0, beta1 = beta1, truepars = truepars)
    print(str(popdata))
    saveRDS(popdata,
            file = paste0(rootdir, "/output/simulation/popdata_usesizes_",
                          use.sizes, "_", outcome.type, ".rds"))

