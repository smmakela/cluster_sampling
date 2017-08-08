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
  -m --size_model <sz_model>               Model for how the population cluster szies are generated [default: multinomial]
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

    # Store the docopt options as variables we can use in the code
    opts <- docopt(doc) 
    opts_names <- names(opts)
    opts_names <- opts_names[-grep("--", opts_names)]
    print(str(opts))
    print(opts_names)
    for (j in 1:length(opts_names)) {
      curr_name <- opts_names[j]
      cat("Currently on:", curr_name, "\n")
      # All arguments get read in as characters, so convert the numeric ones to numeric
      # First check if the argument contains a comma. If it does, the contents
      #   are of the form "i, j", where i and j are numbers, and we need to
      #   convert them into a numeric vector
      # Then check if the argument can be coerced to numeric -- if not, it's a
      #   string and we don't need to do anything
      # If the argument *can* be converted to numeric, do that
      if (as.numeric(regexpr(",", opts[[curr_name]])) > 0) {
        assign(curr_name, as.numeric(unlist(strsplit(opts[[curr_name]], ","))))
      } else if (is.na(as.numeric(opts[[curr_name]]))) {
        assign(curr_name, opts[[curr_name]])
        if (curr_name == "seed" & opts[[curr_name]] == "NULL") {
          seed <- NULL
        }
      } else {
        assign(curr_name, as.numeric(opts[[curr_name]]))
      }
      print(str(get(curr_name)))
    }
    if (grepl("ff", size_model)) {
      J <- 77 # number of cities in FF data
    } else {
      J <- numclusters
    }

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
    #if (length(clustersize_range) != 2) {
    #  stop("clustersize_range must be of length 2")
    #}
    #if (clustersize_range[1] >= clustersize_range[2]) {
    #  stop(paste0("clustersize_range must be of the form
    #               (min_clustersize, max_clustersize)"))
    #}
    if (length(unitcovar_range) != 2) {
      stop("unitcovar_range must be of length 2")
    }
    if (unitcovar_range[1] >= unitcovar_range[2]) {
      stop("unitcovar_range must be of the form (min_unitcovar, max_unitcovar)")
    }
  
  #############################################################################
  ### Generate data -- pop cluster sizes are generated from:
  ### 1. FF observed city sizes
  ### 2. 500 * multinomial draw
  ### 3. Poisson
  #############################################################################
    # Set the seed for the random number generator to reproduce the results
    if (!is.null(seed)) {
      set.seed(seed) 
    }

    # Draw cluster sizes using the specified size_model
    if (size_model == "multinomial") {
      alpha_vec <- rep(10, times = J) # alpha for dirichlet distribution
      # take dirichlet draws by scaling gamma draws
      gamma_draws <- sort(rgamma(alpha_vec, shape = alpha_vec, scale = 1))
      # make sure we don't get any cluster sizes smaller than 325 (max number of
      # births we are gonna sample)
      while (min(round(100*gamma_draws)) < 325) {
        gamma_draws <- sort(rgamma(alpha_vec, shape = alpha_vec, scale = 1))
      }
      dir_draws <- gamma_draws / sum(gamma_draws)
      # use multinomial to draw which gamma sizes we'll use
      mn_draws <- rmultinom(n = 1, size = J, prob = dir_draws)
      # round gamma draws and add 1 (to avoid zero), pick the ones specified
      # by mn_draws
      Mj <- as.numeric(round(rep(100*gamma_draws, times = mn_draws)))
      stratum_id <- rep(1, length(Mj)) # filler
      # Collect all into data frame
      pop_dat <- data.frame(cluster_id = c(1:J), Mj, stratum_id)
    } else if (size_model == "poisson") {
      Mj <- rpois(J, 500)
      # make sure we don't get any cluster sizes smaller than 325 (max number of
      # births we are gonna sample)
      while (min(Mj) < 325) {
        Mj <- rpois(J, 500)
      }
      stratum_id <- rep(1, length(Mj))
      # Collect all into data frame
      pop_dat <- data.frame(cluster_id = c(1:J), Mj, stratum_id)
    } else { # size_model == "ffstrat" or "ff"
      codedir <- "/vega/stats/users/smm2253/cluster_sampling/src/simulation"
      pop_dat <- read.csv(paste0(codedir, "/observed_city_pops.csv"), header = TRUE)
      pop_dat$Mj <- as.integer(round(pop_dat$population/1000))
      pop_dat$stratum_id <- pop_dat$stratum # have their own numbering scheme
      # make a shorter version that only has 2 strata
      pop_dat$stratum_id[pop_dat$stratum_id != 999] <- 1
      pop_dat$stratum_id[pop_dat$stratum_id == 999] <- 2
      # renumber 1-num_strat -- checked that this makes the non-extreme stratum
      # be number 9, which is what we want (not needed for short stratum id version)
      #stratum_id <- as.integer(factor(stratum_id)) 
      pop_dat$tot_births <- pop_dat$ubrth + pop_dat$mbrth
    }

    # Calculate log cluster size, generate unit-level covariate
    pop_dat$logMj_c <- log(pop_dat$Mj) - mean(log(pop_dat$Mj))
    if (grepl("ff", size_model)) {
      x <- base::sample(c(unitcovar_range[1]:unitcovar_range[2]),
                        sum(pop_dat$tot_births), replace = TRUE)
    } else {
      x <- base::sample(c(unitcovar_range[1]:unitcovar_range[2]),
                        sum(pop_dat$Mj),
                        replace = TRUE)
    }
    x <- x - mean(x)

    # Draw hyperparameters, varying slopes and coefficients, and outcomes
    if (outcome_type == "continuous") {
      alpha0 <- rnorm(1)
      alpha1 <- rnorm(1)
      if (use_sizes == 1) {
        gamma0 <- rnorm(1)
        gamma1 <- rnorm(1)
      } else {
        gamma0 <- 0
        gamma1 <- 0
      }
      sigma_beta0 <- abs(rnorm(1, 0, 0.5))
      sigma_beta1 <- abs(rnorm(1, 0, 0.5))
      sigma_y <- abs(rnorm(1, 0, 0.75))
        
      pop_dat$beta0 <- rnorm(n = J, mean = alpha0 + gamma0 * pop_dat$logMj_c,
                             sd = sigma_beta0)
      pop_dat$beta1 <- rnorm(n = J, mean = alpha1 + gamma1 * pop_dat$logMj_c,
                             sd = sigma_beta1)
      if (grepl("ff", size_model)) {
        beta0_rep <- rep(pop_dat$beta0, pop_dat$tot_births)
        beta1_rep <- rep(pop_dat$beta1, pop_dat$tot_births)
      } else {
        beta0_rep <- rep(pop_dat$beta0, pop_dat$Mj)
        beta1_rep <- rep(pop_dat$beta1, pop_dat$Mj)
      }
      ymean <- beta0_rep + beta1_rep * x 
      y <- rnorm(ymean, mean = ymean, sd = sigma_y)
    } else {
      alpha0 <- rnorm(1)
      if (use_sizes == 1) {
        gamma0 <- rnorm(1)
      }  else {
        gamma0 <- 0
      }
      sigma_beta0 <- abs(rnorm(1, 0, 0.5))
        
      pop_dat$beta0 <- rnorm(n = J, mean = alpha0 + gamma0 * pop_dat$logMj_c, sd = sigma_beta0)
      if (grepl("ff", size_model)) {
        beta0_rep <- rep(pop_dat$beta0, pop_dat$tot_births)
      } else {
        beta0_rep <- rep(pop_dat$beta0, pop_dat$Mj)
      }
      y_prob <- inv_logit(beta0_rep) 
      y <- rbinom(y_prob, size = 1, prob = y_prob)
      # make x = 0 everywhere so we can still use svy_ests.R
      x_new <- rep(0, length(x))
      x <- x_new
    }

    # Make data frame of pop data -- pop_dat and pop_data were bad name choices
    # but not changing them anymore
    cluster_level_data <- pop_dat
    if (grepl("ff", size_model)) {
      pop_data <- data.frame(y, x,
                             logMj_c = rep(cluster_level_data$logMj_c,
                                           cluster_level_data$tot_births),
                             stratum_id = rep(cluster_level_data$stratum_id, 
                                              cluster_level_data$tot_births),
                             Mj = rep(cluster_level_data$Mj,
                                      cluster_level_data$tot_births),
                             tot_births = rep(cluster_level_data$tot_births,
                                              cluster_level_data$tot_births))
      pop_data$cluster_id <- rep(c(1:J), times = cluster_level_data$tot_births)
      pop_data$unit_id <- unlist(lapply(cluster_level_data$tot_births, seq_len))
    } else {
      pop_data <- data.frame(y, x,
                             logMj_c = rep(cluster_level_data$logMj_c,
                                           cluster_level_data$Mj),
                             Mj = rep(cluster_level_data$Mj,
                                      cluster_level_data$Mj))
      pop_data$cluster_id <- rep(c(1:J), times = cluster_level_data$Mj)
      pop_data$unit_id <- unlist(lapply(cluster_level_data$Mj, seq_len))
    }

    # Make list of things to save
    ybar_true <- mean(pop_data$y)
    if (outcome_type == "continuous") {
      truepars <- data.frame(alpha0, gamma0, alpha1, gamma1, sigma_beta0,
                             sigma_beta1, sigma_y, ybar_true)
    } else {
      truepars <- data.frame(alpha0, gamma0, sigma_beta0, ybar_true)
      beta1 <- NA
    }
   # if (grepl("ff", size_model)) {
   #   popdata <- list(pop_data = pop_data, J = J, Mj = cluster_level_data$Mj,
   #                   logMj_c = cluster_level_data$logMj_c,
   #                   tot_births = cluster_level_data$tot_births,
   #                   stratum_id = cluster_level_data$stratum_id,
   #                   beta0 = cluster_level_data$beta0,
   #                   beta1 = cluster_level_data$beta1, truepars = truepars)
   # } else {
   #   popdata <- list(pop_data = pop_data, J = J, Mj = cluster_level_data$Mj,
   #                   logMj_c = cluster_level_data$logMj_c,
   #                   beta0 = cluster_level_data$beta0,
   #                   beta1 = cluster_level_data$beta1, truepars = truepars)
   # }
    popdata <- list(pop_data = pop_data, truepars = truepars,
                    cluster_level_data = cluster_level_data)
    print(str(popdata))
    saveRDS(popdata,
            file = paste0(rootdir, "/output/simulation/popdata_usesizes_",
                          use_sizes, "_", outcome_type, "_", size_model, ".rds"))

