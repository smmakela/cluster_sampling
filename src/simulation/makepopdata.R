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
  ### READ ME
  #############################################################################
  # This code has been changed to better distinguish between the cluster size,
  # Nj, and the measure of size, Mj. For size_model = poisson and multinomial,
  # Nj = Mj, so they are the same. BUT, crucially, the ARE NOT THE SAME for
  # the ff populations! In that case, Mj is the city population, but the
  # relevant Nj is the total number of births.

  #############################################################################
  ### Setup of libraries and directories
  #############################################################################
    libdir <- "/vega/stats/users/smm2253/rpackages"
    .libPaths(libdir)
    library(dplyr)
    library(data.table)
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
  ### Useful functions
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
  ### 1. multinomial draws with gamma-xformed probs (dirichlet)
  ### 2. Poisson
  ### 3. FF observed city sizes
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
      pop_cluster_data <- data.frame(cluster_id = c(1:J), Mj, stratum_id)
      pop_cluster_data$Nj <- pop_cluster_data$Mj
    } else if (size_model == "poisson") {
      Mj <- rpois(J, 500)
      # make sure we don't get any cluster sizes smaller than 325 (max number of
      # births we are gonna sample)
      while (min(Mj) < 325) {
        Mj <- rpois(J, 500)
      }
      stratum_id <- rep(1, length(Mj))
      # Collect all into data frame
      pop_cluster_data <- data.frame(cluster_id = c(1:J), Mj, stratum_id)
      pop_cluster_data$Nj <- pop_cluster_data$Mj
    } else { # size_model == "ffstrat" or "ff"
      codedir <- "/vega/stats/users/smm2253/cluster_sampling/src/simulation"
      pop_cluster_data <- read.csv(paste0(codedir, "/observed_city_pops.csv"), header = TRUE)
      pop_cluster_data$Mj <- as.integer(round(pop_cluster_data$population/100))
      pop_cluster_data$stratum_id <- pop_cluster_data$stratum # have their own numbering scheme
      # make a shorter version that only has 2 strata, along with the number of
      # units to sample in each
      pop_cluster_data$stratum_id[pop_cluster_data$stratum_id != 999] <- 1
      pop_cluster_data$stratum_id[pop_cluster_data$stratum_id == 999] <- 2
      pop_cluster_data$num_units_to_sample <- ifelse(pop_cluster_data$stratum_id == 2, 100, 325)
      # renumber 1-num_strat -- checked that this makes the non-extreme stratum
      # be number 9, which is what we want (not needed for short stratum id version)
      #stratum_id <- as.integer(factor(stratum_id)) 
      pop_cluster_data$tot_births <- pop_cluster_data$ubrth + pop_cluster_data$mbrth
      #pop_cluster_data$Nj <- as.integer(round(pop_cluster_data$tot_births))
      pop_cluster_data$Nj <- pop_cluster_data$Mj
      pop_cluster_data$cluster_id <- c(1:J)
    }
    pop_cluster_data <- dplyr::arrange(pop_cluster_data, cluster_id)

    # Calculate log cluster size, generate unit-level covariate
    pop_cluster_data$logMj_c <- log(pop_cluster_data$Mj) - mean(log(pop_cluster_data$Mj))
    x <- base::sample(c(unitcovar_range[1]:unitcovar_range[2]),
                      sum(pop_cluster_data$Nj), replace = TRUE)
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
        
      pop_cluster_data$beta0 <- rnorm(n = J, mean = alpha0 + gamma0 * pop_cluster_data$logMj_c,
                                        sd = sigma_beta0)
      pop_cluster_data$beta1 <- rnorm(n = J, mean = alpha1 + gamma1 * pop_cluster_data$logMj_c,
                                        sd = sigma_beta1)
      beta0_rep <- rep(pop_cluster_data$beta0, pop_cluster_data$Nj)
      beta1_rep <- rep(pop_cluster_data$beta1, pop_cluster_data$Nj)
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
        
      pop_cluster_data$beta0 <- rnorm(n = J, mean = alpha0 + gamma0 * pop_cluster_data$logMj_c, sd = sigma_beta0)
      beta0_rep <- rep(pop_cluster_data$beta0, pop_cluster_data$Nj)
      y_prob <- inv_logit(beta0_rep) 
      y <- rbinom(y_prob, size = 1, prob = y_prob)
      # make x = 0 everywhere so we can still use svy_ests.R
      x_new <- rep(0, length(x))
      x <- x_new
    }

    # Make data frame of pop data (at the unit level)
    pop_data <- data.frame(y, x, cluster_id = rep(c(1:J), times = pop_cluster_data$Nj))
    pop_data <- dplyr::left_join(pop_data, pop_cluster_data, by = "cluster_id")
    ybar_true <- mean(pop_data$y)

    # Calculate mean, sum of x and sum of y by cluster and add to pop_cluster_data
    # (sum of x and y are used in svy_ests for design-based estimates)
    tmp <- pop_data %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(xbar_pop = mean(x),
                       sum_x_pop = sum(x),
                       sum_y_pop = sum(y))
    pop_cluster_data <- dplyr::left_join(pop_cluster_data, tmp, by = "cluster_id")

    # For FF: to save room, sample units in each cluster now, either 325 or 100
    # depending on whether the cluster was originally a high- or low-sample
    # city. NOTE we are not sampling the clusters yet, just the units to
    # save ourselves some storage space
    if (grepl("ff", size_model)) {
      pop_data <- setDT(pop_data)
      pop_data <- pop_data[, .SD[sample(.N, num_units_to_sample)], by = cluster_id]
    }
    # Create unit id (unique within cluster)
    pop_data <- pop_data %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::mutate(unit_id = row_number())

    # Make list of things to save
    if (outcome_type == "continuous") {
      truepars <- data.frame(alpha0, gamma0, alpha1, gamma1, sigma_beta0,
                             sigma_beta1, sigma_y, ybar_true)
    } else {
      truepars <- data.frame(alpha0, gamma0, sigma_beta0, ybar_true)
      beta1 <- NA
    }
    popdata <- list(pop_data = pop_data, truepars = truepars,
                    pop_cluster_data = pop_cluster_data)
    print(str(popdata))
    saveRDS(popdata,
            file = paste0(rootdir, "/output/simulation/popdata_usesizes_",
                          use_sizes, "_", outcome_type, "_", size_model, ".rds"))

