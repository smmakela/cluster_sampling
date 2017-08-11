  ##########################################
  ### Make inputs to stan -- data list, parameters, functions
  ##########################################
  expose_stan_functions(stanmod)
  parlist <- c("alpha0", "gamma0", "alpha1", "gamma1",
               "sigma_beta0", "sigma_beta1", "sigma_y")
  betapars <- c("beta0", "beta1")
print("THREE")
  if (grepl("strat", stanmod_name)) {
    kappapars <- c("kappa0", "kappa1")
    if (grepl("lognormal", stanmod_name)) {
      sbpars <- c("mu", "phi")
    }
    if (grepl("negbin", stanmod_name)) {
      sbpars <- c("mu", "phi")
    }
  }
  if (grepl("bb", stanmod_name)) {
    standata <- list(J = J,
                     K = K,
                     n = n,
                     N = N,
                     M = M,
                     x = x,
                     y = y,
                     cluster_id = cluster_id,
                     Nj_sample = Nj_sample,
                     log_Nj_sample = log_Nj_sample,
                     Nj_unique = Nj_unique,
                     M_counts = M_counts)
  } else if (grepl("cluster_inds_only", stanmod_name)) {
    standata <- list(J = J,
                     K = K,
                     n = n,
                     x = x,
                     y = y,
                     cluster_id = cluster_id)
    parlist <- c("alpha0", "alpha1", "sigma_beta0", "sigma_beta1", "sigma_y")
  } else if (grepl("knowsizes", stanmod_name)) {
    standata <- list(J = J,
                     K = K,
                     n = n,
                     x = x,
                     y = y,
                     cluster_id = cluster_id,
                     Nj_sample = Nj_sample,
                     log_Nj_sample = log_Nj_sample)
  } else if (grepl("lognormal", stanmod_name)) {
    standata <- list(J = J,
                     K = K,
                     n = n,
                     x = x, 
                     y = y,
                     cluster_id = cluster_id,
                     Nj_sample = Nj_sample,
                     log_Nj_sample = log_Nj_sample)
    if (!grepl("strat", stanmod_name)) {
      parlist <- c(parlist, "mu_star", "sigma", "mu", "mu_star_scaled", "sigma_scaled")
    }
  } else if (grepl("negbin", stanmod_name)) {
    standata <- list(J = J,
                     K = K,
                     n = n,
                     x = x, 
                     y = y,
                     cluster_id = cluster_id,
                     Nj_sample = Nj_sample,
                     log_Nj_sample = log_Nj_sample)
    if (size_model == "ff") {
      standata <- c(standata, list(Mj_sample = Mj_
    if (!grepl("strat", stanmod_name)) {
      parlist <- c(parlist, "mu", "phi")
    }
  } else {
    stop("Invalid stan model")
  }
  #if (size_model == "ffstrat") {
  if (grepl("strat", stanmod_name)) {
print("FIVE")
    stratum_matrix <- model.matrix(~ 0 + factor(stratum_id), sampled_cluster_data)
    stratum_matrix_pop <- model.matrix(~ 0 + factor(stratum_id), cluster_data_pop)
    standata <- c(standata, list(stratum_id = stratum_id, S = 2,
                                 stratum_matrix = stratum_matrix))
    #parlist <- c(parlist, "mu", "phi")
  } #else { # don't include these for ff because they're vectors in it
    #parlist <- c(parlist, "mu_star", "phi_star", "mu", "phi")
  #}

  # don't need x in standata when outcome_type is binary (for now), and also
  # don't need alpha1, gamma1, sigma_beta1, etc in parlist
  if (outcome_type == "binary") {
    standata$x <- NULL
    betapars <- "beta0"
    kappapars <- "kappa0"
    plist <- c("alpha1", "gamma1", "sigma_beta1", "sigma_y")
    parlist <- parlist[!(parlist %in% plist)]
  }

  print(str(standata))
  rm(pop_data)
  rm(sample_data)

