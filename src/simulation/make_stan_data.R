make_stan_data <- function(stanmod_name, sample_data, pop_cluster_data,
                           sampled_cluster_data) {
  # stanmod_name -- string of the name of the stan model to be used
  # sample_data -- unit-level data frame of the sampled data
  # pop_cluster_data -- cluster-level data frame for all pop clusters
  # sampled_cluster_data -- cluster-level data frame for sampled clusters

  # SORT pop and sample data by cluster_id
  pop_cluster_data <- dplyr::arrange(pop_cluster_data, cluster_id)
  sample_data <- dplyr::arrange(sample_data, cluster_id)

  # Some constants
  J <- length(pop_cluster_data$cluster_id)
  K <- length(sampled_cluster_data$cluster_id)

  # Take things directly from sample_data and pop_cluster_data
  x <- sample_data$x
  y <- sample_data$y
  cluster_id <- sample_data$cluster_id
  n <- nrow(sample_data)
  xbar_pop <- pop_cluster_data$xbar_pop

  # Get sample data at cluster level
  if (grepl("strat", stanmod_name)) {
    stratum_id <- sampled_cluster_data$stratum_id
  }
  sampled_cluster_data <- dplyr::arrange(sampled_cluster_data, cluster_id)
  Mj_sample <- sampled_cluster_data$Mj
  Nj_sample <- sampled_cluster_data$Nj
  log_Mj_sample <- sampled_cluster_data$logMj_c

  ##########################################
  ### Make inputs to stan -- data list, parameters, functions
  ##########################################
  expose_stan_functions(stanmod)
  parlist <- c("alpha0", "gamma0", "alpha1", "gamma1",
               "sigma_beta0", "sigma_beta1", "sigma_y")
  betapars <- c("beta0", "beta1")
  
  switch(stanmod_name,
         bb = {
            M_tot <- sum(pop_cluster_data$Mj)
            n_dat <- sampled_cluster_data %>%
              dplyr::filter(is_certainty_cluster == FALSE) %>%
              dplyr::arrange(cluster_id) %>%
              dplyr::group_by(cluster_id, Mj) %>%
              dplyr::summarise(n = n_distinct(cluster_id))
            num_uniq_sz <- nrow(n_dat)      # number of unique cluster sizes
            vec_uniq_sz <- n_dat$Mj # vector of unique cluster sizes
            cts_uniq_sz <- n_dat$n   # counts of unique cluster sizes
            standata <- list(J = J,
                             K = K,
                             n = n,
                             M_tot = M_tot,
                             num_uniq_sz = num_uniq_sz,
                             x = x,
                             y = y,
                             cluster_id = cluster_id,
                             Mj_sample = Mj_sample,
                             log_Mj_sample = log_Mj_sample,
                             vec_uniq_sz = vec_uniq_sz,
                             cts_uniq_sz = cts_uniq_sz)
         },
         bb_ff = {
            M_tot <- sum(pop_cluster_data$Mj)
            n_dat <- sampled_cluster_data %>%
              dplyr::filter(is_certainty_cluster == FALSE) %>%
              dplyr::arrange(cluster_id) %>%
              dplyr::group_by(cluster_id, Mj) %>%
              dplyr::summarise(n = n_distinct(cluster_id))
            num_uniq_sz <- nrow(n_dat)      # number of unique cluster sizes
            vec_uniq_sz <- n_dat$Mj # vector of unique cluster sizes
            cts_uniq_sz <- n_dat$n   # counts of unique cluster sizes
            standata <- list(J = J,
                             K = K,
                             n = n,
                             M_tot = M_tot,
                             num_uniq_sz = num_uniq_sz,
                             x = x,
                             y = y,
                             cluster_id = cluster_id,
                             Mj_sample = Mj_sample,
                             log_Mj_sample = log_Mj_sample,
                             vec_uniq_sz = vec_uniq_sz,
                             cts_uniq_sz = cts_uniq_sz)
         },
         bb_binary = {
            M_tot <- sum(pop_cluster_data$Mj)
            n_dat <- sampled_cluster_data %>%
              dplyr::filter(is_certainty_cluster == FALSE) %>%
              dplyr::arrange(cluster_id) %>%
              dplyr::group_by(cluster_id, Mj) %>%
              dplyr::summarise(n = n_distinct(cluster_id))
            num_uniq_sz <- nrow(n_dat)      # number of unique cluster sizes
            cts_uniq_sz <- n_dat$n   # counts of unique cluster sizes
            vec_uniq_sz <- n_dat$Mj # vector of unique cluster sizes
            standata <- list(J = J,
                             K = K,
                             n = n,
                             M_tot = M_tot,
                             num_uniq_sz = num_uniq_sz,
                             y = y,
                             cluster_id = cluster_id,
                             Mj_sample = Mj_sample,
                             log_Mj_sample = log_Mj_sample,
                             vec_uniq_sz = vec_uniq_sz,
                             cts_uniq_sz = cts_uniq_sz)
            betapars <- "beta0"
         },
         cluster_inds_only = {
            standata <- list(J = J,
                             K = K,
                             n = n,
                             x = x,
                             y = y,
                             cluster_id = cluster_id)
            parlist <- c("alpha0", "alpha1",
                         "sigma_beta0", "sigma_beta1", "sigma_y")
         },
         cluster_inds_only_binary = {
            standata <- list(J = J,
                             K = K,
                             n = n,
                             y = y,
                             cluster_id = cluster_id)
            parlist <- c("alpha0", "sigma_beta0")
            betapars <- "beta0"
         },
         knowsizes = {
            standata <- list(J = J,
                             K = K,
                             n = n,
                             x = x,
                             y = y,
                             cluster_id = cluster_id,
                             Mj_sample = Mj_sample,
                             log_Mj_sample = log_Mj_sample)
         },
         knowsizes_binary = {
            standata <- list(J = J,
                             K = K,
                             n = n,
                             y = y,
                             cluster_id = cluster_id,
                             Mj_sample = Mj_sample,
                             log_Mj_sample = log_Mj_sample)
            betapars <- "beta0"
         },
         lognormal = {
            standata <- list(J = J,
                             K = K,
                             n = n,
                             x = x,
                             y = y,
                             cluster_id = cluster_id,
                             Mj_sample = Mj_sample,
                             log_Mj_sample = log_Mj_sample)
            sbpars <- c("mu", "sigma")
            parlist <- c(parlist, c("mu", "sigma"))
         },
         lognormal_ff = {
            standata <- list(J = J,
                             K = K,
                             n = n,
                             x = x,
                             y = y,
                             cluster_id = cluster_id,
                             Mj_sample = Mj_sample,
                             log_Mj_sample = log_Mj_sample)
            sbpars <- c("mu", "sigma")
            parlist <- c(parlist, c("mu", "sigma"))
         },
         lognormal_binary = {
            standata <- list(J = J,
                             K = K,
                             n = n,
                             y = y,
                             cluster_id = cluster_id,
                             Mj_sample = Mj_sample,
                             log_Mj_sample = log_Mj_sample)
            betapars <- "beta0"
            sbpars <- c("mu", "sigma")
            parlist <- c(parlist, c("mu", "sigma"))
         },
         negbin = {
            standata <- list(J = J,
                             K = K,
                             n = n,
                             x = x,
                             y = y,
                             cluster_id = cluster_id,
                             Mj_sample = Mj_sample,
                             log_Mj_sample = log_Mj_sample)
            sbpars <- c("mu", "phi")
            parlist <- c(parlist, c("mu", "phi"))
         },
         negbin_ff = {
            standata <- list(J = J,
                             K = K,
                             n = n,
                             x = x,
                             y = y,
                             cluster_id = cluster_id,
                             Mj_sample = Mj_sample,
                             Nj_sample = Nj_sample,
                             log_Mj_sample = log_Mj_sample)
            parlist <- c(parlist, c("mu", "phi"))
         },
         negbin_binary = {
            standata <- list(J = J,
                             K = K,
                             n = n,
                             y = y,
                             cluster_id = cluster_id,
                             Mj_sample = Mj_sample,
                             log_Mj_sample = log_Mj_sample)
            betapars <- "beta0"
            sbpars <- c("mu", "phi")
            parlist <- c(parlist, c("mu", "phi"))
            print("parlist after negbin_bin in make stan data")
            print(parlist)
         },
         negbin_ffstrat = {
            is_strat2 <- as.numeric(sampled_cluster_data$stratum_id == 2)
            standata <- list(J = J,
                             K = K,
                             n = n,
                             x = x,
                             y = y,
                             cluster_id = cluster_id,
                             Mj_sample = Mj_sample,
                             Nj_sample = Nj_sample,
                             log_Mj_sample = log_Mj_sample,
                             is_strat2 = is_strat2)
            sbpars <- c("mu", "phi")
            parlist <- c(parlist, c("mu", "phi", "a", "b", "kappa0", "kappa1"))
         },
         stop("Invalid value of stanmod_name")
  )
  
  #  if (grepl("strat", stanmod_name)) {
  #    stratum_matrix <- model.matrix(~ 0 + factor(stratum_id), sampled_cluster_data)
  #    stratum_matrix_pop <- model.matrix(~ 0 + factor(stratum_id), pop_cluster_data)
  #    standata <- c(standata, list(stratum_id = stratum_id, S = 2,
  #                                 stratum_matrix = stratum_matrix))
  #    #parlist <- c(parlist, "mu", "phi")
  #  } #else { # don't include these for ff because they're vectors in it
  #    #parlist <- c(parlist, "mu_star", "phi_star", "mu", "phi")
  #  #}
  #
  # don't need alpha1, gamma1, sigma_beta1, etc in parlist when y is binary
  if (outcome_type == "binary") {
    betapars <- "beta0"
    kappapars <- "kappa0"
  print("parlist before")
  print(parlist)
  print("parlist %in% c(sigma_y, alpha1, gamma1, sigma_beta1)")
  print(parlist %in% c("sigma_y", "alpha1", "gamma1", "sigma_beta1"))
    parlist <- setdiff(parlist, c("sigma_y", "alpha1", "gamma1", "sigma_beta1"))
  print("parlist after")
  print(parlist)
  }

  print("parlist")
  print(parlist)
  print(str(standata))
  rm(sample_data)
  
  to_return <- list(standata = standata, parlist = parlist, betapars = betapars)
  if (exists("sbpars")) {
    to_return <- c(to_return, list(sbpars = sbpars))
  }
  if (exists("kappapars")) {
    to_return <- c(to_return, list(kappapars = kappapars))
  }
 
  
  return(to_return)
}

