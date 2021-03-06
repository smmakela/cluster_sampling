# Author: Susanna Makela
# Date: 21 Apr 2014
# Purpose: run stan

runstan <- function(num_clusters, num_units, use_sizes, outcome_type, size_model, rootdir,
                    simno, stanmod, stanmod_name, sim_data, num_iter, num_chains) {
  # num_clusters -- number of clusters to sample
  # num_units -- number of units to sample
  # use_sizes -- 0/1 for whether cluster sizes used in pop data
  # outcome_type -- whether outcome is continuous or binary
  # rootdir -- root directory where src, output, etc folders are
  # simno -- current iteration; used so that multiple instances aren't trying to write to the same file
  # stanmod -- compiled stan model
  # stanmod_name -- string for name of stan model so we know which parts of the code to run
  # sim_data -- data to use for simulation
  # num_iter -- number of iterations stan should run for
  # num_chains -- number of chains to run in stan

  ##########################################
  ### Load pop and sample data
  ##########################################
  rstan_options(auto_write = TRUE)
  options(mc_cores = parallel::detectCores())

  if (num_units <= 1) {
    nunits <- paste(100*num_units, "pct", sep = "")
  } else {
    nunits <- num_units
  }
print("str(sim_data)")
print(str(sim_data))
  for (j in names(sim_data)) {
    assign(j, sim_data[[j]])
  }
  rm(sim_data)
  ybar_true <- mean(pop_data$y)

  model_name <- gsub("_binary", "", stanmod_name)

  ##########################################
  ### Make data for stan
  ##########################################
  print("making stan data")
  print(Sys.time())

  if (num_units != 999) {
    K <- num_clusters
  } else {
    K <- J
  }

  # SORT pop and sample data by cluster_id
  pop_data <- dplyr::arrange(pop_data, cluster_id)
  sample_data <- dplyr::arrange(sample_data, cluster_id)

  # Take things directly from sample_data
  x <- sample_data$x
  y <- sample_data$y
  cluster_id <- sample_data$cluster_id
  n <- nrow(sample_data)

  #Nj_pop <- Mj
  #N <- sum(Nj_pop)
  #n <- sum(pop_data$insample)
  #tmp <- dplyr::summarise(group_by(pop_data, cluster_id),
  #                        xbar_pop = mean(x))
  #tmp <- dplyr::arrange(tmp, cluster_id)
  #xbar_pop <- tmp$xbar_pop
  #x <- sample_data$x
  #y <- sample_data$y
  #cluster_id <- sample_data$cluster_id

  # Get pop data at cluster level -- calculate mean of x, grouping by cluster_id,
  # stratum_id, Mj, logMj_c (cluster_id is the most detailed)
  #if (size_model == "ffstrat") {
print("ONE")
  if (grepl("strat", stanmod_name)) {
    cluster_data_pop <- pop_data %>%
      dplyr::select(cluster_id, stratum_id, Mj, logMj_c, x) %>%
      dplyr::group_by(cluster_id, stratum_id, Mj, logMj_c) %>%
      dplyr::summarise(xbar_pop = mean(x)) %>%
      dplyr::arrange(cluster_id)
  } else {
    cluster_data_pop <- pop_data %>%
      dplyr::select(cluster_id, Mj, logMj_c, x) %>%
      dplyr::group_by(cluster_id, Mj, logMj_c) %>%
      dplyr::summarise(xbar_pop = mean(x)) %>%
      dplyr::arrange(cluster_id)
  }
print("str(cluster_data_pop)")
print(str(cluster_data_pop))
  Nj_pop <- cluster_data_pop$Mj
  N <- sum(Nj_pop)
  xbar_pop <- cluster_data_pop$xbar_pop

  # Get sample data at cluster level
print("TWO")
  if (grepl("strat", stanmod_name)) {
    stratum_id <- sampled_cluster_data$stratum_id
  }
  sampled_cluster_data <- dplyr::arrange(sampled_cluster_data, cluster_id)
  Nj_sample <- sampled_cluster_data$Mj
  log_Nj_sample <- sampled_cluster_data$logMj_c

  # special data for bb
  n_dat <- sampled_cluster_data %>%
    dplyr::arrange(cluster_id) %>%
    dplyr::group_by(cluster_id, Mj) %>%
    dplyr::summarise(n = n_distinct(cluster_id))
  M <- nrow(n_dat)      # number of unique cluster sizes
  M_counts <- n_dat$n   # counts of unique cluster sizes
  Nj_unique <- n_dat$Mj # vector of unique cluster sizes

  # delete pop and sample data to save memory
  #rm(pop_data)
  #rm(sample_data)
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
print("FOUR")
    standata <- list(J = J,
                     K = K,
                     n = n,
                     x = x, 
                     y = y,
                     cluster_id = cluster_id,
                     Nj_sample = Nj_sample,
                     log_Nj_sample = log_Nj_sample)
    #if (size_model != "ffstrat") {
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

  ##########################################
  ### Run stan
  ##########################################
  print(Sys.time())

  fit <- sampling(stanmod, data = standata,
                  iter = num_iter, chains = num_chains,
                  control = list(stepsize = 0.001, adapt_delta = 0.999))
  print("done fitting stan model")
  print(Sys.time())
  print(warnings())

  ##########################################
  ### See if there are any divergent transitions -- if so need to rerun 
  ##########################################
    samp_params <- get_sampler_params(fit)
    num_div_trans <- 0
    for (s in 1:length(samp_params)) {                                
      num_div_trans <- num_div_trans + sum(samp_params[[s]][((num_iter/2)+1):num_iter, "divergent__"])
    }
    ad_val <- 0.999
    if (num_div_trans > 0) {
      counter <- 1
      while (counter <= 3 & num_div_trans > 0) {
        ad_val <- ad_val + 9 * (10^(-1*(counter + 3)))
        cat("----------------------------------------------------------------\n")
        cat("Divergent transitions happened, rerunning with adapt_delta =", ad_val, "\n")
        cat("num_clusters =", num_clusters, "num_units=", num_units, "\n")
        fit <- sampling(stanmod, data = standata,
                        iter = num_iter, chains = num_chains,
                        control = list(stepsize = 0.001, adapt_delta = ad_val))
        print("done fitting stan model")
        print(Sys.time())
        print(warnings())
        samp_params <- get_sampler_params(fit)
        num_div_trans <- 0
        for (s in 1:length(samp_params)) {                                
          num_div_trans <- num_div_trans + sum(samp_params[[s]][((num_iter/2)+1):num_iter, "divergent__"])
        }
        counter <- counter + 1
      } # end while
      # if we still have divergent transitions, try the ncp version (only for continuous)
      if (num_div_trans > 0) {
        cat("----------------------------------------------------------------\n")
        cat("Divergent transitions remain, running NCP version\n")
        cat("num_clusters =", num_clusters, "num_units=", num_units, "\n")
        stanmod <- readRDS(paste0(rootdir, "/src/analysis/", stanmod_name, "_ncp.rds"))
        expose_stan_functions(stanmod)
        fit <- sampling(stanmod, data = standata,
                      iter = num_iter, chains = num_chains,
                      control = list(stepsize = 0.001, adapt_delta = 0.999))
        print("done fitting stan model")
        print(Sys.time())
        print(warnings())
        samp_params <- get_sampler_params(fit)
        num_div_trans <- 0
        for (s in 1:length(samp_params)) {                                
          num_div_trans <- num_div_trans + sum(samp_params[[s]][((num_iter/2)+1):num_iter, "divergent__"])
        }
        # try raising adapt_delta if necessary
        ad_val <- 0.999
        if (num_div_trans > 0) {
          counter <- 1
          while (counter <= 3 & num_div_trans > 0) {
            ad_val <- ad_val + 9 * (10^(-1*(counter + 3)))
            cat("----------------------------------------------------------------\n")
            cat("Divergent transitions happened for NCP, rerunning with adapt_delta =", ad_val, "\n")
            cat("num_clusters =", num_clusters, "num_units=", num_units, "\n")
            fit <- sampling(stanmod, data = standata,
                            iter = num_iter, chains = num_chains,
                            control = list(stepsize = 0.001, adapt_delta = ad_val))
            print("done fitting stan model")
            print(Sys.time())
            print(warnings())
            samp_params <- get_sampler_params(fit)
            num_div_trans <- 0
            for (s in 1:length(samp_params)) {                                
              num_div_trans <- num_div_trans + sum(samp_params[[s]][((num_iter/2)+1):num_iter, "divergent__"])
            }
            counter <- counter + 1
          } # end while
          # if we *still* have divergent transitions, give up
          if (num_div_trans > 0) {
            cat("Unable to get rid of divergent transitions :( \n")
            cat("num_clusters =", num_clusters, "num_units=", num_units, "\n")
            out_msg <- paste0("There were ", num_div_trans, " divergent transitions_")
            saveRDS(out_msg,
                    paste0(rootdir, "output/simulation/stan_results_usesizes_",
                           use_sizes, "_", outcome_type, "_", size_model, "_",
                           model_name, "_nclusters_", num_clusters,
                           "_nunits_", nunits, "_sim_", simno, ".rds"))
            saveRDS(fit,
                    paste0(rootdir, "/output/simulation/stanfit_usesizes_",
                           use_sizes, "_", outcome_type, "_", size_model,
                           "_", model_name, "_nclusters_", num_clusters,
                           "_nunits_", nunits, "_simno_", simno,".rds"))
            return(NULL)
          } # end if for num_div_trans > 0 after ncp
        } # end if for raising adapt_delta in ncp
      } # end if for doing ncp
    } # end if for initial num_div_trans > 0

    if (simno == 1) {
      saveRDS(fit, file = paste0(rootdir, "/output/simulation/stanfit_usesizes_",
                                 use_sizes, "_", outcome_type, "_", size_model,
                                  "_", model_name, "_nclusters_", num_clusters,
                                 "_nunits_", nunits, "_simno_", simno,".rds"))
      print("done saving stanfit object")
      cat("filename:", paste0(rootdir, "/output/simulation/stanfit_usesizes_",
                                 use_sizes, "_", outcome_type, "_", size_model,
                                  "_", model_name, "_nclusters_", num_clusters,
                                 "_nunits_", nunits, "_simno_", simno,".rds"))
      print(Sys.time())
    }

  # Make sure things converged in terms of Rhat too
  rhat_check_df <- data.frame(summary(fit)$summary)
  rhat_check <- rhat_check_df$Rhat
  if (max(rhat_check) >= 1.1) {
    out_msg <- paste0("Rhat check failed! Max Rhat is", max(rhat_check))
    cat("Rhat check for alpha0 failed! Rhat is", max(rhat_check), "\n")
    print(rhat_check_df)
    print(summary(rhat_check_df[!grepl("\\[", rownames(rhat_check_df)), ]))
    saveRDS(out_msg,
            paste0(rootdir, "output/simulation/stan_results_usesizes_",
                   use_sizes, "_", outcome_type, "_", size_model, "_",
                   model_name, "_nclusters_", num_clusters,
                   "_nunits_", nunits, "_sim_", simno, ".rds"))
    saveRDS(fit, file = paste0(rootdir, "/output/simulation/stanfit_usesizes_",
                               use_sizes, "_", outcome_type, "_", size_model,
                                "_", model_name, "_nclusters_", num_clusters,
                               "_nunits_", nunits, "_simno_", simno,".rds"))
    return(NULL)
  }

  # Print summary for the basic parameters
  print(fit, pars = parlist)
  if (grepl(stanmod_name, "lognormal")) {
    cat("sample mean, sd of log(Nj_sample):", mean(log(Nj_sample)), sd(log(Nj_sample)), "\n")
  }
  #if (outcome_type == "continuous") {
  #  if (use_sizes == 1) {
  #    print(fit, pars = c("alpha0", "gamma0", "alpha1", "gamma1",
  #                        "sigma_beta0", "sigma_beta1", "sigma_y"))
  #  } else {
  #    print(fit, pars = c("alpha0", "gamma0", "sigma_beta1", "sigma_y"))
  #  }
  #} else {
  #  if (use_sizes == 1) {
  #    print(fit, pars = c("alpha0", "gamma0", "alpha1", "gamma1",
  #                        "sigma_beta0", "sigma_beta1"))
  #  } else {
  #    print(fit, pars = c("alpha0", "gamma0", "sigma_beta1"))
  #  }
  #}

  ##########################################
  ### Extract samples, format param estimates to just keep necessary info 
  ##########################################
  param_samps <- data.frame(rstan::extract(fit, pars = parlist))

  # since the betas have to be passed in as a vector, deal with them separately
  beta_samps_orig <- data.frame(rstan::extract(fit, pars = betapars))
  nc <- nchar("beta0")
  beta_samps_orig %>%
      tidyr::gather(key = pname, value = value) %>%
      dplyr::mutate(par = substr(pname, 1, nc),
                    cluster_id = as.numeric(substr(pname, nc+2, nchar(pname))),
                    pname = NULL) %>%
      dplyr::group_by(par, cluster_id) %>%
      dplyr::mutate(draw_num = row_number()) %>%
      tidyr::spread(key = par, value = value) -> beta_samps

  # same for the kappas and sb pars, if we're doing ff
  #if (size_model == "ffstrat") {
  if (grepl("strat", stanmod_name)) {
    kappa_samps_orig <- data.frame(rstan::extract(fit, pars = kappapars))
    # this is for when we only have two strata
    kappa_samps <- kappa_samps_orig %>%
      dplyr::mutate(draw_num = row_number(),
                    stratum_id = 2)
    kappa_tmp <- kappa_samps
    kappa_tmp$stratum_id <- 1
    kappa_tmp$kappa0 <- 0.0
    kappa_tmp$kappa1 <- 0.0
    kappa_samps <- rbind(kappa_samps, kappa_tmp)
    #nck <- nchar("kappa0")
    #kappa_samps <- kappa_samps_orig %>%
    #    tidyr::gather(key = pname, value = value) %>%
    #    dplyr::mutate(par = substr(pname, 1, nck),
    #                  cluster_id = as.numeric(substr(pname, nck+2, nchar(pname))),
    #                  pname = NULL) %>%
    #    dplyr::group_by(par, cluster_id) %>%
    #    dplyr::mutate(draw_num = row_number()) %>%
    #    tidyr::spread(key = par, value = value) 
    if (grepl("lognormal", stanmod_name) || grepl("negbin", stanmod_name)) {
      sb_samps_orig <- data.frame(rstan::extract(fit, pars = sbpars))
      sb_samps <- sb_samps_orig %>%
          dplyr::mutate(draw_num = row_number()) %>%
          tidyr::gather(key = pname, value = value, -draw_num) %>%
          tidyr::separate(col = pname, into = c("par", "stratum_id"), sep = "\\.") %>%
          dplyr::mutate(stratum_id = as.numeric(stratum_id)) 
      sb_samps <- sb_samps %>%
          tidyr::spread(key = par, value = value)
    }
  }
  # same for phi_star, if we're doing bb
  if (grepl("bb", stanmod_name)) {
    nc2 <- nchar("phi_star")
    phi_star_samps_orig <- data.frame(rstan::extract(fit, pars = "phi_star"))
    print(str(phi_star_samps_orig))
    phi_star_samps_orig %>%
      tidyr::gather(key = pname, value = value) %>%
      dplyr::mutate(par = substr(pname, 1, nc2),
                    cluster_id = as.numeric(substr(pname, nc2+2, nchar(pname))),
                    pname = NULL) %>%
      dplyr::group_by(par, cluster_id) %>%
      dplyr::mutate(draw_num = row_number()) %>%
      tidyr::spread(key = par, value = value) -> phi_star_samps
    print(str(phi_star_samps))
  }

  print("done making samps")
  print(Sys.time())

  tt <- summary(fit)$summary
  #if (size_model != "ffstrat") {
  if (!grepl("strat", stanmod_name)) {
    par_ests <- tt[parlist, ]
  } else {
    ii <- which(grepl(sbpars[1], rownames(tt)) | grepl(sbpars[2], rownames(tt)))
    jj <- which(grepl(kappapars[1], rownames(tt)) | grepl(kappapars[2], rownames(tt)))
    par_ests <- tt[c(parlist, rownames(tt)[c(ii, jj)]), ]
  }
  rm(fit)
  print("done making par_ests")
  print(Sys.time())

  print("making par_ests")
  print(Sys.time())
  par_ests_rownames <- attr(par_ests, "dimnames")[[1]]
  par_ests_rownames <- gsub("\\[", "", par_ests_rownames) 
  par_ests_rownames <- gsub("\\]", "", par_ests_rownames) 
  par_ests_colnames <- attr(par_ests, "dimnames")[[2]]
  par_ests <- data.frame(par_ests, row_names = par_ests_rownames)
  colnames(par_ests) <- par_ests_colnames

  print("printing par_ests")
  print(Sys.time())      
  print(str(par_ests))

  ##########################################
  ### Draw new cluster sizes (if needed) and ybar_new
  ##########################################
  num_draws <- nrow(param_samps)
  Nj_new_df <- expand.grid(cluster_id = c(1:J), draw_num = c(1:num_draws))
  Nj_new_df$Nj_new <- NA
  N_new <- rep(NA, times = num_draws)
  ybar_new <- rep(NA, times = num_draws)
  if (stanmod_name == "bb") {
    for (s in 1:num_draws) {
      beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
      phi_star_samps_sub <- dplyr::filter(phi_star_samps, draw_num == s)
      Nj_new <- Nj_new_bb_rng(J, K, M, Nj_sample, Nj_unique,
                              phi_star_samps_sub$phi_star)
      Nj_new_df$Nj_new[Nj_new_df$draw_num == s] <- Nj_new
      N_new[s] <- sum(Nj_new)
      ybar_new[s] <- ybar_new_bb_rng(J, K, xbar_pop,
                                     beta_samps_sub$beta0,
                                     beta_samps_sub$beta1,
                                     param_samps$alpha0[s],
                                     param_samps$gamma0[s],
                                     param_samps$alpha1[s],
                                     param_samps$gamma1[s],
                                     param_samps$sigma_beta0[s],
                                     param_samps$sigma_beta1[s],
                                     param_samps$sigma_y[s], Nj_new)
    } # end num_draws loop
    #a_value <- 50 # value of a to use for constrained BB
    #Nj_new_df_2 <- Nj_new_df
    #Tx_mis <- sum(Nj_pop) - sum(Nj_sample)
    #Nj_tots <- Nj_new_df %>%
    #  dplyr::group_by(draw_num) %>%
    #  dplyr::summarise(tot = sum(Nj_new)) %>%
    #  dplyr::mutate(abs_diff = abs(Tx_mis - tot)) %>%
    #  dplyr::arrange(abs_diff)
    #keep_inds <- Nj_tots$draw_num[
  } else if (stanmod_name == "bb_binary") {
    for (s in 1:num_draws) {
      beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
      phi_star_samps_sub <- dplyr::filter(phi_star_samps, draw_num == s)
      Nj_new <- Nj_new_bb_rng(J, K, M, Nj_sample, Nj_unique,
                              phi_star_samps_sub$phi_star)
      Nj_new_df$Nj_new[Nj_new_df$draw_num == s] <- Nj_new
      N_new[s] <- sum(Nj_new)
      ybar_new[s] <- ybar_new_bb_rng(J, K,
                                     beta_samps_sub$beta0,
                                     param_samps$alpha0[s],
                                     param_samps$gamma0[s],
                                     param_samps$sigma_beta0[s],
                                     Nj_new)
    } # end num_draws loop
  } else if (stanmod_name == "cluster_inds_only") {
    for (s in 1:num_draws) {
      beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
      ybar_new[s] <- ybar_new_inds_rng(J, K, xbar_pop,
                                       beta_samps_sub$beta0,
                                       beta_samps_sub$beta1,
                                       param_samps$alpha0[s],
                                       param_samps$alpha1[s],
                                       param_samps$sigma_beta0[s],
                                       param_samps$sigma_beta1[s],
                                       param_samps$sigma_y[s], Nj_pop)
    } # end num_draws loop
  } else if (stanmod_name == "cluster_inds_only_binary") {
    for (s in 1:num_draws) {
      beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
      ybar_new[s] <- ybar_new_inds_rng(J, K,
                                       beta_samps_sub$beta0,
                                       param_samps$alpha0[s],
                                       param_samps$sigma_beta0[s],
                                       Nj_pop)
    } # end num_draws loop
  } else if (stanmod_name == "knowsizes") {
    for (s in 1:num_draws) {
      beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
      ybar_new[s] <- ybar_new_know_rng(J, K, xbar_pop,
                                       beta_samps_sub$beta0,
                                       beta_samps_sub$beta1,
                                       param_samps$alpha0[s],
                                       param_samps$gamma0[s],
                                       param_samps$alpha1[s],
                                       param_samps$gamma1[s],
                                       param_samps$sigma_beta0[s],
                                       param_samps$sigma_beta1[s],
                                       param_samps$sigma_y[s], Nj_pop)
    } # end num_draws loop
  } else if (stanmod_name == "knowsizes_binary") {
    for (s in 1:num_draws) {
      beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
      ybar_new[s] <- ybar_new_know_rng(J, K,
                                       beta_samps_sub$beta0,
                                       param_samps$alpha0[s],
                                       param_samps$gamma0[s],
                                       param_samps$sigma_beta0[s],
                                       Nj_pop)
    } # end num_draws loop
  } else if (stanmod_name == "lognormal") {
    for (s in 1:num_draws) {
      beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
      Nj_new <- Nj_new_ln_rng(J, K, Nj_sample, param_samps$mu[s],
                              param_samps$sigma[s])
      Nj_new_df$Nj_new[Nj_new_df$draw_num == s] <- Nj_new
      N_new[s] <- sum(Nj_new)
      ybar_new[s] <- ybar_new_ln_rng(J, K, xbar_pop,
                                     beta_samps_sub$beta0,
                                     beta_samps_sub$beta1,
                                     param_samps$alpha0[s],
                                     param_samps$gamma0[s],
                                     param_samps$alpha1[s],
                                     param_samps$gamma1[s],
                                     param_samps$sigma_beta0[s],
                                     param_samps$sigma_beta1[s],
                                     param_samps$sigma_y[s], Nj_new)
    } # end num_draws loop
  } else if (stanmod_name == "lognormal_binary") {
    for (s in 1:num_draws) {
      beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
      Nj_new <- Nj_new_ln_rng(J, K, Nj_sample, param_samps$mu[s],
                              param_samps$sigma[s])
      Nj_new_df$Nj_new[Nj_new_df$draw_num == s] <- Nj_new
      N_new[s] <- sum(Nj_new)
      ybar_new[s] <- ybar_new_ln_rng(J, K,
                                     beta_samps_sub$beta0,
                                     param_samps$alpha0[s],
                                     param_samps$gamma0[s],
                                     param_samps$sigma_beta0[s],
                                     Nj_new)
    } # end num_draws loop
  } else if (stanmod_name == "negbin") {
    for (s in 1:num_draws) {
      beta_samps_sub <- beta_samps %>%
        dplyr::filter(draw_num == s) %>%
        dplyr::arrange(cluster_id)
      Nj_new <- Nj_new_nb_rng(J, K, Nj_sample, param_samps$mu[s],
                              param_samps$phi[s])
      Nj_new_df$Nj_new[Nj_new_df$draw_num == s] <- Nj_new
      N_new[s] <- sum(Nj_new)
      ybar_new[s] <- ybar_new_nb_rng(J, K, xbar_pop,
                                     beta_samps_sub$beta0,
                                     beta_samps_sub$beta1,
                                     param_samps$alpha0[s],
                                     param_samps$gamma0[s],
                                     param_samps$alpha1[s],
                                     param_samps$gamma1[s],
                                     param_samps$sigma_beta0[s],
                                     param_samps$sigma_beta1[s],
                                     param_samps$sigma_y[s], Nj_new)
    } # end num_draws loop
  } else if (stanmod_name == "negbin_binary") {
    for (s in 1:num_draws) {
      beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
      Nj_new <- Nj_new_nb_rng(J, K, Nj_sample, param_samps$mu[s],
                              param_samps$phi[s])
      Nj_new_df$Nj_new[Nj_new_df$draw_num == s] <- Nj_new
      N_new[s] <- sum(Nj_new)
      ybar_new[s] <- ybar_new_nb_rng(J, K,
                                     beta_samps_sub$beta0,
                                     param_samps$alpha0[s],
                                     param_samps$gamma0[s],
                                     param_samps$sigma_beta0[s],
                                     Nj_new)
    } # end num_draws loop
  } else if (stanmod_name == "negbin_strat" & outcome_type == "continuous") {
    for (s in 1:num_draws) {
      beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
      kappa_samps_sub <- dplyr::filter(kappa_samps, draw_num == s & stratum_id == 2)
      sb_samps_sub <- dplyr::filter(sb_samps, draw_num == s)
      # make strat_inds so we only take the stratum id's for the nonsampled
      # clusters
      strat_inds <- cluster_data_pop$cluster_id %in% c((K+1):J)
      Nj_new <- Nj_new_nb_rng(J, K, S=2,
                              cluster_data_pop$stratum_id[strat_inds],
                              Nj_sample, sb_samps_sub$mu, sb_samps_sub$phi)
      Nj_new_df$Nj_new[Nj_new_df$draw_num == s] <- Nj_new
      N_new[s] <- sum(Nj_new)
      ybar_new[s] <- ybar_new_nb_rng(J, K, S=2, xbar_pop,
                                     beta_samps_sub$beta0,
                                     beta_samps_sub$beta1,
                                     param_samps$alpha0[s],
                                     param_samps$gamma0[s],
                                     param_samps$alpha1[s],
                                     param_samps$gamma1[s],
                                     kappa_samps_sub$kappa0,
                                     kappa_samps_sub$kappa1,
                                     stratum_matrix_pop,
                                     param_samps$sigma_beta0[s],
                                     param_samps$sigma_beta1[s],
                                     param_samps$sigma_y[s], Nj_new)
    } # end num_draws loop
  } else if (stanmod_name == "negbin_ff2" & outcome_type == "continuous") {
    for (s in 1:num_draws) {
      beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
      kappa_samps_sub <- dplyr::filter(kappa_samps, draw_num == s)
      # make strat_inds so we only take the stratum id's for the nonsampled
      # clusters
      strat_inds <- cluster_data_pop$cluster_id %in% c((K+1):J)
      Nj_new <- Nj_new_nb_rng(J, K, Nj_sample,
                              param_samps$mu[s], param_samps$phi[s])
      Nj_new_df$Nj_new[Nj_new_df$draw_num == s] <- Nj_new
      N_new[s] <- sum(Nj_new)
      ybar_new[s] <- ybar_new_nb_rng(J, K, S=9, xbar_pop,
                                     beta_samps_sub$beta0,
                                     beta_samps_sub$beta1,
                                     param_samps$alpha0[s],
                                     param_samps$gamma0[s],
                                     param_samps$alpha1[s],
                                     param_samps$gamma1[s],
                                     kappa_samps_sub$kappa0,
                                     kappa_samps_sub$kappa1,
                                     stratum_matrix_pop,
                                     param_samps$sigma_beta0[s],
                                     param_samps$sigma_beta1[s],
                                     param_samps$sigma_y[s], Nj_new)
    } # end num_draws loop
  } else {
    stop("Invalid stan model")
  }
  print("done drawing Nj_new, ybar_new")
  print(Sys.time())      

yn_list <- list(ybar_new = ybar_new, Nj_new_df = Nj_new_df)
  saveRDS(yn_list,
          paste0(rootdir, "output/simulation/yn_list_usesizes_",
                 use_sizes, "_", outcome_type, "_", size_model, "_", model_name,
                 "_nclusters_", num_clusters, "_nunits_", nunits,
                 "_sim_", simno, ".rds"))

  ##########################################
  ### Plot draws of cluster sizes vs truth
  ##########################################
  if ((grepl("bb", stanmod_name) || grepl("lognormal", stanmod_name) ||
       grepl("negbin", stanmod_name)) && simno == 1) {
    Nj_new_df$in_sample <- Nj_new_df$cluster_id <= K
    tmpdf <- data.frame(cluster_id = c(1:J), draw_num = 9999,
                        Nj_new = Nj_pop, in_sample = FALSE)
    Nj_new_df <- rbind(Nj_new_df, tmpdf)
    Nj_new_df$draw_num <- as.integer(Nj_new_df$draw_num)
    Nj_new_df$is_truth <- ifelse(Nj_new_df$draw_num == 9999, "truth", "draws")
    Nj_sam_df <- data.frame(Nj_sample)
    Nj_pop_df <- data.frame(Nj_pop)
    plt <- ggplot(Nj_new_df, aes(x = Nj_new)) +
      geom_line(aes(group = draw_num, colour = is_truth), stat = "density") +
      geom_line(data = Nj_sam_df, stat = "density",
                aes(x = Nj_sample, colour = "sample")) +
      scale_colour_manual("", values = c("truth" = "black", "draws" = "grey50",
                                         "sample" = "grey80")) +
      scale_x_continuous(limits = c(0, max(Nj_new))) +
      xlab("Cluster size") +
      ylab("Density") +
      ggtitle(paste0("Model: ", stanmod_name, ", ", size_model)) +
      theme_bw()
    ggsave(plt, file = paste0(rootdir, "/output/figures/Nj_draws_usesizes_",
                              use_sizes, "_", outcome_type, "_", size_model, "_",
                              model_name, "_nclusters_", num_clusters,
                              "_nunits_", nunits, "_sim_", simno, ".png"),
           width = 10, height = 8)
  }
  print("done plotting Nj_new")
  print(Sys.time())      

  ##########################################
  ### Plot draws of ybar_new
  ##########################################
  if (simno == 1) {
    tmpdf <- data.frame(ybar_new, truth = ybar_true)
    plt <- ggplot(tmpdf, aes(x = ybar_new)) +
      geom_vline(aes(xintercept = truth), colour = "grey80") +
      geom_line(stat = "density") +
      xlab("Draws of ybar_new") +
      ylab("Density") +
      ggtitle(paste0("Model: ", stanmod_name, ", ", size_model)) +
      theme_bw()
    ggsave(plt, file = paste0(rootdir, "/output/figures/ybar_new_draws_usesizes_",
                              use_sizes, "_", outcome_type, "_", size_model, "_",
                              model_name, "_nclusters_", num_clusters,
                              "_nunits_", nunits, "_sim_", simno, ".png"),
           width = 10, height = 8)
  }
  print("done making plots")
  print(Sys.time())      

  ##########################################
  ### Summarise Nj_new, ybar_new and add to parameter estimates
  ##########################################
  print("summarizing draws")
  if (grepl("cluster_inds_only", stanmod_name) ||
      grepl("knowsizes", stanmod_name)) {
    Nj_new_means <- NA
    tmp <- data.frame(ybar_new, draw_num = c(1:num_draws))
    tmp %>%
      dplyr::summarise(mean = mean(ybar_new),
                       sd = sd(ybar_new),
                       p025 = quantile(ybar_new, 0.025),
                       p25 = quantile(ybar_new, 0.25),
                       p50 = quantile(ybar_new, 0.50),
                       p75 = quantile(ybar_new, 0.75),
                       p975 = quantile(ybar_new, 0.975)) -> draw_summ
    draw_summ$truth <- ybar_true
  } else {
    Nj_new_df %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(Nj_new = mean(Nj_new)) -> Nj_new_means
    Nj_new_means$truth <- Nj_pop

    true_vals <- data.frame(param = c("N_new", "ybar_new"),
                            truth = c(N, ybar_true))
    tmp <- data.frame(N_new, ybar_new, draw_num = c(1:num_draws))
    tmp %>%
      tidyr::gather(key = param, value = value, -draw_num) %>%
      dplyr::group_by(param) %>%
      dplyr::summarise(mean = mean(value),
                       sd = sd(value),
                       p025 = quantile(value, 0.025),
                       p25 = quantile(value, 0.25),
                       p50 = quantile(value, 0.50),
                       p75 = quantile(value, 0.75),
                       p975 = quantile(value, 0.975)) %>%
      dplyr::left_join(., true_vals, by = "param") -> draw_summ
  }
  print(Sys.time())      
  print(draw_summ)
  print(warnings())
  results_list <- list(par_ests = par_ests, Nj_new_means = Nj_new_means,
                       draw_summ = draw_summ, Nj_sample = Nj_sample)
  saveRDS(results_list,
          paste0(rootdir, "output/simulation/stan_results_usesizes_",
                 use_sizes, "_", outcome_type, "_", size_model, "_", model_name,
                 "_nclusters_", num_clusters, "_nunits_", nunits,
                 "_sim_", simno, ".rds"))
  cat("results saved to:", 
      paste0(rootdir, "output/simulation/stan_results_usesizes_",
             use_sizes, "_", outcome_type, "_", size_model, "_", model_name,
             "_nclusters_", num_clusters, "_nunits_", nunits,
             "_sim_", simno, ".rds"), "\n")

  return(NULL)
}
