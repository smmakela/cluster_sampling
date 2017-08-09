
  expose_stan_functions(stanmod)
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

