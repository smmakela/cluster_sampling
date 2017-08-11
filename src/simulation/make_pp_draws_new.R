make_pp_draws <- function(stanmod, stanmod_name, param_samps, beta_samps, 
                          phi_star_samps = NULL, sb_samps = NULL,
                          pop_cluster_data, sampled_cluster_data, sample_data) {
  # stanmod -- actual compiled stan model
  # stanmod_name -- string containing name of stan model
  # param_samps -- data frame of posterior draws of parameters
  # beta_samps -- data frame of posterior draws of betas
  # phi_star_samps -- data frame of posterior draws of phi_star
  # sb_samps -- data frame of posterior draws of sb_pars
  # pop_cluster_data -- cluster-level data frame of all pop clusters

  J <- nrow(pop_cluster_data)
  K <- sum(pop_cluster_data$is_sampled_cluster)
  pop_cluster_data <- dplyr::arrange(pop_cluster_data, cluster_id)
  sample_data <- dplyr::arrange(pop_cluster_data, cluster_id)

  expose_stan_functions(stanmod)
  num_draws <- nrow(param_samps)
  Mj_new_df <- expand.grid(cluster_id = c(1:J), draw_num = c(1:num_draws))
  Mj_new_df$Mj_new <- NA
  M_tot_new <- rep(NA, times = num_draws)
  ybar_new <- rep(NA, times = num_draws)

  switch(
    bb = {
      M_tot <- sum(pop_cluster_data$Mj)
      n_dat <- sampled_cluster_data %>%
        dplyr::arrange(cluster_id) %>%
        dplyr::group_by(cluster_id, Mj) %>%
        dplyr::summarise(n = n_distinct(cluster_id))
      num_uniq_sz <- nrow(n_dat)      # number of unique cluster sizes
      vec_uniq_sz <- n_dat$Mj # vector of unique cluster sizes
      cts_uniq_sz <- n_dat$n   # counts of unique cluster sizes
      for (s in 1:num_draws) {
        beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
        phi_star_samps_sub <- dplyr::filter(phi_star_samps, draw_num == s)
        Mj_new <- Mj_new_bb_rng(J, K, M_tot, sample_data$Mj, vec_uniq_sx,
                                phi_star_samps_sub$phi_star)
        Mj_new_df$Mj_new[Mj_new_df$draw_num == s] <- Mj_new
        M_tot_new[s] <- sum(Mj_new)
        ybar_new[s] <- ybar_new_bb_rng(J, K, pop_cluster_data$xbar_pop,
                                       beta_samps_sub$beta0,
                                       beta_samps_sub$beta1,
                                       param_samps$alpha0[s],
                                       param_samps$gamma0[s],
                                       param_samps$alpha1[s],
                                       param_samps$gamma1[s],
                                       param_samps$sigma_beta0[s],
                                       param_samps$sigma_beta1[s],
                                       param_samps$sigma_y[s], Mj_new)
      } # end num_draws loop
    },
    bb_binary = {
      M_tot <- sum(pop_cluster_data$Mj)
      n_dat <- sampled_cluster_data %>%
        dplyr::arrange(cluster_id) %>%
        dplyr::group_by(cluster_id, Mj) %>%
        dplyr::summarise(n = n_distinct(cluster_id))
      num_uniq_sz <- nrow(n_dat)      # number of unique cluster sizes
      vec_uniq_sz <- n_dat$Mj # vector of unique cluster sizes
      cts_uniq_sz <- n_dat$n   # counts of unique cluster sizes
      for (s in 1:num_draws) {
        beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
        phi_star_samps_sub <- dplyr::filter(phi_star_samps, draw_num == s)
        Mj_new <- Mj_new_bb_rng(J, K, M, sample_data$Mj, vec_uniq_sz,
                                phi_star_samps_sub$phi_star)
        Mj_new_df$Mj_new[Mj_new_df$draw_num == s] <- Mj_new
        M_tot_new[s] <- sum(Mj_new)
        ybar_new[s] <- ybar_new_bb_rng(J, K,
                                       beta_samps_sub$beta0,
                                       param_samps$alpha0[s],
                                       param_samps$gamma0[s],
                                       param_samps$sigma_beta0[s],
                                       Mj_new)
      } # end num_draws loop
    },
    cluster_inds_only = {
      for (s in 1:num_draws) {
        beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
        ybar_new[s] <- ybar_new_inds_rng(J, K, pop_cluster_data$xbar_pop,
                                         beta_samps_sub$beta0,
                                         beta_samps_sub$beta1,
                                         param_samps$alpha0[s],
                                         param_samps$alpha1[s],
                                         param_samps$sigma_beta0[s],
                                         param_samps$sigma_beta1[s],
                                         param_samps$sigma_y[s],
                                         pop_cluster_data$Nj_pop)
      } # end num_draws loop
    },
    cluster_inds_only_binary = {
      for (s in 1:num_draws) {
        beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
        ybar_new[s] <- ybar_new_inds_rng(J, K,
                                         beta_samps_sub$beta0,
                                         param_samps$alpha0[s],
                                         param_samps$sigma_beta0[s],
                                         pop_cluster_data$Nj_pop)
      } # end num_draws loop
    },
    knowsizes = {
      for (s in 1:num_draws) {
        beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
        ybar_new[s] <- ybar_new_know_rng(J, K, pop_cluster_data$xbar_pop,
                                         beta_samps_sub$beta0,
                                         beta_samps_sub$beta1,
                                         param_samps$alpha0[s],
                                         param_samps$gamma0[s],
                                         param_samps$alpha1[s],
                                         param_samps$gamma1[s],
                                         param_samps$sigma_beta0[s],
                                         param_samps$sigma_beta1[s],
                                         param_samps$sigma_y[s],
                                         pop_cluster_data$Nj_pop)
      } # end num_draws loop
    },
    knowsizes_binary = {
      for (s in 1:num_draws) {
        beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
        ybar_new[s] <- ybar_new_know_rng(J, K,
                                         beta_samps_sub$beta0,
                                         param_samps$alpha0[s],
                                         param_samps$gamma0[s],
                                         param_samps$sigma_beta0[s],
                                         pop_cluster_data$Nj_pop)
      } # end num_draws loop
    },
    lognormal = {
      for (s in 1:num_draws) {
        beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
        Mj_new <- Mj_new_ln_rng(J, K, sample_data$Mj, param_samps$mu[s],
                                param_samps$sigma[s])
        Mj_new_df$Mj_new[Mj_new_df$draw_num == s] <- Mj_new
        M_tot_new[s] <- sum(Mj_new)
        ybar_new[s] <- ybar_new_ln_rng(J, K, pop_cluster_data$xbar_pop,
                                       beta_samps_sub$beta0,
                                       beta_samps_sub$beta1,
                                       param_samps$alpha0[s],
                                       param_samps$gamma0[s],
                                       param_samps$alpha1[s],
                                       param_samps$gamma1[s],
                                       param_samps$sigma_beta0[s],
                                       param_samps$sigma_beta1[s],
                                       param_samps$sigma_y[s], Mj_new)
      } # end num_draws loop
    },
    lognormal_binary = {
      for (s in 1:num_draws) {
        beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
        Mj_new <- Mj_new_ln_rng(J, K, sample_data$Mj, param_samps$mu[s],
                                param_samps$sigma[s])
        Mj_new_df$Mj_new[Mj_new_df$draw_num == s] <- Mj_new
        M_tot_new[s] <- sum(Mj_new)
        ybar_new[s] <- ybar_new_ln_rng(J, K,
                                       beta_samps_sub$beta0,
                                       param_samps$alpha0[s],
                                       param_samps$gamma0[s],
                                       param_samps$sigma_beta0[s],
                                       Mj_new)
      } # end num_draws loop
    },
    negbin = {
      for (s in 1:num_draws) {
        beta_samps_sub <- beta_samps %>%
          dplyr::filter(draw_num == s) %>%
          dplyr::arrange(cluster_id)
        Mj_new <- Mj_new_nb_rng(J, K, sample_data$Mj, param_samps$mu[s],
                                param_samps$phi[s])
        Mj_new_df$Mj_new[Mj_new_df$draw_num == s] <- Mj_new
        M_tot_new[s] <- sum(Mj_new)
        ybar_new[s] <- ybar_new_nb_rng(J, K, pop_cluster_data$xbar_pop,
                                       beta_samps_sub$beta0,
                                       beta_samps_sub$beta1,
                                       param_samps$alpha0[s],
                                       param_samps$gamma0[s],
                                       param_samps$alpha1[s],
                                       param_samps$gamma1[s],
                                       param_samps$sigma_beta0[s],
                                       param_samps$sigma_beta1[s],
                                       param_samps$sigma_y[s], Mj_new)
      } # end num_draws loop
    },
    negbin_binary = {
      for (s in 1:num_draws) {
        beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
        Mj_new <- Mj_new_nb_rng(J, K, sample_data$Mj, param_samps$mu[s],
                                param_samps$phi[s])
        Mj_new_df$Mj_new[Mj_new_df$draw_num == s] <- Mj_new
        M_tot_new[s] <- sum(Mj_new)
        ybar_new[s] <- ybar_new_nb_rng(J, K,
                                       beta_samps_sub$beta0,
                                       param_samps$alpha0[s],
                                       param_samps$gamma0[s],
                                       param_samps$sigma_beta0[s],
                                       Mj_new)
      } # end num_draws loop
    },
    negbin_ff = {
    },
    negbin_ffstrat = {
      for (s in 1:num_draws) {
        beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
        kappa_samps_sub <- dplyr::filter(kappa_samps, draw_num == s & stratum_id == 2)
        sb_samps_sub <- dplyr::filter(sb_samps, draw_num == s)
        # make strat_inds so we only take the stratum id's for the nonsampled
        # clusters
        strat_inds <- pop_cluster_data$cluster_id %in% c((K+1):J)
        Mj_new <- Mj_new_nb_rng(J, K, S=2,
                                pop_cluster_data$stratum_id[strat_inds],
                                sample_data$Mj, sb_samps_sub$mu, sb_samps_sub$phi)
        Mj_new_df$Mj_new[Mj_new_df$draw_num == s] <- Mj_new
        M_tot_new[s] <- sum(Mj_new)
        ybar_new[s] <- ybar_new_nb_rng(J, K, S=2, pop_cluster_data$xbar_pop,
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
                                       param_samps$sigma_y[s], Mj_new)
      } # end num_draws loop
    },
    stop("Invalid stan model")
  ) # close switch()
   # } else if (stanmod_name == "negbin_ff2" & outcome_type == "continuous") {
   #   for (s in 1:num_draws) {
   #     beta_samps_sub <- dplyr::filter(beta_samps, draw_num == s)
   #     kappa_samps_sub <- dplyr::filter(kappa_samps, draw_num == s)
   #     # make strat_inds so we only take the stratum id's for the nonsampled
   #     # clusters
   #     strat_inds <- pop_cluster_data$cluster_id %in% c((K+1):J)
   #     Mj_new <- Mj_new_nb_rng(J, K, sample_data$Mj,
   #                             param_samps$mu[s], param_samps$phi[s])
   #     Mj_new_df$Mj_new[Mj_new_df$draw_num == s] <- Mj_new
   #     M_tot_new[s] <- sum(Mj_new)
   #     ybar_new[s] <- ybar_new_nb_rng(J, K, S=9, pop_cluster_data$xbar_pop,
   #                                    beta_samps_sub$beta0,
   #                                    beta_samps_sub$beta1,
   #                                    param_samps$alpha0[s],
   #                                    param_samps$gamma0[s],
   #                                    param_samps$alpha1[s],
   #                                    param_samps$gamma1[s],
   #                                    kappa_samps_sub$kappa0,
   #                                    kappa_samps_sub$kappa1,
   #                                    stratum_matrix_pop,
   #                                    param_samps$sigma_beta0[s],
   #                                    param_samps$sigma_beta1[s],
   #                                    param_samps$sigma_y[s], Mj_new)
   #   } # end num_draws loop

    to_return <- list(Mj_new_df = Mj_new_df, ybar_new = ybar_new)
} # end function

