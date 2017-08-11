# Author: Susanna Makela
# Date: Jan 28, 2016
# Purpose: use the survey package to estimate ybar from the sampled data
svy_ests <- function(J, num_clusters, num_units, use_sizes, outcome_type,
                     size_model, sim_data, rootdir, simno) {
  # Inputs:
  #   J -- number of clusters in pop
  #   num_clusters -- number of sampled clusters
  #   num_units -- number of sampled units
  #   use_sizes -- whether cluster sizes are related to outcome values (0/1)
  #   outcome_type -- whether outcomes are continuous or binary
  #   size_model -- model used to generate pop cluster sizes
  #   sim_data -- list output from sampledata function, has pop and sampled data
  #   rootdir -- directory where everything is stored
  #   simno -- which simulation this is

  #############################################################################
  # Useful functions
  #############################################################################
  # y_mat_sum calculates the summand in V_i_hat in equation 8_6_4 (p 314 of
  # Sarndal et al)
    y_mat_sum <- function(i) {
      samdat <- sample_data[sample_data$cluster_id == i, ]
      y_k_pi_k <- samdat$y / samdat$PI_k_given_i
      y_l_pi_l <- y_k_pi_k
      #y_k <- sample_data$y[sample_data$cluster_id == i]
      #y_l <- y_k
      #pi_k <- PI_k_given_i[i]
      #pi_l <- pi_k
      delta_check_kl_given_i <- matrix((pi_kl_i[i] - PI_k_given_i[i]^2) / pi_kl_i[i],
                                       nrow = n_vec$n[i], ncol = n_vec$n[i])
      diag(delta_check_kl_given_i) <- 1 - PI_k_given_i[i]
      res_mat <- delta_check_kl_given_i * outer(y_k_pi_k, y_l_pi_l)
      #res_mat <- delta_check_kl_given_i * outer(y_k / pi_k, y_l / pi_l)
      res <- sum(res_mat)
      return(res)
    }

  # ge_mat_sum calculates the summand in V_CEi_hat in result 8_9_2 for a given
  # cluster i (p 326 of Sarndal et al)
    ge_mat_sum <- function(i) {
      # e_ks defined in equation 8_9_4
      e_ks <- sample_data$y[sample_data$cluster_id == i] -
                sample_data$y_k_hat[sample_data$cluster_id == i]
      e_ls <- e_ks
      # g_ks defined in equation 8_9_12
      # ** the vector of covariates x_k_bold is (1, x), since we have an 
      #    intercept and the covariate x in the data
      # ** t_xi_bold and t_xi_bold_hat both involve sums of x_k_bold within
      #    cluster i_ as shown on the bottom of p 323 just below 8_9_30,
      #    t_xi_bold is the sum over ALL units in cluster i, while
      #    t_xi_bold_hat is the sum over all SAMPLED units in cluster i
      sum_x_1 <- cluster_data_pop$Nj_pop[cluster_data_pop$cluster_id <= num_clusters]
      sum_x_2 <- cluster_data_pop$sum_x_i[cluster_data_pop$cluster_id <= num_clusters]
      t_xi_bold <- cbind(sum_x_1, sum_x_2)
      # in t_xi_bold_hat, we can do the division outside of the summation in
      # cluster_data$sum_x_i because PI_k_given_i is constant within clusters
      sum_x_1 <- n_vec$n
      sum_x_2 <- cluster_data$sum_x_i
      t_xi_bold_hat <- cbind(sum_x_1, sum_x_2) / PI_k_given_i
      t_diff <- colSums((t_xi_bold - t_xi_bold_hat) / PI_i)
      # colSums removes the dims from t_diff, so make it into a row vector
      # (we could've just made it have two cols, but since 8_9_12 has
      # t(t_diff), want to keep things consistent with the book)
      t_diff <- matrix(t_diff, ncol = 1) 
      x_k_bold <- rbind(1, sample_data$x[sample_data$cluster_id == i])
      g_ks <- 1 + t(t_diff) %*% solve(T_hat) %*% x_k_bold / sigma2_est # 8_9_12
      g_ks <- as.vector(g_ks) # NEED so that omat has right dimensions
      g_ls <- g_ks
      pi_k <- PI_k_given_i[i]
      pi_l <- pi_k
      delta_check_kl_given_i <- matrix((pi_kl_i[i] - PI_k_given_i[i]^2) / pi_kl_i[i],
                                       nrow = n_vec$n[i], ncol = n_vec$n[i])
      #diag(delta_check_kl_given_i) <- (pi_kk_i[i] - PI_k_given_i[i]^2) / pi_kk_i[i]
      diag(delta_check_kl_given_i) <- 1 - PI_k_given_i[i]
      omat <- outer(g_ks * e_ks / pi_k, g_ls * e_ls / pi_l)
      res <- sum(delta_check_kl_given_i * outer(g_ks * e_ks / pi_k, g_ls * e_ls / pi_l))
      return(res)
    }

  #############################################################################
  # Load pop and simno data
  #############################################################################
  # load pop data -- use this to get the total sample size so we can calculate
  # the sampling probs (here we are assuming we know M_j in all clusters, so
  # obviously we can calculate the sampling probs
    if (num_units <= 1) {
      nunits <- paste(num_units*100, "pct", sep = "")
    } else {
      nunits <- num_units
    }
    # pull everything out of sim_data and assign to its name
    for (j in names(sim_data)) {
      assign(j, sim_data[[j]])
    }

    # Sort things by cluster
    pop_cluster_data <- dplyr::arrange(pop_cluster_data, cluster_id)
    sampled_cluster_data <- dplyr::arrange(sampled_cluster_data, cluster_id)
    sample_data <- dplyr::arrange(sample_data, cluster_id)
 
    sizetot <- sum(pop_cluster_data$Mj)
    ybar_true <- mean(pop_data$y)
    rm(sim_data)

    Mj_sample <- sampled_cluster_data$Mj
    Mj_pop <- pop_cluster_data$Mj

    # fpc is the number of clusters in the pop
    sample_data$fpc <- J 
    # prob of selecting the cluster
    sample_data$prob <- num_clusters*Mj_sample/sizetot 
    # prob of selecting each unit within the cluster
    if (num_units > 1) {
      sample_data$prob2 <- num_units/Mj_sample 
    } else {
      sample_data$prob2 <- num_units
    }
    sample_data$wt <- 1/sample_data$prob

  #############################################################################
  # Make vectors/matrices of inclusion probabilities that we'll need later,
  # and other useful quantities
  #############################################################################
  # for the PI's, rename so we can use PI_i etc for the sample values
    PI_i_all <- PI_i
    PI_i <- PI_i[1:num_clusters]
    PI_ij_all <- PI_ij
    PI_ij <- PI_ij[1:num_clusters, 1:num_clusters]
  # make vectors of unit/joint inclusion probabilities for units in clusters
  # pi_kl_i and pi_kk_i are vectors of length num_clusters    Mj_sample <- Mj[1:num_clusters]
    if (num_units <= 1) { # sample size = num_units * Mj_sample
      PI_k_given_i <- num_units
      PI_k_given_i <- rep(num_units, num_clusters)
      ns <- num_units * Mj_sample # vector of sample sizes
      pi_kl_i <- (ns * (ns - 1)) / (Mj_sample * (Mj_sample - 1))
    } else {
      PI_k_given_i <- num_units / Mj_sample
      pi_kl_i <- num_units * (num_units - 1) / (Mj_sample * (Mj_sample - 1))
    }
    pi_kk_i <- PI_k_given_i

  # covariance of cluster inclusion probabilities
    delta_check_ij <- (PI_ij - outer(PI_i, PI_i)) / PI_ij
    diag(delta_check_ij) <- 1 - PI_i

  # cluster-level data for sample and pop
    cluster_data <- sample_data %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(sum_y_i = sum(y), sum_x_i = sum(x))
    cluster_data_pop <- pop_data %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(sum_y_i = sum(y), sum_x_i = sum(x), Nj_pop = mean(Mj))
    sigma2_k <- var(sample_data$y)

  # number of sampled units in each sampled column
    n_vec <- sample_data %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(n=n())

  # add PI_i, PI_k_given_i, PI_k to sample_data
    pi_dat <- data_frame(PI_i, PI_k_given_i, cluster_id = c(1:num_clusters))
    sample_data <- left_join(sample_data, pi_dat, by = "cluster_id")
    sample_data$PI_k <- sample_data$PI_i * sample_data$PI_k_given_i

  #############################################################################
  # HAJEK ESTIMATE -- result 8_6_1, p314 of Sarndal et al
  #############################################################################
  # Sarndal doesn't actually call this the Hajek estimator,
  # but the mean estimate is the same as with the Hajekestimator function___ 
  # Mean estimate
    t_yi_star <- cluster_data$sum_y_i / PI_k_given_i
    num <- sum(t_yi_star / PI_i[1:num_clusters])
    den <- sum(Mj_sample / PI_i[1:num_clusters])
    ybar_hat_8_6_1 <- num / den

  # Variance estimate
    d_i <- t_yi_star - Mj_sample * ybar_hat_8_6_1
    d_j <- d_i
    pi_i <- PI_i
    pi_j <- pi_i
    term1 <- sum(delta_check_ij * outer(d_i / pi_i, d_j / pi_j))
    term2 <- sum(sapply(c(1:num_clusters), y_mat_sum) / PI_i)
    denom <- sum(Mj_sample / PI_i)^2
    var_8_6_1 <- (1 / denom) * (term1 + term2)
    ybar_se_8_6_1 <- sqrt(var_8_6_1)

  # Using Hajekestimator function
  # variance est is same as for HT but with y_i - ybar_hat, according to p 5 of
  # http://jkim_public_iastate_edu/teaching/book8_pdf
    PI_kl <- matrix(NA, nrow = sum(n_vec$n), ncol = sum(n_vec$n),
                        dimnames = list(rep(c(1:num_clusters), times = n_vec$n),
                                        rep(c(1:num_clusters), times = n_vec$n)))
    for (i in 1:num_clusters) {
      for (j in i:num_clusters) {
        if (i == j) {
          rs <- which(rownames(PI_kl) == i)
          cs <- which(colnames(PI_kl) == i)
          PI_kl[rs, cs] <- pi_kl_i[i]
          diag(PI_kl[rs, cs]) <- pi_kk_i[i]
        } else {
          rs <- which(rownames(PI_kl) == i)
          cs <- which(colnames(PI_kl) == j)
          PI_kl[rs, cs] <- PI_ij[i,j] * PI_k_given_i[i] * PI_k_given_i[j]
          rst <- which(rownames(PI_kl) == j)
          cst <- which(colnames(PI_kl) == i)
          PI_kl[rst, cst] <- PI_kl[rs, cs]
        }
      }
    }

    ybar_hat_hajek <- Hajekestimator(sample_data$y, pik = sample_data$PI_k,
                                     type = "mean")
    # note this SE is wrong but I'm not sure why -- probably something with the
    # pikl's   
    # add sum(sample_data$PI_k) in front because is like multiplying by
    # 1 / N_hat, where N_hat = 1 / sum(pi_k) for all units k in the sample 
    ybar_se_hajek <- 1/(sum(1/sample_data$PI_k)) *
                     sqrt(varHT(sample_data$y - ybar_hat_hajek,
                                pikl = PI_kl, method = 1))
  # Use the survey package -- this is the old way and gives the same answers
  # as the method from result 8_6_1
    #des <- svydesign(id = ~new_cluster_id, fpc = ~fpc, weights = ~wt, data = sample_data)
    des <- svydesign(id = ~cluster_id+unit_id, fpc = ~prob+prob2,
                     data = sample_data, pps = "brewer")
    #des2 <- svydesign(id = ~cluster_id+unit_id, fpc = ~prob+prob2,
    #                  data = sample_data, pps = ppsmat(PI_ij), variance = "YG")

    # estimate pop mean, pull out std err
    tt <- svymean(~y, des)
    old_ybar_hat_hajek <- as.numeric(tt[1])
    old_ybar_se_hajek <- as.numeric(sqrt(attr(tt, "var")))

  #############################################################################
  # GREG ESTIMATE -- section 8_9, p322 of Sarndal et al
  #############################################################################
  # Equation 8_9_30 -- this is a generalization of Result 8_6_1, which we used
  # above_ The variance estimator is still given by 8_9_27, but we use the
  # formula given in 8_9_30
  # NOTE: this only makes sense for the continuous outcome where we use x as a
  # covariate! Otherwise 
  if (outcome_type == "continuous") {
    # Mean estimate
    m1 <- lm(y ~ x, data = sample_data)
    sigma2_est <- summary(m1)$sigma^2 # estimate of sigma^2
    PI_k <- PI_i * PI_k_given_i
    PI_k_rep <- rep(PI_k, times = n_vec$n)
    W <- diag(1 / (sigma2_est * PI_k_rep))
    X <- cbind(rep(1, nrow(sample_data)), sample_data$x)
    T_hat <- t(X) %*% W %*% X
    # basic weighted regression estimator
    Beta_hat <- solve(T_hat) %*% t(X) %*% W %*% sample_data$y 
    y_k_hat <- X %*% Beta_hat
    sample_data$y_k_hat <- as.vector(y_k_hat)
    # make t_yhat_ir as in 8_9_6
    # first term of 8_9_6 -- predictions for ALL units **in the i-th cluster**
    t1 <- pop_data %>%
      dplyr::filter(cluster_id <= num_clusters) %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(y_k_hat = sum(Beta_hat[1] + Beta_hat[2]*x)) %>%
      dplyr::select(y_k_hat)
    # second term of 8_9_6 -- weighted residuals for sampled units
    t2 <- sample_data %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(wt_diff = sum((y - y_k_hat) / PI_k_given_i)) %>%
      dplyr::select(wt_diff)
    # add to get t_yhat_ir
    t_yhat_ir <- t1 + t2
    # now calculate 8_9_30
    ybar_hat_8_9_30 <- sum(t_yhat_ir / PI_i) / sum(Mj_sample / PI_i)
    # note that in 8_9_30, t_yhat_i is the sum of y_hat for ALL units in the
    # sampled clusters
    t_yhat_i <- pop_data %>%
      dplyr::filter(cluster_id <= num_clusters) %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(y_k_hat = sum(Beta_hat[1] + Beta_hat[2]*x )) %>%
      dplyr::select(y_k_hat)
    denom <- sum(Mj_sample / PI_i)
    ybar_hat_8_9_8 <- (sum(t_yhat_i / PI_i) +
                       sum((sample_data$y - sample_data$y_k_hat) / sample_data$PI_k_rep)) /
                       denom

    # Variance estimate
    V_CEi_hat <- sapply(c(1:num_clusters), ge_mat_sum)
    V_CSSU_hat <- sum(V_CEi_hat / PI_i^2)
    t_yi_hat <- cluster_data$sum_y_i / PI_k_given_i
    t_yj_hat <- t_yi_hat
    pi_i <- PI_i
    pi_j <- pi_i
    res_mat <- delta_check_ij * outer(t_yi_hat / pi_i, t_yj_hat / pi_j)
    res <- sum(res_mat)
    V_i_hat <- sapply(c(1:num_clusters), y_mat_sum)
    V_CPSU_hat <- res - sum((1/PI_i) * ((1/PI_i) - 1) * V_i_hat)
    denom <- sum(Mj_sample / PI_i)^2
    var_8_9_30 <- (V_CPSU_hat + V_CSSU_hat) / denom # 8_9_27
    ybar_se_8_9_30 <- sqrt(var_8_9_30)
  # Using the survey package
    ptot <- sum(pop_data$x)
    pop_totals <- c(`(Intercept)`=nrow(sample_data), x = ptot)
    if (num_units > 1) { # self-weighting
      des2 <- calibrate(des, formula = ~x, pop_totals)
    } else { # make weights be the same by cluster
      des2 <- calibrate(des, formula = ~x, pop_totals, aggregate = 1) 
    }
    tt2 <- svymean(~y, des2)
    old_ybar_hat_greg <- as.numeric(tt2[1])
    old_ybar_se_greg <- as.numeric(sqrt(attr(tt2,"var")))
  } else { # end if outcome_type == "continuous"
    # for binary outcome, no greg estimator
    ybar_hat_8_9_30 <- NA
    ybar_se_8_9_30 <- NA
    ybar_hat_8_9_8 <- NA
    old_ybar_hat_greg <- NA
    old_ybar_se_greg <- NA
  }

  #############################################################################
  # Print and save results
  #############################################################################
    print("**************************************************")
    print(paste0("ybar true: ", round(ybar_true, digits = 2)))
    print(paste0("Hajek estimate (se): ", round(ybar_hat_hajek, digits = 2),
                 " (", round(ybar_se_hajek, digits = 2), ")"))
    print(paste0("Hajek 8_6_1 estimate (se): ",
                 round(ybar_hat_8_6_1, digits = 2),
                 " (", round(ybar_se_8_6_1, digits = 2), ")"))
    print(paste0("GREG 8_9_30 estimate (se): ",
                 round(ybar_hat_8_9_30, digits = 2),
                 " (", round(ybar_se_8_9_30, digits = 2), ")"))
    print(paste0("GREG 8_9_8 estimate: ",
                 round(ybar_hat_8_9_8, digits = 2)))
    print(paste0("OLD Hajek estimate (se): ", round(old_ybar_hat_hajek, digits = 2),
                 " (", round(old_ybar_se_hajek, digits = 2), ")"))
    print(paste0("OLD GREG estimate (se): ", round(old_ybar_hat_greg, digits = 2),
                 " (", round(old_ybar_se_greg, digits = 2), ")"))
    print("**************************************************")

  # store results so that statistics are wide, model names are long
    res <- data_frame(ybar_true,
                      ybar_hat_hajek = ybar_hat_8_6_1,
                      ybar_se_hajek = ybar_se_8_6_1,
                      ybar_hat_greg = ybar_hat_8_9_30,
                      ybar_se_greg = ybar_se_8_9_30)
    res %>%
      tidyr::gather(key = tmpname, value = value, -ybar_true) %>%
      tidyr::extract(col = tmpname, into = c("stat_name", "model_name"),
                    regex = "^([^_]*_[^_]*)_(.*)$") %>%
      tidyr::spread(stat_name, value) -> res
    res$model_name <- gsub("\\_", "", res$model_name)
    res$ybar_hat_lci50 <- res$ybar_hat - res$ybar_se
    res$ybar_hat_uci50 <- res$ybar_hat + res$ybar_se
    res$ybar_hat_lci95 <- res$ybar_hat - 1.96*res$ybar_se
    res$ybar_hat_uci95 <- res$ybar_hat + 1.96*res$ybar_se

    saveRDS(res,
            paste0(rootdir, "output/simulation/svy_ests_usesizes_",
                   use_sizes, "_", outcome_type, "_", size_model, "_nclusters_",
                   num_clusters, "_nunits_", nunits, "_simno_", simno, ".rds"))
    return(NULL)
}

