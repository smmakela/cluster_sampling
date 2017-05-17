# Author: Susanna Makela
# Date: Jan 28, 2016
# Purpose: use the survey package to estimate ybar from the sampled data
svy_ests <- function(J, num.clusters, num.units, use.sizes, outcome.type,
                     rootdir, simno) {
  # Inputs:
  #   J -- number of clusters in pop
  #   num.clusters -- number of sampled clusters
  #   num.units -- number of sampled units
  #   use.sizes -- whether cluster sizes are related to outcome values (0/1)
  #   outcome.type -- whether outcomes are continuous or binary
  #   rootdir -- directory where everything is stored
  #   simno -- which simulation this is

  #############################################################################
  # Useful functions
  #############################################################################
  # y_mat_sum calculates the summand in V_i_hat in equation 8.6.4 (p 314 of
  # Sarndal et al)
    y_mat_sum <- function(i) {
      y_k <- sample.data$y[sample.data$cluster.id == i]
      y_l <- y_k
      pi_k <- PI_k_given_i[i]
      pi_l <- pi_k
      delta_check_kl_given_i <- matrix((pi_kl_i[i] - PI_k_given_i[i]^2) / pi_kl_i[i],
                                       nrow = n.vec$n[i], ncol = n.vec$n[i])
      diag(delta_check_kl_given_i) <- 1 - PI_k_given_i[i]
      res_mat <- delta_check_kl_given_i * outer(y_k / pi_k, y_l / pi_l)
      res <- sum(res_mat)
      return(res)
    }

  # ge_mat_sum calculates the summand in V_CEi_hat in result 8.9.2 for a given
  # cluster i (p 326 of Sarndal et al)
    ge_mat_sum <- function(i) {
      # e_ks defined in equation 8.9.4
      e_ks <- sample.data$y[sample.data$cluster.id == i] -
                sample.data$y_k_hat[sample.data$cluster.id == i]
      e_ls <- e_ks
      # g_ks defined in equation 8.9.12
      # ** the vector of covariates x_k_bold is (1, x), since we have an 
      #    intercept and the covariate x in the data
      # ** t_xi_bold and t_xi_bold_hat both involve sums of x_k_bold within
      #    cluster i. as shown on the bottom of p 323 just below 8.9.8,
      #    t_xi_bold is the sum over ALL units in cluster i, while
      #    t_xi_bold_hat is the sum over all SAMPLED units in cluster i
      sum_x_1 <- cluster.data.pop$Nj_pop[cluster.data.pop$cluster.id <= num.clusters]
      sum_x_2 <- cluster.data.pop$sum_x_i[cluster.data.pop$cluster.id <= num.clusters]
      t_xi_bold <- cbind(sum_x_1, sum_x_2)
      # in t_xi_bold_hat, we can do the division outside of the summation in
      # cluster.data$sum_x_i because PI_k_given_i is constant within clusters
      sum_x_1 <- n.vec$n
      sum_x_2 <- cluster.data$sum_x_i
      t_xi_bold_hat <- cbind(sum_x_1, sum_x_2) / PI_k_given_i
      t_diff <- colSums((t_xi_bold - t_xi_bold_hat) / PI_i)
      # colSums removes the dims from t_diff, so make it into a row vector
      # (we could've just made it have two cols, but since 8.9.12 has
      # t(t_diff), want to keep things consistent with the book)
      t_diff <- matrix(t_diff, ncol = 1) 
      x_k_bold <- rbind(1, sample.data$x[sample.data$cluster.id == i])
      g_ks <- 1 + t(t_diff) %*% solve(T_hat) %*% x_k_bold / sigma2_est # 8.9.12
      g_ks <- as.vector(g_ks) # NEED so that omat has right dimensions
      g_ls <- g_ks
      pi_k <- PI_k_given_i[i]
      pi_l <- pi_k
      delta_check_kl_given_i <- matrix((pi_kl_i[i] - PI_k_given_i[i]^2) / pi_kl_i[i],
                                       nrow = n.vec$n[i], ncol = n.vec$n[i])
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
    if (num.units <= 1) {
      nunits <- paste(num.units*100, "pct", sep = "")
    } else {
      nunits <- num.units
    }
    popdata <- readRDS(paste0(rootdir, "/output/simulation/popdata_usesizes_",
                       use.sizes, "_", outcome.type, ".rds"))
    pop.data <- popdata[["pop.data"]]
    sizetot <- sum(pop.data$Mj)
    ybar.true <- mean(pop.data$y)
    truepars <- popdata[["truepars"]]
    print(truepars)
    rm(popdata)

  # load simno data
    simdata <- readRDS(paste0(rootdir, "output/simulation/simdata_usesizes_",
                              use.sizes, "_", outcome.type, "_nclusters_",
                              num.clusters, "_nunits_", nunits, "_simno_", simno,
                              ".rds"))
    # pull out of list and assign
    for (j in names(simdata)) {
      assign(j, simdata[[j]])
    }
    Nj_sample <- Mj[1:num.clusters]

    # fpc is the number of clusters in the pop
    sample.data$fpc <- J 
    # prob of selecting the cluster
    sample.data$prob <- num.clusters*sample.data$Mj/sizetot 
    # prob of selecting each unit within the cluster
    if (num.units > 1) {
      sample.data$prob2 <- num.units/sample.data$Mj 
    } else {
      sample.data$prob2 <- num.units
    }
    sample.data$wt <- 1/sample.data$prob
    rm(simdata)

  #############################################################################
  # Make vectors/matrices of inclusion probabilities that we'll need later,
  # and other useful quantities
  #############################################################################
  # for the PI's, rename so we can use PI_i etc for the sample values
    PI_i_all <- PI_i
    PI_i <- PI_i[1:num.clusters]
print("PI_i")
print(PI_i)
    PI_ij_all <- PI_ij
    PI_ij <- PI_ij[1:num.clusters, 1:num.clusters]
  # make vectors of unit/joint inclusion probabilities for units in clusters
  # pi_kl_i and pi_kk_i are vectors of length num.clusters    Nj_sample <- Mj[1:num.clusters]
    if (num.units <= 1) { # sample size = num.units * Nj_sample
      PI_k_given_i <- num.units
      PI_k_given_i <- rep(num.units, num.clusters)
      ns <- num.units * Nj_sample # vector of sample sizes
      pi_kl_i <- (ns * (ns - 1)) / (Nj_sample * (Nj_sample - 1))
    } else {
      PI_k_given_i <- num.units / Nj_sample
      pi_kl_i <- num.units * (num.units - 1) / (Nj_sample * (Nj_sample - 1))
    }
    pi_kk_i <- PI_k_given_i

  # covariance of cluster inclusion probabilities
    delta_check_ij <- (PI_ij - outer(PI_i, PI_i)) / PI_ij
    diag(delta_check_ij) <- 1 - PI_i

  # cluster-level data for sample and pop
    cluster.data <- sample.data %>%
      dplyr::group_by(cluster.id) %>%
      dplyr::summarise(sum_y_i = sum(y), sum_x_i = sum(x))
    cluster.data.pop <- pop.data %>%
      dplyr::group_by(cluster.id) %>%
      dplyr::summarise(sum_y_i = sum(y), sum_x_i = sum(x), Nj_pop = mean(Mj))
    sigma2_k <- var(sample.data$y)

  # number of sampled units in each sampled column
    n.vec <- sample.data %>%
      dplyr::group_by(cluster.id) %>%
      dplyr::summarise(n=n())

  #############################################################################
  # HAJEK ESTIMATE -- result 8.6.1, p314 of Sarndal et al
  #############################################################################
  # Sarndal doesn't actually call this the Hajek estimator,
  # but the mean estimate is the same as with the Hajekestimator function... 
  # Mean estimate
    t_yi_star <- cluster.data$sum_y_i / PI_k_given_i
    num <- sum(t_yi_star / PI_i[1:num.clusters])
    den <- sum(Nj_sample / PI_i[1:num.clusters])
    ybar_hat_8.6.1 <- num / den

  # Variance estimate
    d_i <- t_yi_star - Nj_sample * ybar_hat_8.6.1
    d_j <- d_i
    pi_i <- PI_i
    pi_j <- pi_i
    term1 <- sum(delta_check_ij * outer(d_i / pi_i, d_j / pi_j))
    term2 <- sum(sapply(c(1:num.clusters), y_mat_sum) / PI_i)
    denom <- sum(Nj_sample / PI_i)^2
    var_8.6.1 <- (1 / denom) * (term1 + term2)
    ybar_se_8.6.1 <- sqrt(var_8.6.1)

  # Using Hajekestimator function -- variance estimate here is WRONG but
  # I haven't fixed it
  # variance est is same as for HT but with y_i - ybar_hat, according to p 5 of
  # http://jkim.public.iastate.edu/teaching/book8.pdf
    PI_i_rep <- rep(PI_i * PI_k_given_i, times = n.vec$n)
    PI_kl <- matrix(NA, nrow = sum(n.vec$n), ncol = sum(n.vec$n),
                        dimnames = list(rep(c(1:num.clusters), times = n.vec$n),
                                        rep(c(1:num.clusters), times = n.vec$n)))
    for (i in 1:num.clusters) {
      rs <- which(rownames(PI_kl) == i)
      cs <- which(colnames(PI_kl) == i)
      PI_kl[rs, cs] <- pi_kl_i[i]
      diag(PI_kl[rs, cs]) <- pi_kk_i[i]
      for (j in i:num.clusters) {
        rs <- which(rownames(PI_kl) == i)
        cs <- which(colnames(PI_kl) == j)
        PI_kl[rs, cs] <- PI_k_given_i[i] * PI_k_given_i[j]
        rst <- which(rownames(PI_kl) == j)
        cst <- which(colnames(PI_kl) == i)
        PI_kl[rst, cst] <- PI_kl[rs, cs]
      }
    }

    #PI_ij_rep <- matrix(NA, nrow = sum(n.vec$n), ncol = sum(n.vec$n),
    #                    dimnames = list(rep(c(1:num.clusters), times = n.vec$n),
    #                                    rep(c(1:num.clusters), times = n.vec$n)))
    #for (i in 1:num.clusters) {
    #  for (j in i:num.clusters) {
    #    rs <- which(rownames(PI_ij_rep) == i)
    #    cs <- which(colnames(PI_ij_rep) == j)
    #    PI_ij_rep[rs, cs] <- PI_ij[i, j]
    #    rst <- which(rownames(PI_ij_rep) == j)
    #    cst <- which(colnames(PI_ij_rep) == i)
    #    PI_ij_rep[rst, cst] <- PI_ij[i, j]
    #  }
    #}
    ybar_hat_hajek <- Hajekestimator(sample.data$y, PI_i_rep, type = "mean")
    # note this SE is wrong but I'm not sure why -- probably something with the
    # pikl's    
    ybar_se_hajek <- sqrt(varHT(sample.data$y - ybar_hat_hajek,
                                pikl = PI_kl, method = 1))

  # Use the survey package -- this is the old way and gives the same answers
  # as the method from result 8.6.1
    #des <- svydesign(id = ~new.cluster.id, fpc = ~fpc, weights = ~wt, data = sample.data)
    des <- svydesign(id = ~cluster.id+unit.id, fpc = ~prob+prob2,
                     data = sample.data, pps = "brewer")
    #des2 <- svydesign(id = ~cluster.id+unit.id, fpc = ~prob+prob2,
    #                  data = sample.data, pps = ppsmat(PI_ij), variance = "YG")

    # estimate pop mean, pull out std err
    tt <- svymean(~y, des)
    old.ybar_hat_hajek <- as.numeric(tt[1])
    old.ybar_se_hajek <- as.numeric(sqrt(attr(tt,"var")))

  #############################################################################
  # GREG ESTIMATE -- section 8.9, p322 of Sarndal et al
  #############################################################################
  # Equation 8.9.30 -- this is a generalization of Result 8.6.1, which we used
  # above. The variance estimator is still given by 8.9.27, but we use the
  # formula given in 8.9.30
  # NOTE: this only makes sense for the continuous outcome where we use x as a
  # covariate! Otherwise 
  if (outcome.type == "continuous") {
    # Mean estimate
    m1 <- lm(y ~ x, data = sample.data)
    sigma2_est <- summary(m1)$sigma^2 # estimate of sigma^2
    PI_k <- PI_i * PI_k_given_i
    PI_k_rep <- rep(PI_k, times = n.vec$n)
    W <- diag(1 / (sigma2_est * PI_k_rep))
    X <- cbind(rep(1, nrow(sample.data)), sample.data$x)
    T_hat <- t(X) %*% W %*% X
    # basic weighted regression estimator
    Beta_hat <- solve(T_hat) %*% t(X) %*% W %*% sample.data$y 
    y_k_hat <- X %*% Beta_hat
    sample.data$y_k_hat <- as.vector(y_k_hat)
    # 
    # note that in 8.9.8, t_yhat_i is the sum of y_hat for ALL units in the
    # sampled clusters
    t_yhat_i <- pop.data %>%
      dplyr::filter(cluster.id <= num.clusters) %>%
      dplyr::group_by(cluster.id) %>%
      dplyr::summarise(y_k_hat = sum(Beta_hat[1] + Beta_hat[2]*x )) %>%
      dplyr::select(y_k_hat)
    denom <- sum(Nj_sample / PI_i)
    ybar_hat_8.9.8 <- (sum(t_yhat_i / PI_i) +
                       sum((sample.data$y - sample.data$y_k_hat) / PI_k_rep)) /
                       denom

    # Variance estimate
    V_CEi_hat <- sapply(c(1:num.clusters), ge_mat_sum)
    V_CSSU_hat <- sum(V_CEi_hat / PI_i^2)
    t_yi_hat <- cluster.data$sum_y_i / PI_k_given_i
    t_yj_hat <- t_yi_hat
    pi_i <- PI_i
    pi_j <- pi_i
    res_mat <- delta_check_ij * outer(t_yi_hat / pi_i, t_yj_hat / pi_j)
print("res_mat")
print(res_mat)
print("delta_check_ij")
print(delta_check_ij)
print("sum(diag(res_mat))")
print(sum(diag(res_mat)))
print("sum(notdiag(res_mat))")
print(sum(res_mat) - sum(diag(res_mat)))
print("sum((res_mat))")
print(sum(res_mat))
    res <- sum(res_mat)
    V_i_hat <- sapply(c(1:num.clusters), y_mat_sum)
    V_CPSU_hat <- res - sum((1/PI_i) * ((1/PI_i) - 1) * V_i_hat)
    denom <- sum(Nj_sample / PI_i)^2
    var_8.9.8 <- (V_CPSU_hat + V_CSSU_hat) / denom # 8.9.27
    ybar_se_8.9.8 <- sqrt(var_8.9.8)
print("PRINTING VARS:")
print(V_CEi_hat)
print(V_CSSU_hat)
print(V_i_hat)
print(V_CPSU_hat)
print(res)
print(sum((1/PI_i) * ((1/PI_i) - 1) * V_i_hat))
print((1/PI_i) * ((1/PI_i) - 1) * V_i_hat)
  # Using the survey package
    ptot <- sum(pop.data$x)
    pop.totals <- c(`(Intercept)`=nrow(sample.data), x = ptot)
    if (num.units > 1) { # self-weighting
      des2 <- calibrate(des, formula = ~x, pop.totals)
    } else { # make weights be the same by cluster
      des2 <- calibrate(des, formula = ~x, pop.totals, aggregate = 1) 
    }
    tt2 <- svymean(~y, des2)
    old.ybar_hat_greg <- as.numeric(tt2[1])
    old.ybar_se_greg <- as.numeric(sqrt(attr(tt2,"var")))
  } else { # end if outcome.type == "continuous"
    # for binary outcome, no greg estimator
    ybar_hat_8.9.8 <- NA
    ybar_se_8.9.8 <- NA
    old.ybar_hat_greg <- NA
    old.ybar_se_greg <- NA
  }

  #############################################################################
  # Print and save results
  #############################################################################
    print("**************************************************")
    print(paste0("ybar true: ", round(ybar.true, digits = 2)))
    print(paste0("HT estimate (se): ", round(ybar_hat_hajek, digits = 2),
                 " (", round(ybar_se_hajek, digits = 2), ")"))
    print(paste0("HT 8.6.1 estimate (se): ",
                 round(ybar_hat_8.6.1, digits = 2),
                 " (", round(ybar_se_8.6.1, digits = 2), ")"))
    print(paste0("GREG 8.9.8 estimate (se): ",
                 round(ybar_hat_8.9.8, digits = 2),
                 " (", round(ybar_se_8.9.8, digits = 2), ")"))
    print(paste0("OLD HT estimate (se): ", round(old.ybar_hat_hajek, digits = 2),
                 " (", round(old.ybar_se_hajek, digits = 2), ")"))
    print(paste0("OLD GREG estimate (se): ", round(old.ybar_hat_greg, digits = 2),
                 " (", round(old.ybar_se_greg, digits = 2), ")"))
    print("**************************************************")

  # store results so that statistics are wide, model names are long
    res <- data.frame(ybar.true,
                      ybar_hat_hajek = ybar_hat_8.6.1,
                      ybar_se_hajek = ybar_se_8.6.1,
                      ybar_hat_greg = ybar_hat_8.9.8,
                      ybar_se_greg = ybar_se_8.9.8)
    res %>%
      tidyr::gather(key = tmpname, value = value, -ybar.true) %>%
      tidyr::extract(col = tmpname, into = c("stat.name", "model.name"),
                    regex = "(^[^_]+[_][^_]+)(_.*)") %>%
      tidyr::spread(stat.name, value) -> res
    res$model.name <- gsub("\\_", "", res$model.name)
    res$ybar_hat_lci50 <- res$ybar_hat - res$ybar_se
    res$ybar_hat_uci50 <- res$ybar_hat + res$ybar_se
    res$ybar_hat_lci95 <- res$ybar_hat - 1.96*res$ybar_se
    res$ybar_hat_uci95 <- res$ybar_hat + 1.96*res$ybar_se

    saveRDS(res,
            paste0(rootdir, "output/simulation/svy_ests_usesizes_",
                   use.sizes, "_", outcome.type, "_nclusters_", num.clusters,
                   "_nunits_", nunits, "_simno_", simno, ".rds"))
    return(NULL)
}

