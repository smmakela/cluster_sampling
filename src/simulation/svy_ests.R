# Author: Susanna Makela
# Date: Jan 28, 2016
# Purpose: use the survey package to estimate ybar from the sampled data
svy_ests <- function(J, num.clusters, num.units, use.sizes, outcome.type, rootdir, sim) {
  # Inputs:
  #   J -- number of clusters in pop
  #   num.clusters -- number of sampled clusters
  #   num.units -- number of sampled units
  #   use.sizes -- whether cluster sizes are related to outcome values (0/1)
  #   outcome.type -- whether outcomes are continuous or binary
  #   rootdir -- directory where everything is stored
  #   sim -- which simulation this is

  # load pop data -- use this to get the total sample size so we can calculate the sampling probs (here we are assuming we know M_j in all clusters, so obviously we can calculate the sampling probs
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

  # load sim data
    simdata <- readRDS(paste0(rootdir, "output/simulation/simdata_usesizes_",
                              use.sizes, "_", outcome.type, "_nclusters_",
                              num.clusters, "_nunits_", nunits, "_simno_", sim,
                              ".rds"))
    # pull out of list and assign
    for (j in names(simdata)) {
      assign(j, simdata[[j]])
    }
    # for the PI's, rename so we can use PI_i etc for the sample values
    PI_i_all <- PI_i
    PI_i <- PI_i[1:num.clusters]
    PI_ij_all <- PI_ij
    PI_ij <- PI_ij[1:num.clusters, 1:num.clusters]
  # inclusion and joint inclusion probabilities for units in clusters
  # pi_kl_i and pi_kk_i are vectors of length num.clusters -- the matrix
  # PI_kl_given_i gets made later
    Nj_sample <- Mj[1:num.clusters]
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

    sample.data$fpc <- J # fpc is the number of clusters in the pop
    sample.data$prob <- num.clusters*sample.data$Mj/sizetot # prob of selecting the cluster
    if (num.units > 1) {
      sample.data$prob2 <- num.units/sample.data$Mj # prob of selecting each unit within the cluster
    } else {
      sample.data$prob2 <- num.units
    }
    sample.data$wt <- 1/sample.data$prob
    rm(simdata)


  # Hajek ESTIMATE
    # use sampling package to manually calculate hajek estimator; variance est
    # is same as for HT but with y_i - ybar_hat_as on p 5 of
    # http://jkim.public.iastate.edu/teaching/book8.pdf
    n.vec <- sample.data %>%
      dplyr::group_by(cluster.id) %>%
      dplyr::summarise(n=n())
    PI_i_rep <- rep(PI_i * PI_k_given_i, times = n.vec$n)
    PI_ij_rep <- matrix(NA, nrow = sum(n.vec$n), ncol = sum(n.vec$n),
                        dimnames = list(rep(c(1:num.clusters), times = n.vec$n),
                                        rep(c(1:num.clusters), times = n.vec$n)))
    for (i in 1:num.clusters) {
      for (j in i:num.clusters) {
        rs <- which(rownames(PI_ij_rep) == i)
        cs <- which(colnames(PI_ij_rep) == j)
        PI_ij_rep[rs, cs] <- PI_ij[i, j]
        rst <- which(rownames(PI_ij_rep) == j)
        cst <- which(colnames(PI_ij_rep) == i)
        PI_ij_rep[rst, cst] <- PI_ij[i, j]
      }
    }
    ybar_hat_hajek <- Hajekestimator(sample.data$y, PI_i_rep, type = "mean")    
    ybar_se_hajek <- sqrt(varHT(sample.data$y - ybar_hat_hajek,
                                pikl = PI_ij_rep, method = 1))
  # describe survey design
    #des <- svydesign(id = ~new.cluster.id, fpc = ~fpc, weights = ~wt, data = sample.data)
    des <- svydesign(id = ~cluster.id+unit.id, fpc = ~prob+prob2,
                     data = sample.data, pps = "brewer")
    #des2 <- svydesign(id = ~cluster.id+unit.id, fpc = ~prob+prob2,
    #                  data = sample.data, pps = ppsmat(PI_ij), variance = "YG")

  # estimate pop mean, pull out std err
    tt <- svymean(~y, des)
    old.ybar_hat_hajek <- as.numeric(tt[1])
    old.ybar_se_hajek <- as.numeric(sqrt(attr(tt,"var")))
    #print(svymean(~y, des2))
  # if we got an estimate that's way off from the truth, record the info here
    #if (abs((ybar_hat_- ybar_true)/ybar_true) > 2) {
    #  write.table(sample, file = paste(rootdir, "/Results/Simplify/vary_K/svydesign_usesizes_", use.sizes, "_nclusters_", num.clusters,
    #                         "_nunits_", nunits, "_sim_", sim, ".txt", sep = "")) 
    #  save(des, file = paste(rootdir, "/Results/Simplify/vary_K/svydesign_usesizes_", use.sizes, "_nclusters_", num.clusters,
    #                         "_nunits_", nunits, "_sim_", sim, ".txt", sep = "")) 
    #}

  # Sarndal estimators
  # some functions and constants to make our lives easier
    # covariance of cluster inclusion probabilities
    delta_check_ij <- (PI_ij - outer(PI_i, PI_i)) / PI_ij
    y_mat_sum <- function(i) {
      y_k <- sample.data$y[sample.data$cluster.id == i]
      y_l <- y_k
      pi_k <- PI_k_given_i[i]
      pi_l <- pi_k
      delta_check_kl_given_i <- matrix((pi_kl_i[i] - PI_k_given_i[i]^2) / pi_kl_i[i],
                                       nrow = n.vec$n[i], ncol = n.vec$n[i])
      diag(delta_check_kl_given_i) <- (pi_kk_i[i] - PI_k_given_i[i]^2) / pi_kk_i[i]
      res <- sum(delta_check_kl_given_i * outer(y_k / pi_k, y_l / pi_l))
      return(res)
    }
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
      diag(delta_check_kl_given_i) <- (pi_kk_i[i] - PI_k_given_i[i]^2) / pi_kk_i[i]
      omat <- outer(g_ks * e_ks / pi_k, g_ls * e_ls / pi_l)
      res <- sum(delta_check_kl_given_i * outer(g_ks * e_ks / pi_k, g_ls * e_ls / pi_l))
      return(res)
    }
    cluster.data <- sample.data %>%
      dplyr::group_by(cluster.id) %>%
      dplyr::summarise(sum_y_i = sum(y), sum_x_i = sum(x))
    cluster.data.pop <- pop.data %>%
      dplyr::group_by(cluster.id) %>%
      dplyr::summarise(sum_y_i = sum(y), sum_x_i = sum(x), Nj_pop = mean(Mj))
    sigma2_k <- var(sample.data$y)
  # 8.6.1 -- I think this is actually the Hajek estimator??? they just don't
  #   call it that... but the mean estimate is the same as with the
  #   Hajekestimator function above
    # estimate
    t_yi_star <- cluster.data$sum_y_i / PI_k_given_i
    num <- sum(t_yi_star / PI_i[1:num.clusters])
    den <- sum(Nj_sample / PI_i[1:num.clusters])
    ybar_hat_8.6.1 <- num / den
    # variance
    d_i <- t_yi_star - Nj_sample * ybar_hat_8.6.1
    d_j <- d_i
    pi_i <- PI_i
    pi_j <- pi_i
    term1 <- sum(delta_check_ij * outer(d_i / pi_i, d_j / pi_j))
    term2 <- sum(sapply(c(1:num.clusters), y_mat_sum) / PI_i)
    denom <- sum(Nj_sample / PI_i)^2
    var_8.6.1 <- (1 / denom) * (term1 + term2)
    ybar_se_8.6.1 <- sqrt(var_8.6.1)
    
  # 8.9.8 -- the formulas here are for estimating the total, but we can just
  #   divide the estimator by sum(1 / PI_i) and the variance estimator by
  #   sum(N_i / PI_i)^2 to estimate the mean
  # NOTE: this only makes sense for the continuous outcome where we use x
  if (outcome.type == "continuous") {
    # estimate
    m1 <- lm(y ~ x, data = sample.data)
    sigma2_est <- summary(m1)$sigma^2 # estimate of sigma^2
    PI_k <- PI_i * PI_k_given_i
    PI_k_rep <- rep(PI_k, times = n.vec$n)
    W <- diag(1 / (sigma2_est * PI_k_rep))
    X <- cbind(rep(1, nrow(sample.data)), sample.data$x)
    T_hat <- t(X) %*% W %*% X
    Beta_hat <- solve(T_hat) %*% t(X) %*% W %*% sample.data$y # basic weighted estimator
    y_k_hat <- X %*% Beta_hat
    sample.data$y_k_hat <- as.vector(y_k_hat)
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

    # variance
    V_CEi_hat <- sapply(c(1:num.clusters), ge_mat_sum)
    V_i_hat <- sapply(c(1:num.clusters), y_mat_sum)
    V_CSSU_hat <- sum(V_CEi_hat / PI_i^2)
    t_yi_hat <- cluster.data$sum_y_i / PI_k_given_i
    t_yj_hat <- t_yi_hat
    pi_i <- PI_i
    pi_j <- pi_i
    res <- sum(delta_check_ij * outer(t_yi_hat / pi_i, t_yj_hat / pi_j))
    V_CPSU_hat <- res - sum((1/PI_i) * ((1/PI_i) - 1) * V_i_hat)
    denom <- sum(Nj_sample / PI_i)^2
    var_8.9.8 <- (V_CPSU_hat + V_CSSU_hat) / denom # 8.9.27
    ybar_se_8.9.8 <- sqrt(var_8.9.8)

  # old way
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
    ybar_hat_8.9.8 <- NA
    ybar_se_8.9.8 <- NA
    old.ybar_hat_greg <- NA
    old.ybar_se_greg <- NA
  }
  # if we got an estimate that's way off from the truth, record the info here
    #d1 <- abs((ybar_hat_- ybar_true)/ybar_true)
    #d2 <- abs((ybar_hat_ - ybar_true)/ybar_true)
    #if ((d1 > 2) | (d2 > 2)) {
    #  write.csv(sample.data, file = paste(rootdir, "/Results/Simplify/vary_K/sampdat_weird_svy_usesizes_", use.sizes,
    #                                      "_nclusters_", num.clusters, "_nunits_", nunits, "_sim_", sim, ".txt", sep = ""))
    #}

  # print
    print("**************************************************")
    print(paste0("ybar true: ", round(ybar.true, digits = 2)))
    print(paste0("HT estimate (se): ", round(ybar_hat_hajek, digits = 2),
                 " (", round(ybar_se_hajek, digits = 2), ")"))
    print(paste0("GREG 8.6.1 estimate (se): ",
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
    #res <- data.frame(ybar.true,
    #                  ybar_hat_hajek = old.ybar_hat_hajek,
    #                  ybar_se_hajek = old.ybar_se_hajek,
    #                  ybar_hat_greg = old.ybar_hat_greg,
    #                  ybar_se_greg = old.ybar_se_greg)
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
    res$ybar_se <- NULL
  saveRDS(results.list,
          paste0(rootdir, "output/simulation/results_usesizes_",
                 use.sizes, "_", outcome.type, "_", sn, 
                 "_nclusters_", num.clusters,
                 "_nunits_", nunits, "_sim_", simno, ".rds"))
    return(res)
}
