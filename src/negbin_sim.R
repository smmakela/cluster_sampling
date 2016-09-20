negbin_sim <- function(J, K, k, p, stan.model, n.iter, n.chains, n.sims) {
  # Purpose: Simulate cluster sizes from negative binomial, do PPS sampling, and
  #   fit stan model to see how well it fits
  #
  # Author: Susanna Makela
  #
  # Arguments:
  #   J -- number of clusters in population
  #   K -- number of clusters to sample
  #   k -- parameter for negative binomial population distribution
  #   p -- parameter for negative binomial population distribution
  #   stan.mod -- compiled stan model
  #   n.iter -- number of iterations to run stan for
  #   n.chains -- number of chains to use in stan model
  #   n.sims -- number of posterior predictive simulations to do
  #
  # Returns:
  #   Nj.pop -- population cluster sizes
  #   Nj.sims -- data frame of "true" missing cluster sizes and posterior
  #     predictive simulations
  #   summary.df -- data frame of summary statistics of Nj.sims
  #   p.vals -- posterior p-values
  
  
  library(dplyr)
  
  # Generate population data
  Nj.pop <- rnbinom(n = J, size = k, prob = p)
  pop.data <- data.frame(ids = c(1:J), Nj.pop)
  
  # Do PPS sampling
  sample.inds <- rspps(Nj.pop, c(1:J), K)
  Nj.sample <- Nj.pop[sample.inds]
  sample.data <- data.frame(sample.ids = c(1:K),
                            orig.ids = sample.inds, Nj.sample)
  
  # Calculate mean, var to get empirical priors for mu, phi
  mean_Nj_sample <- mean(Nj.sample)
  var_Nj_sample <- var(Nj.sample)
  mu_empirical <- mean_Nj_sample
  phi_empirical <- (mu_empirical^2 / (var_Nj_sample - mu_empirical)) - 1
  
  # Run stan model
  stan.data <- list(JminusK = J - K, K = K, Nj_sample = Nj.sample,
                    mu_empirical = mu_empirical,
                    phi_empirical = phi_empirical)
  stan.fit <- sampling(stan.model, data = stan.data,
                       iter = n.iter, chains = n.chans)
  samps <- extract(stan.fit, permute = FALSE)
  res.summary <- summary(stan.fit)$summary
  cat("###################### Stan results ###############################\n")
  cat(round(res.summary, digits = 2))
  
  # Calculate test statistics in "true" data
  true.min <- min(missing.data$Nj.pop)
  true.max <- max(missing.data$Nj.pop)
  true.mean <- mean(missing.data$Nj.pop)
  true.median <- median(missing.data$Nj.pop)
  true.q1 <- quantile(missing.data$Nj.pop, probs = 0.25)
  true.q3 <- max(missing.data$Nj.pop, probs = 0.75)
  summary.df <- data.frame(value = c(true.min, true.max, true.mean,
                                     true.median, true.q1, true.q3),
                           stat = c("min", "max", "mean",
                                    "median", "Q1", "Q3"),
                           simno = 0)
  
  # Do posterior predictive sampling
  k.est <- res.summary["k", "mean"]
  p.est <- res.summary["p", "mean"]
  for (j in 1:n.sims) {
    Nj.draw <- rnbinom(n = J - K, size = k.est + 1, prob = p.est)
    tmp2 <- data.frame(Nj = Nj.draw, simno = j)
    Nj.sims <- rbind(Nj.sims, tmp2)
    sim.min <- min(Nj.draw)
    sim.max <- max(Nj.draw)
    sim.mean <- mean(Nj.draw)
    sim.median <- median(Nj.draw)
    sim.q1 <- quantile(Nj.draw, probs = 0.25)
    sim.q3 <- quantile(Nj.draw, probs = 0.75)
    dd <- data.frame(value = c(sim.min, sim.max, sim.mean,
                               sim.median, sim.q1, sim.q3),
                     stat = c("min", "max", "mean", "median", "Q1", "Q3"),
                     simno = j)
    summary.df <- rbind(summary.df, dd)
  }
  
  for (j in 1:n.sims) {
    tmp <- Nj.sims[Nj.sims$simno == j, ]
    
  }
  summary.df <- cbind(summary.df, truth = c(true.min, true.max, true.mean,
                                            true.median, true.q1, true.q3))
  
  # Calculate p-values
  p.vals <- dplyr::summarise(group_by(summary.df, stat),
                             p_value = mean(stat > truth))
}
