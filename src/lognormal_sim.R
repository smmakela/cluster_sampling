negbin_sim <- function(Nj.pop, J, K, stan.model, n.iter, n.chains, n.sims) {
  # Purpose: Simulate cluster sizes from negative binomial, do PPS sampling, and
  #   fit stan model to see how well it fits
  #
  # Author: Susanna Makela
  #
  # Arguments:
  #   Nj.pop -- cluster sizes in simulated population
  #   J -- number of clusters in population
  #   K -- number of clusters sampled
  #   stan.mod -- compiled stan model
  #   n.iter -- number of iterations to run stan for
  #   n.chains -- number of chains to use in stan model
  #   n.sims -- number of posterior predictive simulations to do
  #
  # Returns:
  #   Nj.sims -- data frame of "true" missing cluster sizes and posterior
  #     predictive simulations (true sizes have simno = 0)
  #   summary.df -- data frame of summary statistics of Nj.sims
  #   p.vals -- posterior p-values
  
  
  library(dplyr)
  
  # Generate population data frame
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
                       iter = n.iter, chains = n.chains)
  samps <- extract(stan.fit, permute = FALSE)
  res.summary <- summary(stan.fit)$summary
  cat("###############################################",
      "Stan results ###############################################\n")
  print(round(res.summary, digits = 2))
  
  # Pull out info for true unsampled values
  missing.data <- pop.data[-sample.inds, ]
  Nj.sims <- data.frame(Nj = missing.data$Nj.pop, simno = 0) # simno = 0 for truth
  
  # Calculate test statistics in "true" data
  true.min <- min(missing.data$Nj.pop)
  true.max <- max(missing.data$Nj.pop)
  true.mean <- mean(missing.data$Nj.pop)
  true.median <- median(missing.data$Nj.pop)
  true.q1 <- quantile(missing.data$Nj.pop, probs = 0.25)
  true.q3 <- max(missing.data$Nj.pop, probs = 0.75)
  
  # Do posterior predictive sampling
  k.est <- res.summary["k", "mean"]
  p.est <- res.summary["p", "mean"]
  summary.df <- data.frame()
  for (j in 1:n.sims) {
    Nj.draw <- rnbinom(n = J - K, size = k.est + 1, prob = p.est)
    tmp <- data.frame(Nj = Nj.draw, simno = j)
    Nj.sims <- rbind(Nj.sims, tmp)
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
  
  # Add true values to summary.df
  true.vals <- data.frame(truth = c(true.min, true.max, true.mean,
                                    true.median, true.q1, true.q3),
                          stat = c("min", "max", "mean",
                                   "median", "Q1", "Q3"))
  summary.df <- left_join(summary.df, true.vals, by = "stat")
  
  # Calculate p-values
  p.vals <- dplyr::summarise(group_by(summary.df, stat),
                             p_value = mean(value > truth))
  
  # Return
  toreturn <- list(Nj.sims = Nj.sims, summary.df = summary.df, p.vals = p.vals)
  return(toreturn)
  
}
