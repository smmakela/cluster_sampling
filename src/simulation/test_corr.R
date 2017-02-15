library(rstan)
library(ggplot2)
library(tidyr)
library(dplyr)
options(mc.cores = parallel::detectCores())
corr.model <- stan_model("/Users/susanna/projects/cluster_sampling/src/simulation/knowsizes_corr.stan")
params <- c("alpha0", "gamma0", "alpha1", "gamma1",
            "sigma_beta0", "sigma_beta1", "sigma_y", "ybar")
  J_pop <- 1000
  J_sam <- 30
  Mj_pop <- sample(c(10:500), size = J_pop, replace = TRUE)
  logMj_pop <- log(Mj_pop) - mean(log(Mj_pop))
  inds <- rep(c(1:J_pop), Mj_pop)
  N_pop <- sum(Mj_pop)
  alpha0 <- rnorm(1, 0, 1)
  gamma0 <- rnorm(1, 0, 1)
  alpha1 <- rnorm(1, 0, 1)
  gamma1 <- rnorm(1, 0, 1)
  sigma_beta0 <- abs(rnorm(1, 0, 0.5))
  sigma_beta1 <- abs(rnorm(1, 0, 0.5))
  beta0 <- rnorm(J_pop, mean = alpha0 + gamma0*logMj_pop, sigma_beta0)
  beta1 <- rnorm(J_pop, mean = alpha1 + gamma1*logMj_pop, sigma_beta1)
  sigma_y <- abs(rnorm(1, 0, 0.5))
  x <- runif(N_pop, min = 15, max = 45)
  x <- x - mean(x)
  y <- rep(0, times = N_pop)
  for (i in 1:N_pop) {
    y[i] <- rnorm(1, mean = beta0[inds[i]] + beta1[inds[i]] * x[i], sd = sigma_y)
  }
  sampled.clusters <- sample.int(J_pop, size = J_sam, prob = Mj_pop, replace = FALSE)
  Mj_sam <- Mj_pop[sampled.clusters]
  logMj_sam <- logMj_pop[sampled.clusters]
  df <- data.frame(inds, x, y)
  sample.data <- dplyr::filter(df, inds %in% sampled.clusters)
  sample.data$cluster_id_long <- as.integer(as.factor(sample.data$inds))
  N_sam <- nrow(sample.data)
  standata <- list(J_pop = J_pop, J_sam = J_sam, N_sam = N_sam,
                   cluster_id_long = sample.data$cluster_id_long,
                   y = sample.data$y, x = cbind(rep(1, nrow(sample.data)), sample.data$x),
                   Mj_sam = Mj_sam, logMj_sam = logMj_sam,
                   u = cbind(rep(1, length(logMj_sam)), logMj_sam))
  stan.fit <- sampling(corr.model, data = standata, iter = 1000, chains = 4)
  