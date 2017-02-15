N <- 1000
alpha <- rnorm(1, 0, 1)
beta <- rnorm(1, 0, 1)
sigma <- 1
x <- runif(N, min = 15, max = 45)
x <- x - mean(x)
y <- rnorm(N, mean = alpha + beta * x, sd = sigma)
standata <- list(N = N, y = y, x = x)
stan.fit <- stan(file = "/Users/susanna/projects/cluster_sampling/src/simulation/simple.stan",
                 data = standata, iter = 1000, chains = 4)
print(stan.fit)
cat("True values:\n", "alpha =", alpha, "\n", "beta =", beta, "\n", "sigma =", sigma, "\n")


library(rstan)
library(ggplot2)
library(tidyr)
library(dplyr)
simple3.model <- stan_model("/Users/susanna/projects/cluster_sampling/src/simulation/simple3.stan")
simple3.nomatt.model <- stan_model("/Users/susanna/projects/cluster_sampling/src/simulation/simple3_nomatt.stan")
n.sims <- 1
all.res <- data.frame()
beta.res <- data.frame()
params <- c("alpha0", "gamma0", "alpha1", "gamma1",
            "sigma_beta0", "sigma_beta1", "sigma_y", "ybar")
for (s in 1:n.sims) {
  cat("Currently on sim", s, "of", n.sims, "sims\n")
  simno <- s
  N <- 1000
  J <- 30
  inds <- sample.int(n = J, size = N, replace = TRUE, prob = c(1:J))
  while(length(unique(inds)) != 30) {
    inds <- sample.int(n = J, size = N, replace = TRUE, prob = c(1:J))
  }
  inds <- sort(inds)
  tmp <- data.frame(table(inds))
  tmp$inds <- as.numeric(tmp$inds)
  Mj <- rep(0, J)
  Mj[tmp$inds] <- tmp$Freq
  logMj <- log(Mj) - mean(log(Mj))
  alpha0 <- rnorm(1, 0, 1)
  gamma0 <- rnorm(1, 0, 1)
  alpha1 <- rnorm(1, 0, 1)
  gamma1 <- rnorm(1, 0, 1)
  sigma_beta0 <- abs(rnorm(1, 0, 0.5))
  sigma_beta1 <- abs(rnorm(1, 0, 0.5))
  beta0 <- rnorm(J, mean = alpha0 + gamma0*logMj, sigma_beta0)
  beta1 <- rnorm(J, mean = alpha1 + gamma1*logMj, sigma_beta1)
  sigma_y <- abs(rnorm(1, 0, 0.5))
  x <- runif(N, min = 15, max = 45)
  x <- x - mean(x)
  y <- rep(0, times = N)
  for (i in 1:N) {
    y[i] <- rnorm(1, mean = beta0[inds[i]] + beta1[inds[i]] * x[i], sd = sigma_y)
  }
  truth <- c(alpha0, gamma0, alpha1, gamma1, sigma_beta0, sigma_beta1, sigma_y, mean(y))
  tmp <- data.frame(inds, x)
  tmp <- dplyr::summarise(group_by(tmp, inds), xbar = mean(x))
  xbar <- rep(0, J)
  xbar[tmp$inds] <- tmp$xbar
  standata <- list(N = N, J = J, inds = inds, Mj = Mj, logMj = logMj, y = y, x = x, xbar = xbar)
  stan.fit <- sampling(simple3.model, data = standata,
                       iter = 1000, chains = 4, refresh = -1)
  stan.fit.nomatt <- sampling(simple3.nomatt.model, data = standata,
                              iter = 1000, chains = 4, refresh = -1)
  tmp <- summary(stan.fit)
  tmp.nomatt <- summary(stan.fit.nomatt)
  #rm(stan.fit)
  #rm(stan.fit.nomatt)
  tmp <- tmp$summary
  tmp.nomatt <- tmp.nomatt$summary
  posterior.mean <- c(tmp[params, "mean"], tmp.nomatt[params, "mean"])
  p025 <- c(tmp[params, "2.5%"], tmp.nomatt[params, "2.5%"])
  p975 <- c(tmp[params, "97.5%"], tmp.nomatt[params, "97.5%"])
  p25 <- c(tmp[params, "25%"], tmp.nomatt[params, "25%"])
  p75 <- c(tmp[params, "75%"], tmp.nomatt[params, "75%"])
  p50 <- c(tmp[params, "50%"], tmp.nomatt[params, "50%"])
  n.eff <- c(tmp[params, "n_eff"], tmp.nomatt[params, "n_eff"])
  Rhat <- c(tmp[params, "Rhat"], tmp.nomatt[params, "Rhat"])
  whichmod <- rep(c("matt", "nomatt"), each = length(params))
  #beta0.est <- tmp[grep("^beta0", rownames(tmp)), "mean"]
  #beta1.est <- tmp[grep("^beta1", rownames(tmp)), "mean"]
  #curr.betas <- data.frame(simno = s, ind = c(1:J), Mj, xbar,
  #                         beta0, beta1, beta0.est, beta1.est)
  #beta.res <- rbind(beta.res, curr.betas)
  curr.res <- data.frame(simno, params, truth, posterior.mean, p025, p975,
                         p25, p75, p50, n.eff, Rhat, whichmod)
  all.res <- rbind(all.res, curr.res)
  #print(stan.fit.nomatt, pars = params)
  #cat("True values:\n", "alpha0 =", round(alpha0, digits = 2), "\n",
  #    "gamma0 =", round(gamma0, digits = 2), "\n",
  #    "alpha1 =", round(alpha1, digits = 2), "\n",
  #    "gamma1 =", round(gamma1, digits = 2), "\n",
  #    "sigma_beta0 =", round(sigma_beta0, digits = 2), "\n",
  #    "sigma_beta1 =", round(sigma_beta1, digits = 2), "\n",
  #    "sigma_y =", round(sigma_y, digits = 2), "\n",
  #    "ybar =", round(mean(y), digits = 2))
}

ggplot(all.res, aes(x = n.eff)) +
  geom_line(stat = "density", aes(colour = whichmod)) +
  facet_wrap(~params)

ggplot(all.res, aes(x = log(Rhat))) +
  geom_line(stat = "density", aes(colour = whichmod)) +
  facet_wrap(~params)



ggplot(all.res, aes(x = posterior.mean, y = truth)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~params, scales = "free") +
  theme_bw()

beta.res$yj.true <- beta.res$beta0 + beta.res$beta1 * beta.res$xbar
beta.res$yj.est <- beta.res$beta0.est + beta.res$beta1.est * beta.res$xbar
beta.res$yj.diff <- beta.res$yj.true - beta.res$yj.est
beta.res$yj.rel.diff <- beta.res$yj.diff/beta.res$yj.true
ggplot(beta.res, aes(x = yj.est, y = yj.true)) +
  geom_point() +
  geom_point(data = beta.res[beta.res$yj.rel.diff > 100, ], colour = "red") +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw()
ggplot(beta.res, aes(x = yj.true, y = yj.diff)) +
  geom_point() +
  theme_bw()
ggplot(beta.res, aes(x = Mj, y = yj.rel.diff)) +
  geom_point() +
  theme_bw()

beta.res$yj.est <- NULL
beta.res$yj.true <- NULL
beta.res.long <- gather(beta.res, key = param, value = value, -c(simno, ind, Mj, xbar))
beta.res.long$type[beta.res.long$param %in% c("beta0", "beta1")] <- "true"
beta.res.long$type[beta.res.long$param %in% c("beta0.est", "beta1.est")] <- "est"
beta.res.long$param <- substr(beta.res.long$param, start = 1, stop = 5)
beta.res.long <- spread(beta.res.long, key = type, value = value)
beta.res.long$diff <- beta.res.long$true - beta.res.long$est

ggplot(beta.res.long, aes(x = Mj, y = diff)) +
  geom_point() +
  facet_wrap(~param) +
  theme_bw()

ggplot(beta.res.long, aes(x = est, y = true)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~param) +
  theme_bw()

all.res %>%
  group_by(params, whichmod) %>%
  summarise(bias = median(truth - posterior.mean),
            rel.bias = mean((truth - posterior.mean)/truth),
            sd = sqrt(mean((posterior.mean - mean(posterior.mean))^2)),
            rmse = sqrt(mean((truth - posterior.mean)^2)),
            cov50 = mean(p25 <= truth & truth <= p75),
            len50 = mean(p75 - p25),
            cov95 = mean(p025 <= truth & truth <= p975),
            len95 = mean(p975 - p025),
            med.n.eff = median(n.eff),
            med.Rhat = median(Rhat)) -> res.summ

print(cbind(res.summ$whichmod, res.summ$params, round(res.summ[, -c(1, 2)], digits = 3)))
