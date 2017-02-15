library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(shinystan)
source("src/rspps.r")
negbin_std <- stan_model(file = "src/cluster_size_model_testing/negbin.stan")
J <- 1000 # number of clusters in pop
K <- 50 # number of clusters to sample
k <- 200 # param for negative binomial
p <- 0.2 # param for negative binomial
mu <- k * (1-p) / p
phi <- k
#alpha <- k
#beta <- p / (1 - p)
# mean(Nj.pop) = k*(1-p)/p = 800
# var(Nj.pop) = k*(1-p)/p^2 = 4000

Nj.pop <- rnbinom(n = J, size = k, prob = p)

# Generate population data frame
pop.data <- data.frame(ids = c(1:J), Nj.pop)

n.sims <- 100
res.all <- data.frame()
ppd.all <- data.frame()
for (sim in 1:n.sims) {
  # Do PPS sampling
  sample.inds <- rspps(Nj.pop, c(1:J), K)
  Nj.sample <- Nj.pop[sample.inds]
  sample.data <- data.frame(sample.ids = c(1:K),
                            orig.ids = sample.inds, Nj.sample)
  
  # stan data
  stan.data <- list(K = K, Nj_sample_minus1 = Nj.sample - 1)
  
  # run stan
  stan.fit <- sampling(negbin_std, data = stan.data,
                       iter = 2000, chains = 4, control = list(adapt_delta = 0.99))
  
  # get output
  fit.summary <- data.frame(summary(stan.fit)$summary)
  fit.summary$param.name <- rownames(fit.summary)
  ff <- fit.summary[c("mu", "phi", "mu_star", "phi_star"),
                    c("param.name", "mean", "X2.5.", "X25.", "X75.", "X97.5.",
                      "n_eff", "Rhat")]
  ff$truth <- c(mu, phi)
  ff$simno <- sim
  res.all <- rbind(res.all, ff)
  
  # do posterior predictive draws
  Nj.pop.new <- rstan::extract(stan.fit, "Nj_new")
  Nj.pop.new <- unlist(Nj.pop.new, use.names = FALSE)
  Nj.sample.new <- rstan::extract(stan.fit, "Nj_new_star")
  Nj.sample.new <- unlist(Nj.sample.new, use.names = FALSE)
  df <- data.frame(Nj = c(Nj.pop, Nj.sample, Nj.pop.new, Nj.sample.new),
                   type = c(rep("pop", length(Nj.pop)),
                            rep("sample", length(Nj.sample)),
                            rep("PPD_pop", length(Nj.pop.new)),
                            rep("PPD_sampled", length(Nj.sample.new))))
  df$simno <- sim
  ppd.all <- rbind(ppd.all, df)
}
names(res.all) <- c("param.name", "postmean", "p025", "p25", "p75", "p975",
                    "n_eff", "Rhat", "truth", "simno")


ggplot(res.all[res.all$Rhat < 10 & res.all$n_eff > 20,], aes(x = postmean)) +
  geom_line(stat = "density") +
  geom_vline(aes(xintercept = truth)) +
  facet_wrap(~ param.name, scales = "free") +
  theme_bw()

res.all %>%
  group_by(param.name) %>%
  dplyr::filter(Rhat < 10 & n_eff > 10) %>%
  summarise(nobs = n(),
            cov50 = mean(p25 <= truth & truth <= p75),
            covg95 = mean(p025 <= truth & truth <= p975),
            bias = mean(truth - postmean),
            rmse = sqrt(mean((truth - postmean)^2)),
            rel.bias = mean((truth - postmean)/truth),
            rel.rmse = sqrt(mean(((truth - postmean)/truth)^2))) -> ff2


ggplot(ppd.all, aes(x = Nj)) +
  geom_line(stat = "density", aes(group = simno), colour = "grey50") +
  geom_vline(xintercept = mu, colour = "grey80") +
  #geom_line(data = )
  facet_wrap(~ type, scales = "free") +
  theme_bw()

ppd.all %>%
  dplyr::filter(type %in% c("PPD_pop", "PPD_sampled")) -> ppd2 #%>%
  #tidyr::spread(key = type, value = Nj) -> ppd2

ppd2 <- rbind(ppd2, data.frame(Nj = c(Nj.pop, Nj.sample),
                               simno = rep(0, times = J+K),
                               type = rep(c("pop", "sample"), times = c(J, K))))

ggplot(ppd2, aes(x = Nj)) +
  geom_line(stat = "density", aes(group = simno, colour = type)) +
  scale_colour_manual(breaks = c(""))
  theme_bw()

ppd2 <- ppd.all[ppd.all$type %in% c("PPD_pop", "PPD_sampled"),]

# coverage is for when you draw parameters from their priors, generate data,
# fit model, and then check whether truth is in interval (cook paper) -- dist of
# posterior *means* doesn't have the same coverage guarantees

