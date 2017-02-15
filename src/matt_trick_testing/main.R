# Purpose: see how reparameterization (the "Matt trick") is affected by the
# magnitude of group-level variance

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(tidyr)
library(dplyr)

# What should the data look like? should it have correlation between the coeffs
# or no? I think no b/c we want to compare the centered vs noncentered
# parameterization within a given model, not across models...

# Make a list of data frames for the 4 data scenarios


J.list <- c(5, 15, 30) # number of groups
tau.list <- c(0.001, 0.01, 0.1, 1, 10, 25) # group-level variance
mu <- 0.0
beta1 <- 0.5
sigma.y <- 0.2
n.obs.per.grp <- 30

#plot.pars <- c("mu", "tau", "beta1", "sigma_y")
plot.pars <- c("mu", "tau")

# varying intercept, no individual-level predictor
var_int_only_cp <- stan_model("src/matt_trick_testing/var_int_cp.stan")
var_int_only_ncp <- stan_model("src/matt_trick_testing/var_int_ncp.stan")
# varying intercept plus individual-level predictor
var_int_fix_slope_cp <- stan_model("src/matt_trick_testing/var_int_fixed_slope_cp.stan")
var_int_fix_slope_ncp <- stan_model("src/matt_trick_testing/var_int_fixed_slope_ncp.stan")
# varying intercept *and* slope (independent), plus individual-level predictor
var_int_slope_no_corr_cp <- stan_model("src/matt_trick_testing/var_int_slope_no_corr_cp.stan")
var_int_slope_no_corr_ncp <- stan_model("src/matt_trick_testing/var_int_slope_no_corr_ncp.stan")
# varying intercept *and* slope (correlated), plus individual-level predictor
var_int_slope_with_corr_cp <- stan_model("src/matt_trick_testing/var_int_slope_with_corr_cp.stan")
var_int_slope_with_corr_cp <- stan_model("src/matt_trick_testing/var_int_slope_with_corr_ncp.stan")

allres <- data.frame()
plot.posteriors <- data.frame()

#res.list <- vector(mode = "list", length = length(J.list) * length(tau.list))
df.all <- data.frame()
n.sims <- 1
#pdf(file = "/Users/susanna/projects/cluster_sampling/matt_test_pairs_plots.pdf",
#    width = 10, height = 8, onefile = TRUE)
for (s in 1:n.sims) {
  cat("========================== sim", s, "================================\n")
  for (k in 1:length(J.list)) {
    J <- J.list[k]
    n <- J*n.obs.per.grp # number of data points
    x <- rnorm(n, mean = 0, sd = 1)
    for (m in 1:length(tau.list)) {
      tau <- tau.list[m]
      beta0 <- rnorm(J, mean = mu, sd = tau)
      truepars <- data.frame(params = c("mu", "tau"),
                             true.val = c(mu, tau))
      truepars$params <- as.character(truepars$params)
      beta0.rep <- rep(beta0, each = n.obs.per.grp)
      #y <- beta0.rep + beta1 * x + rnorm(n, mean = 0, sd = sigma.y)
      y <- beta0.rep + rnorm(n, mean = 0, sd = sigma.y)
      df <- data.frame(x, y, grp.id = rep(c(1:J), each = n.obs.per.grp))
      df$J <- J
      df$tau <- tau
      df.all <- rbind(df.all, df)
#   ggplot(df.all, aes(x = y)) +
#     geom_line(stat = "density", aes(colour = factor(grp.id), group = factor(grp.id))) +
#     facet_grid(J ~ tau) +
#     theme_bw()
  #stan_data <- list(n = n, J = J, x = df$x, y = df$y, grp_id = df$grp.id)
      stan_data <- list(n = n, J = J, y = df$y, grp_id = df$grp.id, sigma_y = sigma.y)
      no_res <- sampling(var_int_nomatt_model, data = stan_data,
                         control = list(adapt_delta = 0.99))
      with_res <- sampling(var_int_withmatt_model, data = stan_data,
                           control = list(adapt_delta = 0.99))
      if (s == 1) {
        pairs(no_res, pars = c("mu", "tau"))
        title(paste0("J = ", J, ", tau = ", tau, ", NO MATT"), line = 2.5)
        pairs(with_res, pars = c("mu", "tau"))
        title(paste0("J = ", J, ", tau = ", tau, ", WITH MATT"), line = 2.5)
      }
      summ.no <- data.frame(summary(no_res)$summary)
      summ.with <- data.frame(summary(with_res)$summary)
      params <- c(rownames(summ.no), rownames(summ.with))
      post.means <- c(summ.no$mean, summ.with$mean)
      post.sds <- c(summ.no$sd, summ.with$sd)
      lci.95 <- c(summ.no[["X2.5."]], summ.with[["X2.5."]])
      uci.95 <- c(summ.no[["X97.5."]], summ.with[["X97.5."]])
      lci.50 <- c(summ.no[["X25."]], summ.with[["X25."]])
      uci.50 <- c(summ.no[["X75."]], summ.with[["X75."]])
      Rhats <- c(summ.no$Rhat, summ.with$Rhat)
      n_effs <- c(summ.no$n_eff, summ.with$n_eff)
      mod <- rep(c("no", "with"), times = c(nrow(summ.no), nrow(summ.with)))
      resdf <- data.frame(params, post.means, post.sds, lci.95, uci.95,
                          lci.50, uci.50, Rhats, n_effs, mod)
      resdf$J <- J
      resdf$tau <- tau
      resdf$params <- as.character(resdf$params)
      resdf$simno <- s
      resdf <- left_join(resdf, truepars, by = "params")
      allres <- rbind(allres, resdf)
      if (s == 1) {
        samps.no <- data.frame(rstan::extract(no_res, pars = plot.pars))
        samps.no <- tidyr::gather(samps.no, key = "params", value = "samp")
        samps.with <- data.frame(rstan::extract(with_res, pars = plot.pars))
        samps.with <- tidyr::gather(samps.with, key = "params", value = "samp")
        samps.all <- rbind(samps.no, samps.with)
        samps.all$J <- J
        samps.all$tau <- tau
        samps.all$simno <- s
        samps.all$mod <- rep(c("no", "with"),
                             times = c(nrow(samps.no), nrow(samps.with)))
        samps.all$params <- as.character(samps.all$params)
        samps.all <- left_join(samps.all, truepars, "params")
        plot.posteriors <- rbind(plot.posteriors, samps.all)
      }
    }
  }
}
#dev.off()

  ggplot(df.all, aes(x = y)) +
    geom_line(stat = "density", aes(colour = factor(grp.id), group = factor(grp.id))) +
    facet_grid(J ~ tau) +
    theme_bw()

allres_orig <- allres
allres <- dplyr::filter(allres, !grepl("eta\\[", params))

# calculate coverage of 50%, 95% intervals
allres %>%
  dplyr::group_by(params, J, tau, mod) %>%
  dplyr::summarise(covg_50 = mean(lci.50 <= true.val & true.val <= uci.50),
                   covg_95 = mean(lci.95 <= true.val & true.val <= uci.95),
                   bias = mean(true.val - post.means),
                   rmse = sqrt(mean((true.val - post.means)^2))) -> bb
bbsub <- bb[bb$params %in% c("mu", "tau"), ]
bbsub <- gather(bbsub, key = "level", value = "coverage", c(covg_50, covg_95))
bbsub$level <- as.numeric(substr(bbsub$level, 6, 8))
ggplot(bbsub, aes(x = factor(tau), y = coverage)) +
  geom_point(aes(colour = mod, shape = factor(level)), size = 2) +
  geom_hline(yintercept = c(0.5, 0.95), colour = "grey80") +
  scale_shape_manual("coverage", values = c(1, 4)) +
  facet_grid(J ~ params) +
  theme_bw()
ggsave("coverage.pdf", width = 6.5, height = 5)

# ggplot(bbsub, aes(x = factor(tau), y = bias)) +
#   geom_point(aes(colour = mod)) +
#   geom_hline(yintercept = 0, colour = "grey80") +
#   facet_grid(params ~ J, scales = "free") +
#   theme_bw()
# ggsave("coverage.pdf", width = 10, height = 8)

# boxplot of effective sample sizes
ggplot(allres[allres$params %in% plot.pars,],
       aes(x = factor(tau), y = n_effs)) +
  geom_boxplot(aes(colour = mod),
               position = position_dodge(width = 0.5)) +
  facet_grid(params ~ J, scales = "free") +
  scale_colour_discrete(name = "model") +
  xlab("tau") +
  ylab("n_eff") +
  theme_bw()
ggsave("boxplot_neff.pdf", width = 9, height = 5)

# boxplot of Rhats
ggplot(allres[allres$params %in% plot.pars,],
       aes(x = factor(tau), y = Rhats)) +
  geom_boxplot(aes(colour = mod),
               position = position_dodge(width = 0.5)) +
  facet_grid(params ~ J, scales = "free") +
  scale_colour_discrete(name = "model") +
  xlab("tau") +
  ylab("Rhat") +
  theme_bw()
ggsave("boxplot_Rhat.pdf", width = 9, height = 5)

# boxplot of bias
allres$diff <- allres$true.val - allres$post.means
ggplot(allres[allres$params %in% plot.pars,],
       aes(x = factor(tau), y = diff)) +
  geom_boxplot(aes(colour = mod),
               position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  facet_grid(params ~ J, scales = "free") +
  scale_colour_discrete(name = "model") +
  xlab("tau") +
  ylab("truth - posterior mean") +
  theme_bw()
ggsave("boxplot_diff.pdf", width = 8, height = 5)

plot.posteriors$tau_facet <- paste0("tau = ", plot.posteriors$tau)
plot.posteriors$J_facet <- paste0("J = ", plot.posteriors$J)
plot.posteriors$J_facet <- factor(plot.posteriors$J_facet,
                                  levels = c("J = 5", "J = 15", "J = 30"))
ggplot(plot.posteriors[plot.posteriors$params == "tau",], aes(x = samp)) +
  geom_line(stat = "density", aes(colour = mod)) +
  geom_vline(aes(xintercept = true.val), colour = "gray80") +
  facet_wrap(J_facet ~ tau_facet, scales = "free", ncol = 6) +
  ylab("posterior for tau") +
  theme_bw()
ggsave("posterior_densities_tau.pdf", width = 12, height = 6)
ggplot(plot.posteriors[plot.posteriors$params == "mu",], aes(x = samp)) +
  geom_line(stat = "density", aes(colour = mod)) +
  geom_vline(aes(xintercept = true.val), colour = "gray80") +
  facet_wrap(J_facet ~ tau_facet, scales = "free", ncol = 6) +
  ylab("posterior for mu") +
  theme_bw()
ggsave("posterior_densities_mu.pdf", width = 12, height = 6)
