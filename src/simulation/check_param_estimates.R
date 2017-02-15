rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)

allres <- readRDS("/Users/susanna/projects/cluster_sampling/output/compiled_simulation_results_2016_12_12.rds")
param.ests <- allres[[1]]
param.ests %>%
  group_by(param.name, model.name, use.sizes, num.clusters, num.units) %>%
  mutate(num.obs = n()) -> param.ests
param.ests$num.units <- factor(param.ests$num.units,
                               levels = c("5pct", "10pct", "25pct", "50pct", "100pct", "10", "50", "100"))
param.ests$model.name <- as.character(param.ests$model.name)
param.ests$model.name[param.ests$model.name == "knowsizes"] <- "knowsizes.stan"
param.ests$model.name[param.ests$model.name == "with_cluster_sizes"] <- "knowsizes.lmer"
param.ests$model.name[param.ests$model.name == "cluster_inds_only"] <- "cluster_inds_only.stan"
param.ests$model.name[param.ests$model.name == "cluster_indicators_only"] <- "cluster_inds_only.lmer"

param.ests.orig <- param.ests
tmp <- dplyr::filter(param.ests, grepl("*stan", model.name))
tmp$model.name[tmp$model.name == "knowsizes.stan"] <- "knowsizes"
tmp$model.name[tmp$model.name == "cluster_inds_only.stan"] <- "cluster_inds_only"
tmp$method <- "stan"
tmp2 <- dplyr::filter(param.ests, grepl("*lmer", model.name))
tmp2$model.name[tmp2$model.name == "knowsizes.lmer"] <- "knowsizes"
tmp2$model.name[tmp2$model.name == "cluster_inds_only.lmer"] <- "cluster_inds_only"
tmp2$method <- "lmer"
param.ests <- rbind(tmp, tmp2)
param.ests$method <- factor(param.ests$method,
                            levels = c("stan", "lmer"))
#names(param.ests) <- c("param.name", "est", "se_est", "sd", "p025", "p25", "p50",
#                       "p75", "p975", "n_eff", "Rhat", "model.name", "use.sizes",
#                       "outcome.type", "Mj", "ybar_true", "truth", "num.clusters",
#                       "num.units", "simno", "num.obs", "method")
names(param.ests) <- c("param.name", "est", "se_est", "sd", "p025", "p25", "p50",
               "p75", "p975", "n_eff", "Rhat", "model.name", "use.sizes",
               "outcome.type", "Mj", "ybar_true", "truth", "num.clusters",
               "num.units", "simno", "num.obs")
param.list <- c("alpha0", "gamma0", "alpha1", "gamma1", "sigma_beta0",
                "sigma_beta1", "sigma_y", "ybar_new")
ytrue0 <- unique(param.ests$ybar_true[param.ests$use.sizes == 0])
ytrue0 <- ytrue0[!is.na(ytrue0)]
ytrue1 <- unique(param.ests$ybar_true[param.ests$use.sizes == 1])
ytrue1 <- ytrue1[!is.na(ytrue1)]
param.ests$truth[param.ests$use.sizes == 0 & param.ests$param.name == "ybar_new"] <- ytrue0
param.ests$truth[param.ests$use.sizes == 1 & param.ests$param.name == "ybar_new"] <- ytrue1


pdf(file = "/Users/susanna/projects/cluster_sampling/output/param_dot_plots_onesim.pdf",
    width = 10, height = 6, onefile = TRUE)
  for (p in param.list) {
    data.sub <- dplyr::filter(param.ests, param.name == p)
    pp <- ggplot(data.sub, aes(x = factor(num.clusters), y = est,
                           group = model.name, colour = model.name)) +
      geom_point(position = position_dodge(width = 0.9)) +
      geom_errorbar(aes(ymin = p25, ymax = p75),
                    position = position_dodge(width = 0.9)) +
      geom_hline(aes(yintercept = truth)) +
      facet_grid(use.sizes ~ num.units, scales = "free") +
      ggtitle(p) +
      theme_bw()
    print(pp)
}
dev.off()

pdf(file = "/Users/susanna/projects/cluster_sampling/output/param_density_plots_separate_12_12.pdf",
    width = 10, height = 6, onefile = TRUE)
for (u in 0:1) {
  for (p in param.list) {
    data.sub <- dplyr::filter(param.ests, param.name == p & use.sizes == u)
    pp <- ggplot(data.sub, aes(x = est)) +
      geom_line(stat = "density", aes(colour = model.name)) +
      facet_grid(num.units ~ num.clusters, scales = "free") +
      ggtitle(paste0(p, ", use.sizes = ", u)) +
      theme_bw()
    if (p == "ybar_hat") {
      pp <- pp + geom_vline(aes(xintercept = ybar_true))
    } else {
      pp <- pp + geom_vline(aes(xintercept = truth))
    }
    print(pp)
  }
}
dev.off()

pdf(file = "/Users/susanna/projects/cluster_sampling/output/param_density_plots_separate.pdf",
    width = 10, height = 6, onefile = TRUE)
for (u in 0:1) {
  for (p in param.list) {
    data.sub <- dplyr::filter(param.ests, param.name == p & use.sizes == u)
    pp <- ggplot(data.sub, aes(x = est)) +
      geom_line(stat = "density", aes(colour = model.name, linetype = method)) +
      facet_grid(num.units ~ num.clusters, scales = "free") +
      ggtitle(paste0(p, ", use.sizes = ", u)) +
      theme_bw()
    if (p == "ybar_hat") {
      pp <- pp + geom_vline(aes(xintercept = ybar_true))
    } else {
      pp <- pp + geom_vline(aes(xintercept = truth))
    }
    print(pp)
  }
}
dev.off()






pdf(file = "/Users/susanna/projects/cluster_sampling/output/param_density_plots.pdf",
    width = 10, height = 6, onefile = TRUE)
for (p in param.list) {
  dd.sub <- dplyr::filter(dd, param.name == p)
  pp <- ggplot(dd.sub, aes(x = est)) +
    geom_line(stat = "density", aes(linetype = factor(use.sizes),
                                    colour = model.name)) +
    facet_grid(num.units ~ num.clusters, scales = "free") +
    ggtitle(p) +
    theme_bw()
  if (p == "ybar_hat") {
    pp <- pp + geom_vline(aes(xintercept = ybar_true, linetype = factor(use.sizes)))
  } else {
    pp <- pp + geom_vline(aes(xintercept = truth, linetype = factor(use.sizes)))
  }
  print(pp)
}
dev.off()


ybar.ests <- allres[[2]]
ybar.ests %>%
  group_by(model.name, use.sizes, num.clusters, num.units) %>%
  summarise(ybar.hat.sd = sd(ybar.hat),
            ybar.hat.est = mean(ybar.hat),
            bias = mean(ybar.true - ybar.hat),
            rel.bias = mean((ybar.true - ybar.hat)/ybar.true),
            rmse = sqrt(mean((ybar.true - ybar.hat)^2)),
            cov50 = mean(ybar.hat.lci50 <= ybar.true & ybar.true <= ybar.hat.uci50),
            cov95 = mean(ybar.hat.lci95 <= ybar.true & ybar.true <= ybar.hat.uci95),
            len50 = mean(ybar.hat.uci50 - ybar.hat.lci50),
            len95 = mean(ybar.hat.uci95 - ybar.hat.lci95),
            ybar.true = mean(ybar.true),
            num.sims = n()) -> ybar.ests.summ
ybar.ests %>%
  group_by(model.name, use.sizes, num.clusters, num.units, param.name) %>%
  summarise(param.sd = sd(est),
            param.est = mean(est),
            bias = mean(truth - est),
            rel.bias = mean((truth - est)/truth),
            rmse = sqrt(mean((truth - est)^2)),
            cov50 = mean(p25 <= truth & truth <= p75),
            cov95 = mean(p025 <= truth & truth <= p975),
            len50 = mean(ybar.hat.uci50 - ybar.hat.lci50),
            len95 = mean(ybar.hat.uci95 - ybar.hat.lci95),
            ybar.hat = mean(ybar.hat),
            num.sims = n()) -> ybar.ests.summ

ggplot(ss2, aes(x = factor(num.clusters), y = cov50)) +
  geom_point(aes(colour = model.name)) +
  geom_line(aes(colour = model.name, group = model.name)) +
  geom_hline(aes(yintercept = 0.5)) +
  facet_grid(use.sizes ~ num.units) +
  theme_bw()

ggplot(ybar.ests.summ, aes(x = factor(num.clusters), y = ybar.hat.est)) +
  geom_point(aes(colour = model.name),
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = ybar.hat.est - ybar.hat.sd,
                    ymax = ybar.hat.est + ybar.hat.sd,
                    colour = model.name),
                position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = ybar.true)) +
  facet_grid(use.sizes ~ num.units, scales = "free") +
  theme_bw()

ggplot(ss2, aes(x = factor(num.clusters), y = bias)) +
  geom_point(aes(colour = model.name)) +
  geom_line(aes(colour = model.name, group = model.name)) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(use.sizes ~ num.units) +
  theme_bw()

ggplot(ss2, aes(x = factor(num.clusters), y = rel.bias)) +
  geom_point(aes(colour = model.name)) +
  geom_line(aes(colour = model.name, group = model.name)) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(use.sizes ~ num.units) +
  theme_bw()

ggplot(ss2, aes(x = factor(num.clusters), y = rmse)) +
  geom_point(aes(colour = model.name)) +
  geom_line(aes(colour = model.name, group = model.name)) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(use.sizes ~ num.units) +
  theme_bw()

ggplot(ss2, aes(x = factor(num.clusters), y = ybar.hat)) +
  geom_point(aes(colour = model.name),
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = ybar.hat - len95/2,
                    ymax = ybar.hat + len95/2,
                    colour = model.name),
                width = 0,
                position = position_dodge(width = 0.9),
                size = 0.75) +
  geom_errorbar(aes(ymin = ybar.hat - len50/2,
                    ymax = ybar.hat + len50/2,
                    colour = model.name),
                width = 0,
                position = position_dodge(width = 0.9),
                size = 1.3) +
  geom_hline(aes(yintercept = ybar.true)) +
  facet_grid(use.sizes ~ num.units) +
  theme_bw()


ss2.long <- gather(ss2, key = metric, value = value, cov50, cov95, len50, len95)
ss2.long <- extract(ss2.long, col = "metric", into = c("metric", "level"),
                    regex = "([[:alpha:]]*)([[:digit:]]*)")
ss2.long <- spread(ss2.long, key = metric, value = value)
ss2.long$model.level <- paste0(ss2.long$model.name, ss2.long$level)
ggplot(ss2.long, aes(x = factor(num.clusters), y = cov)) +
  geom_point(aes(colour = model.name)) +
  geom_line(aes(colour = model.name, group = model.level, linetype = level)) +
  geom_hline(aes(yintercept = 0.5)) +
  geom_hline(aes(yintercept = 0.95)) +
  facet_grid(use.sizes ~ num.units) +
  theme_bw()

ggplot(ss2, aes(x = factor(num.clusters), y = len50)) +
  geom_point(aes(colour = model.name)) +
  geom_line(aes(colour = model.name, group = model.name)) +
  facet_grid(use.sizes ~ num.units) +
  theme_bw()

ggplot(ss, aes(x = factor(num.clusters), y = ybar.hat)) +
  #geom_point(aes(colour = model.name), alpha = 0.5) +
  geom_boxplot(aes(colour = model.name), position = "dodge") +
  geom_hline(aes(yintercept = ybar.true)) +
  facet_grid(use.sizes ~ num.units) +
  theme_bw()

