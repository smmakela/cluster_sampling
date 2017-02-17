# Purpose: check that the parameter estimates and ybar estimateas from the basic
# models (where we know all the cluster sizes and where we only use cluster
# indicators) are close to the truth

# ###################################
# Prep
# ###################################
rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
allres <- readRDS("/Users/susanna/projects/cluster_sampling/output/compiled_simulation_results_2016_12_13.rds")

# ###################################
# Parameter estimate plotss
# ###################################
param.ests <- allres[[1]]
param.ests %>%
  group_by(param.name, model.name, use.sizes, num.clusters, num.units) %>%
  mutate(num.obs = n()) -> param.ests
param.ests$num.units <- factor(param.ests$num.units,
                               levels = c("5pct", "10pct", "25pct", "50pct", "100pct", "10", "50", "100"))
param.ests <- dplyr::filter(param.ests,
                            model.name %in% c("cluster_inds_only", "knowsizes"))
param.ests$model.name <- factor(param.ests$model.name,
                                levels = c("knowsizes", "cluster_inds_only"))
param.ests$truth <- as.numeric(param.ests$truth)
names(param.ests) <- c("param.name", "est", "se_est", "sd", "p025", "p25", "p50",
                       "p75", "p975", "n_eff", "Rhat", "model.name", "use.sizes",
                       "outcome.type", "M_tot", "truth", "num.clusters",
                       "num.units", "simno", "num.obs")
param.list <- c("alpha0", "gamma0", "alpha1", "gamma1", "sigma_beta0",
                "sigma_beta1", "sigma_y", "ybar_new")

# Density plots of posterior means across many simulations
pdf(file = "/Users/susanna/projects/cluster_sampling/output/param_density_plots_12_12.pdf",
    width = 10, height = 6, onefile = TRUE)
for (u in 0:1) {
  for (p in param.list) {
    data.sub <- dplyr::filter(param.ests, param.name == p & use.sizes == u)
    pp <- ggplot(data.sub, aes(x = est)) +
      geom_line(stat = "density", aes(colour = model.name)) +
      facet_grid(num.units ~ num.clusters, scales = "free") +
      geom_vline(aes(xintercept = truth)) +
      ggtitle(paste0(p, ", use.sizes = ", u)) +
      theme_bw()
    print(pp)
  }
}
dev.off()

# Density plots of posterior means centered at truth across many simulations
pdf(file = "/Users/susanna/projects/cluster_sampling/output/param_centered_density_plots_12_12.pdf",
    width = 10, height = 6, onefile = TRUE)
for (p in param.list) {
  data.sub <- dplyr::filter(param.ests, param.name == p)
  pp <- ggplot(data.sub, aes(x = est - truth)) +
    geom_line(stat = "density",
              aes(colour = model.name, linetype = factor(use.sizes))) +
    facet_grid(num.units ~ num.clusters, scales = "free") +
    geom_vline(xintercept = 0) +
    ggtitle(p) +
    theme_bw()
  print(pp)
}
dev.off()

# Calculate summary statistics
param.ests %>%
  group_by(model.name, use.sizes, num.clusters, num.units, param.name) %>%
  summarise(param.sd = sd(est),
            param.est = mean(est),
            bias = mean(truth - est),
            rel.bias = mean((truth - est)/truth),
            rmse = sqrt(mean((truth - est)^2)),
            cov50 = mean(p25 <= truth & truth <= p75),
            cov95 = mean(p025 <= truth & truth <= p975),
            len50 = mean(p75 - p25),
            len95 = mean(p975 - p025),
            num.sims = n()) -> param.ests.summ

# Plot coverages -- need to make long first
param.ests.summ %>%
  select(-c(param.sd, param.est, bias, rel.bias, rmse, num.sims)) -> ff
ff %>%
  gather(key = varname, value = value, c(contains("50"), contains("95"))) %>%
  extract(varname, into = c("stat", "level"),
           regex = "([[:alpha:]]*)([[:digit:]]*)") %>% 
  mutate(level = as.numeric(level)) %>%
  spread(key = stat, value = value) -> ff

pdf(file = "/Users/susanna/projects/cluster_sampling/output/param_coverage_plots_12_13.pdf",
    width = 10, height = 6, onefile = TRUE)
for (p in param.list) {
  ff2 <- dplyr::filter(ff, param.name == p)
  pp <- ggplot(ff2, aes(x = factor(num.clusters), y = cov)) +
    geom_point(aes(colour = model.name, shape = factor(level)),
               position = position_dodge(width = 0.9)) +
    geom_hline(yintercept = c(0.5, 0.95), colour = "gray") +
    scale_shape_manual(name = "nominal coverage %",
                       values = c(1, 4)) +
    scale_colour_discrete(name = "model") +
    facet_grid(use.sizes ~ num.units) +
    xlab("Number of clusters") +
    ylab("Coverage") +
    ggtitle(p) +
    theme_bw() +
    theme(legend.position = "bottom", legend.box = "horizontal")
  print(pp)
}
dev.off()

# CI lengths
pdf(file = "/Users/susanna/projects/cluster_sampling/output/param_CIlen_plots_12_12.pdf",
    width = 10, height = 6, onefile = TRUE)
for (p in param.list) {
  ff2 <- dplyr::filter(ff, param.name == p)
  pp <- ggplot(ff2, aes(x = factor(num.clusters), y = len)) +
    geom_point(aes(colour = model.name, shape = factor(level)),
               position = position_dodge(width = 0.9)) +
    scale_shape_manual(name = "nominal coverage %",
                       values = c(1, 4)) +
    scale_colour_discrete(name = "model") +
    facet_grid(use.sizes ~ num.units) +
    xlab("Number of clusters") +
    ylab("CI Length") +
    ggtitle(p) +
    theme_bw() +
    theme(legend.position = "bottom", legend.box = "horizontal")
  print(pp)
}
dev.off()

# Bias
pdf(file = "/Users/susanna/projects/cluster_sampling/output/param_bias_plots_12_13.pdf",
    width = 10, height = 6, onefile = TRUE)
for (p in param.list) {
  ff2 <- dplyr::filter(param.ests.summ, param.name == p)
  pp <- ggplot(ff2, aes(x = factor(num.clusters), y = bias)) +
    geom_point(aes(colour = model.name),
               position = position_dodge(width = 0.9)) +
    geom_hline(yintercept = 0, colour = "gray") +
    scale_colour_discrete(name = "model") +
    facet_grid(use.sizes ~ num.units) +
    xlab("Number of clusters") +
    ylab("Bias") +
    ggtitle(p) +
    theme_bw() +
    theme(legend.position = "bottom", legend.box = "horizontal")
  print(pp)
}
dev.off()

# Relative bias
pdf(file = "/Users/susanna/projects/cluster_sampling/output/param_rel_bias_plots_12_13.pdf",
    width = 10, height = 6, onefile = TRUE)
for (p in param.list) {
  ff2 <- dplyr::filter(param.ests.summ, param.name == p)
  pp <- ggplot(ff2, aes(x = factor(num.clusters), y = rel.bias)) +
    geom_point(aes(colour = model.name),
               position = position_dodge(width = 0.9)) +
    geom_hline(yintercept = 0, colour = "gray") +
    scale_colour_discrete(name = "model") +
    facet_grid(use.sizes ~ num.units) +
    xlab("Number of clusters") +
    ylab("Relative bias ((truth - est) / truth)") +
    ggtitle(p) +
    theme_bw() +
    theme(legend.position = "bottom", legend.box = "horizontal")
  print(pp)
}
dev.off()

# RMSE
pdf(file = "/Users/susanna/projects/cluster_sampling/output/param_rmse_plots_12_13.pdf",
    width = 10, height = 6, onefile = TRUE)
for (p in param.list) {
  ff2 <- dplyr::filter(param.ests.summ, param.name == p)
  pp <- ggplot(ff2, aes(x = factor(num.clusters), y = rmse)) +
    geom_point(aes(colour = model.name),
               position = position_dodge(width = 0.9)) +
    geom_hline(yintercept = 0, colour = "gray") +
    scale_colour_discrete(name = "model") +
    facet_grid(use.sizes ~ num.units) +
    xlab("Number of clusters") +
    ylab("RMSE") +
    ggtitle(p) +
    theme_bw() +
    theme(legend.position = "bottom", legend.box = "horizontal")
  print(pp)
}
dev.off()


# ###################################
# Ybar estimates
# ###################################
ybar.ests <- allres[[2]]
ybar.ests$num.units <- factor(ybar.ests$num.units,
                              levels = c("5pct", "10pct", "25pct", "50pct", "100pct", "10", "50", "100"))
ybar.ests <- dplyr::filter(ybar.ests,
                           model.name %in% c("cluster_inds_only", "knowsizes",
                                             "greg", "hajek"))
ybar.ests$model.name <- factor(ybar.ests$model.name,
                               levels = c("knowsizes", "cluster_inds_only",
                                          "greg", "hajek"))

ggplot(ybar.ests, aes(x = factor(num.clusters), y = ybar.hat)) +
  geom_boxplot(aes(colour = model.name), position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = ybar.true), colour = "gray") +
  facet_grid(use.sizes ~ num.units, scales = "free") +
  theme_bw()

for (u in c(0, 1)) {
  gg <- dplyr::filter(ybar.ests, use.sizes == u)
  pp <- ggplot(gg, aes(x = ybar.hat)) +
    geom_line(stat = "density", aes(colour = factor(num.clusters))) +
    geom_vline(aes(xintercept = ybar.true), colour = "gray") +
    scale_color_manual(values = c("#bdd7e7", "#6baed6", "#3182bd", "#08519c")) +
    facet_grid(model.name ~ num.units) +
    ggtitle(paste0("use.sizes = ", u)) +
    theme_bw()
  print(p)
}

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

# Bias
ggplot(ybar.ests.summ, aes(x = factor(num.clusters), y = bias)) +
  geom_point(aes(colour = model.name)) +
  geom_line(aes(colour = model.name, group = model.name)) +
  geom_hline(yintercept = 0, colour = "gray") +
  facet_grid(use.sizes ~ num.units, scales = "free") +
  theme_bw()
# ggplot(ybar.ests.summ, aes(x = factor(num.clusters), y = bias)) +
#   geom_point(aes(colour = model.name),
#              position = position_dodge(width = 0.9)) +
#   geom_hline(yintercept = 0, colour = "gray") +
#   facet_grid(num.units ~ use.sizes) +
#   theme_bw()


# Relative bias
ggplot(ybar.ests.summ, aes(x = factor(num.clusters), y = rel.bias)) +
  geom_point(aes(colour = model.name)) +
  geom_line(aes(colour = model.name, group = model.name)) +
  geom_hline(yintercept = 0, colour = "gray") +
  facet_grid(use.sizes ~ num.units, scales = "free") +
  theme_bw()
# ggplot(ybar.ests.summ, aes(x = factor(num.clusters), y = rel.bias)) +
#   geom_point(aes(colour = model.name),
#              position = position_dodge(width = 0.9)) +
#   geom_hline(yintercept = 0, colour = "gray") +
#   facet_grid(num.units ~ use.sizes) +
#   theme_bw()

# RMSE
ggplot(ybar.ests.summ, aes(x = factor(num.clusters), y = rmse)) +
  geom_point(aes(colour = model.name), size = 2) +
  geom_line(aes(colour = model.name, group = model.name)) +
  geom_hline(yintercept = 0, colour = "gray") +
  facet_grid(use.sizes ~ num.units, scales = "free") +
  theme_bw()
# ggplot(ybar.ests.summ, aes(x = factor(num.clusters), y = rmse)) +
#   geom_point(aes(colour = model.name, shape = factor(use.sizes)),
#              size = 2, position = position_dodge(width = 0.9)) +
#   geom_hline(yintercept = 0, colour = "gray") +
#   scale_shape_manual(values = c(1, 4)) +
#   facet_wrap(~num.units, nrow = 2) +
#   theme_bw()


ggplot(ybar.ests.summ, aes(x = rmse, y = model.name)) +
  geom_point(aes(colour = factor(num.clusters))) +
  facet_grid(num.units ~ use.sizes, scales = "free") +
  theme_bw()


# Coverage -- make long
ybar.ests.summ %>%
  select(-c(ybar.hat.sd, bias, rel.bias, rmse, num.sims)) -> ff
ff %>%
  gather(key = varname, value = value, c(contains("50"), contains("95"))) %>%
  extract(varname, into = c("stat", "level"),
          regex = "([[:alpha:]]*)([[:digit:]]*)") %>% 
  mutate(level = as.numeric(level)) %>%
  spread(key = stat, value = value) -> ff

ggplot(ff, aes(x = factor(num.clusters), y = cov)) +
  geom_point(aes(colour = model.name, shape = factor(level))) +
  geom_line(aes(colour = model.name, group = paste0(model.name, level),
                linetype = factor(level))) +
  geom_hline(yintercept = c(0.5, 0.95), colour = "gray") +
  scale_shape_manual(name = "nominal coverage %",
                     values = c(1, 4)) +
  scale_colour_discrete(name = "model") +
  facet_grid(use.sizes ~ num.units) +
  xlab("Number of clusters") +
  ylab("Coverage") +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "horizontal")

ggplot(ff, aes(x = factor(num.clusters), y = len/abs(ybar.true))) +
  geom_point(aes(colour = model.name, shape = factor(level))) +
  geom_line(aes(colour = model.name, group = paste0(model.name, level),
                linetype = factor(level))) +
  scale_shape_manual(name = "nominal coverage %",
                     values = c(1, 4)) +
  scale_colour_discrete(name = "model") +
  facet_grid(use.sizes ~ num.units, scales = "free") +
  xlab("Number of clusters") +
  ylab("Relative UI length") +
  theme_bw() +
  theme(legend.position = "bottom")




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

