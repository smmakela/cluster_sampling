ff <- readRDS("/Users/susanna/projects/cluster_sampling/output/compiled_simulation_results_2016_11_28_2.rds")
dd <- ff[[1]]
ss <- ff[[2]]
library(dplyr)
library(ggplot2)
ss$num.units <- factor(ss$num.units,
                       levels = c("5pct", "10pct", "25pct", "50pct", "10", "50", "100"))
dd$num.units <- factor(dd$num.units,
                       levels = c("5pct", "10pct", "25pct", "50pct", "10", "50", "100"))
ss$model.name <- factor(ss$model.name,
                        levels = c("knowsizes", "cluster_inds_only", "hajek", "greg"))
dd <- dplyr::filter(dd, model.name %in% c("knowsizes", "cluster_inds_only"))
dd$model.name <- factor(ss$model.name,
                        levels = c("knowsizes", "cluster_inds_only"))
ss %>%
  group_by(model.name, use.sizes, num.clusters, num.units) %>%
  summarise(ybar.hat.se = sd(ybar.hat),
            ybar.true = mean(ybar.true),
            bias = mean(ybar.true - ybar.hat),
            rel.bias = mean((ybar.true - ybar.hat)/ybar.true),
            rmse = sqrt(mean((ybar.true - ybar.hat)^2)),
            cov50 = mean(ybar.hat.lci50 <= ybar.true & ybar.true <= ybar.hat.uci50),
            cov95 = mean(ybar.hat.lci95 <= ybar.true & ybar.true <= ybar.hat.uci95),
            len50 = mean(ybar.hat.uci50 - ybar.hat.lci50),
            len95 = mean(ybar.hat.uci95 - ybar.hat.lci95),
            ybar.hat = mean(ybar.hat),
            num.sims = n()) -> ss2
ss2 <- dplyr::filter(ss2, num.sims > 70)

names(dd) <- c("param.name", "est", "se_est", "sd", "p025", "p25", "p50",
               "p75", "p975", "n_eff", "Rhat", "model.name", "use.sizes",
               "outcome.type", "Mj", "ybar_true", "truth", "num.clusters",
               "num.units", "simno")
param.list <- c("gamma0", "gamma1", "alpha0", "alpha1", "sigma_beta0",
                "sigma_beta1", "sigma_y", "ybar_hat")
ytrue0 <- unique(dd$ybar_true[dd$use.sizes == 0])
ytrue0 <- ytrue0[!is.na(ytrue0)]
ytrue1 <- unique(dd$ybar_true[dd$use.sizes == 1])
ytrue1 <- ytrue1[!is.na(ytrue1)]
dd$ybar_true[dd$use.sizes == 0 & dd$param.name == "ybar_hat"] <- ytrue0
dd$ybar_true[dd$use.sizes == 1 & dd$param.name == "ybar_hat"] <- ytrue1
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

pdf(file = "/Users/susanna/projects/cluster_sampling/output/param_density_plots_separate.pdf",
    width = 10, height = 6, onefile = TRUE)
for (u in 0:1) {
  for (p in param.list) {
    dd.sub <- dplyr::filter(dd, param.name == p & use.sizes == u)
    pp <- ggplot(dd.sub, aes(x = est)) +
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
  
  
dd %>%
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
            num.sims = n()) -> dd2

ggplot(ss2, aes(x = factor(num.clusters), y = cov50)) +
  geom_point(aes(colour = model.name)) +
  geom_line(aes(colour = model.name, group = model.name)) +
  geom_hline(aes(yintercept = 0.5)) +
  facet_grid(use.sizes ~ num.units) +
  theme_bw()

ggplot(ss2, aes(x = factor(num.clusters), y = ybar.hat)) +
  geom_point(aes(colour = model.name),
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = ybar.hat - ybar.hat.se,
                    ymax = ybar.hat + ybar.hat.se,
                    colour = model.name),
                position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = ybar.true)) +
  facet_grid(use.sizes ~ num.units) +
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

