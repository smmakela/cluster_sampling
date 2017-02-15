library(rstan)
library(tidyr)
library(dplyr)
simdat0 <- readRDS("/Users/susanna/projects/cluster_sampling/output/simdata_usesizes_0_continuous_nclusters_999_nunits_999_simno_999_nomatt.rds")
simres0 <- readRDS("/Users/susanna/projects/cluster_sampling/output/ccc_usesizes_0_continuous_knowsizes_ccc_nomatt.rds")
popdata0 <- readRDS("/Users/susanna/projects/cluster_sampling/output/popdata_usesizes_0_continuous_nomatt.rds")
load("/Users/susanna/projects/cluster_sampling/output/stanfit_usesizes_0_nclusters_999_nunits_999_simno_999_knowsizes_ccc_nomatt.RData")
fit0 <- fit
rm(fit)
simdat1 <- readRDS("/Users/susanna/projects/cluster_sampling/output/simdata_usesizes_1_continuous_nclusters_999_nunits_999_simno_999_nomatt.rds")
simres1 <- readRDS("/Users/susanna/projects/cluster_sampling/output/ccc_usesizes_1_continuous_knowsizes_ccc_nomatt.rds")
popdata1 <- readRDS("/Users/susanna/projects/cluster_sampling/output/popdata_usesizes_1_continuous_nomatt.rds")
load("/Users/susanna/projects/cluster_sampling/output/stanfit_usesizes_1_nclusters_999_nunits_999_simno_999_knowsizes_ccc_nomatt.RData")
fit1 <- fit
rm(fit)

param.true.values <- popdata0[["truepars"]]
param.true.values$ybar_new <- mean(popdata0[["pop.data"]]$y)
param.true.values <- gather(param.true.values, key = param.name, value = true.value)
param.true.values$param.name <- as.character(param.true.values$param.name)
param.post.summ <- simres0[["par.ests"]]
param.post.summ$param.name <- as.character(param.post.summ$param.name)
param.post.summ <- left_join(param.post.summ, param.true.values, by = "param.name")
tt <- param.post.summ[, c("param.name", "true.value", "mean", "2.5%", "97.5%", "25%", "75%")]
tt[, c("true.value", "mean", "2.5%", "97.5%", "25%", "75%")] <-
  round(tt[, c("true.value", "mean", "2.5%", "97.5%", "25%", "75%")], digits = 2)
print(tt)
traceplot(fit0, pars = c("alpha0", "gamma0", "alpha1", "gamma1", "sigma_beta0", "sigma_beta1", "sigma_y", "ybar_new"))
print(fit0, pars = c("alpha0", "gamma0", "alpha1", "gamma1", "sigma_beta0", "sigma_beta1", "sigma_y", "ybar_new"))

param.true.values <- popdata1[["truepars"]]
param.true.values$ybar_new <- mean(popdata1[["pop.data"]]$y)
param.true.values <- gather(param.true.values, key = param.name, value = true.value)
param.true.values$param.name <- as.character(param.true.values$param.name)
param.post.summ <- simres1[["par.ests"]]
param.post.summ$param.name <- as.character(param.post.summ$param.name)
param.post.summ <- left_join(param.post.summ, param.true.values, by = "param.name")
tt <- param.post.summ[, c("param.name", "true.value", "mean", "2.5%", "97.5%", "25%", "75%")]
tt[, c("true.value", "mean", "2.5%", "97.5%", "25%", "75%")] <-
  round(tt[, c("true.value", "mean", "2.5%", "97.5%", "25%", "75%")], digits = 2)
print(tt)
traceplot(fit1, pars = c("alpha0", "gamma0", "alpha1", "gamma1", "sigma_beta0", "sigma_y", "ybar_new"))
print(fit1, pars = c("alpha0", "gamma0", "sigma_beta0", "sigma_y", "ybar_new"))

# plot true vs estimated y_new
ff0 <- summary(fit0)
ff0 <- ff0$summary
beta0_new0 <- ff0[grep("^beta0", rownames(ff0)), "mean"]
beta0_true0 <- popdata0[["beta0"]]
beta1_new0 <- ff0[grep("^beta1", rownames(ff0)), "mean"]
beta1_true0 <- popdata0[["beta1"]]
y_new0 <- ff0[grep("y_new*", rownames(ff0)), "mean"]
y_true0 <- dplyr::summarise(group_by(popdata0[["pop.data"]], cluster.id),
                            ybar_true0 = mean(y))
y_true0 <- y_true0$ybar_true0
plot(y_new0, y_true0)
points(beta0_new0, beta0_true0, col = "red")
abline(0,1)
alpha0_new0 <- ff0["alpha0", "mean"]
gamma0_new0 <- ff0["gamma0", "mean"]
plot(beta0_true0, alpha0_new0 + gamma0_new0*popdata0[["logMj_c"]])
abline(0,1)
plot(beta0_new0, beta0_true0)
abline(0,1)
plot(beta1_new0, beta1_true0)
abline(0,1)
plot(beta0_true0, popdata0[["logMj_c"]])
abline(lm(beta0_true0 ~ popdata0[["logMj_c"]]))
abline(alpha0_new0,0)

beta0_true0 <- popdata0[["beta0"]]
beta0_new0 <- ff0[grep("^beta0*", rownames(ff0)), "mean"]
plot(beta0_new0, beta0_true0)
abline(0,1)


ff1 <- summary(fit1)
ff1 <- ff1$summary
y_new1 <- ff1[grep("y_new*", rownames(ff1)), "mean"]
y_true1 <- dplyr::summarise(group_by(popdata1[["pop.data"]], cluster.id),
                            ybar_true1 = mean(y))
y_true1 <- y_true1$ybar_true1
plot(y_new1, y_true1)
abline(0,1)
beta0_true1 <- popdata1[["beta0"]]
beta0_new1 <- ff1[grep("^beta0*", rownames(ff1)), "mean"]
plot(beta0_new1, beta0_true1)
abline(0,1)
alpha0_new1 <- ff1["alpha0", "mean"]
gamma0_new1 <- ff1["gamma0", "mean"]
points(beta0_new1, alpha0_new1 + gamma0_new1*popdata1[["logMj_c"]],
       col = "red")
