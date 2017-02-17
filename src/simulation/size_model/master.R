#!/usr/bin/R Rscript
libdir <- "/vega/stats/users/smm2253/rpackages"
.libPaths(libdir)
rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
Sys.setenv(HOME = rootdir)
codedir <- paste0(rootdir, "src/simulation/size_model/")
resdir <- paste0(rootdir, "output/simulation/size_model/")

library(rstan)
library(ggplot2)
library(tidyr)
library(dplyr)
source(paste0(codedir, "sim.R"))
source(paste0(codedir, "rspps.R"))

J <- 100
K.list <- c(5, 10, 25)
n.iter <- 1000
n.chains <- 4
n.sims <- 100
n.reps <- 100
p.vals.df <- data.frame()
model.name.list <- c("negbin", "lognormal")
#pop.mean <- 300
#pop.sd <- 50
#k <- pop.mean^2 / (
params.list <- list(negbin = list(k = 200, p = 0.1),
                    lognormal = list(mu = 7.5, sigma = 0.1))

# If needed, compile stan models
for (j in 1:length(model.name.list)) {
  stanmod_name <- model.name.list[j]
  code.file <- paste0(codedir, stanmod_name, ".stan")
  code.mtime <- file.mtime(code.file)
  rdat.file <- paste0(codedir, stanmod_name, ".RData")
  rdat.mtime <- file.mtime(rdat.file)
  if (code.mtime > rdat.mtime) {
    stanmod <- stan_model(file = code.file)
    save(stanmod, file = rdat.file)
  }
}

for (model.name in model.name.list) {
  load(paste0(codedir, model.name, ".RData"))
  params <- params.list[[model.name]]
  # Generate population data
  chk <- rep(1, times = J)
  while (sum(chk >= 1/max(K.list)) > 0) {
    if (model.name == "negbin") {
      k <- params[["k"]]
      p <- params[["p"]]
      Nj.pop <- rnbinom(n = J, size = k, prob = p)
    } else {
      mu <- params[["mu"]]
      sigma <- params[["sigma"]]
      Nj.pop <- exp(rnorm(n = J, mean = mu, sd = sigma)) # log(Nj.pop) ~ N(0,1)
    }
    chk <- Nj.pop / sum(Nj.pop)
  }
  pop.data <- data.frame(ids = c(1:J), Nj.pop)
  
  for (K in K.list) {
    for (j in 1:n.reps) {
      cat("On rep", j, "of", n.reps, "for", model.name, "K =", K, "\n")
      res <- sim(J, K, pop.data, params, model.name, stanmod,
                 n.iter, n.chains, n.sims)
      p.vals <- res[["p.vals"]]
      p.vals$rep <- j
      p.vals$model.name <- model.name
      p.vals$K <- K
      p.vals.df <- rbind(p.vals.df, p.vals)
      if (j == 1) { # save the full results (incl stan output) for just the first rep
        saveRDS(res, paste0(resdir, "sim1_", model.name, "_K_", K, ".rds"))
      }
    } # end j loop
  } # end k loop
}
saveRDS(p.vals.df, file = paste0(resdir, "sim_res.rds"))

# density plot of p-values
pp <- ggplot(p.vals.df, aes(x = p_value)) +
  geom_line(aes(colour = factor(K)), stat = "density") +
  facet_grid(model.name ~ stat) +
  theme_bw()
ggsave(filename = paste0(resdir, "/pvalue_densities.pdf"),
       plot = pp, width = 10, height = 6)

# compile results
summary.df <- data.frame()
Nj.sims <- data.frame()
pop.data <- data.frame()
for (m in c("negbin", "lognormal")) {
  for (K in c(5, 10, 25)) {
    ff <- readRDS(paste0(resdir, "sim1_", m, "_K_", K, ".rds"))
    tmp <- data.frame(ff[["summary.df"]], model.name = ff[["model.name"]], K = K)
    summary.df <- rbind(summary.df, tmp)
    tmp <- data.frame(ff[["Nj.sims"]], model.name = ff[["model.name"]], K = K)
    Nj.sims <- rbind(Nj.sims, tmp)
    if (K == 5) {
      tmp <- data.frame(Nj.pop = ff[["Nj.pop"]], model.name = ff[["model.name"]])
      pop.data <- rbind(pop.data, tmp)
    }
  }
}
pp <- ggplot(summary.df, aes(x = value, colour = factor(K))) +
  geom_line(stat = "density") +
  geom_vline(aes(xintercept = truth, colour = factor(K))) +
  facet_wrap(model.name ~ stat, scales = "free", nrow = 2) +
  theme_bw()
ggsave(filename = paste0(resdir, "/test_stat_densities.pdf"),
       plot = pp, width = 10, height = 6)

pp <- ggplot(Nj.sims[Nj.sims$simno != 0, ], aes(x = Nj)) +
  geom_line(stat = "density", aes(group = simno, colour = "Posterior samples")) +
  geom_line(data = pop.data, aes(x = Nj.pop, colour = "Population sizes"),
            stat = "density", size = 1.5) +
  geom_line(data = Nj.sims[Nj.sims$simno == 0, ],
            aes(x = Nj, colour = "Missing sizes"), stat = "density", size = 1.5) +
  facet_wrap(~ model.name, scales = "free") +
  scale_colour_manual("",
                      values = c("Posterior samples" = "grey80",
                                 "Population sizes" = "grey50",
                                 "Missing sizes" = "black")) +
  xlab("Cluster size") +
  ylab("Density") +
  theme_bw()
ggsave(filename = paste0(resdir, "/size_densities.pdf"),
       plot = pp, width = 10, height = 6)

