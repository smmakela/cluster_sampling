args <- commandArgs(FALSE)
print(args)
len <- length(args)
sim <- -1*as.numeric(args[len])
print(sim)

#############################################################################
### Set lib paths, source files 
#############################################################################
  libdir <- "/vega/stats/users/smm2253/rpackages"
  .libPaths(libdir)
  rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
  Sys.setenv(HOME = rootdir)

# Set up libraries, directories, options, etc
library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

codedir <- paste0(rootdir, "/src/simulation/") 

# Load stan models
negbin_file <- paste0(codedir, "negbin_full_centered.stan")
nb_stan_cp <- stan_model(file = negbin_file)
expose_stan_functions(nb_stan_cp)
lognormal_file <- paste0(codedir, "lognormal_full_centered.stan")
ln_stan_cp <- stan_model(file = lognormal_file)
expose_stan_functions(ln_stan_cp)


source(paste0(codedir, "rspps.r"))


all.data <- readRDS(paste0(codedir, "/all_data_with_rng.rds"))
print("str(all.data)")
print(str(all.data))
pop.data <- all.data[["pop.data"]]
true.pars.df <- all.data[["true.pars.df"]]
print("TRUE PARS:")
print(str(true.pars.df))

J <- 1000
K <- 50

# Sample K of the J clusters and 10% of units in each cluster
ptype <- "cp"
ybar.df <- data.frame()
sample.data <- data.frame()
xbar.pop.data <- data.frame()
for (sfx in c("nb", "ln")) {
  Nj.pop.dat <- unique(pop.data[pop.data$mod == sfx, c("cluster.id", "Nj.pop")])
  print("Nj.pop.dat:")
  print(str(Nj.pop.dat))

  # Sample clusters PPS
  sampled.clusters <- rspps(Nj.pop.dat$Nj.pop, c(1:J), K)
  
  # Create map to renumber sampled clusters 1:K and the rest (K+1):J
  unsampled.clusters <- sort(setdiff(c(1:J), sampled.clusters))
  id.map <- data.frame(cluster.id = c(sampled.clusters, unsampled.clusters),
                       new.cluster.id = c(1:J))
  
  # Create sampled data
  pop.data %>%
    dplyr::filter(mod == sfx) %>%
    # merge id map to get renumbered cluster ids
    left_join(., id.map, by = "cluster.id") %>%
    dplyr::rename(orig.cluster.id = cluster.id,
                  cluster.id = new.cluster.id) %>%
    # only keep sampled clusters
    dplyr::filter(orig.cluster.id %in% sampled.clusters) %>%
    # get random sample 10% of each cluster
    dplyr::group_by(cluster.id) %>%
    dplyr::sample_frac(0.1) -> sampdat
  sample.data <- rbind(data.frame(sample.data), data.frame(sampdat))
  
  # Calculate cluster-level info (x_bar)
  pop.data %>%
    dplyr::filter(mod == sfx) %>%
    dplyr::group_by(cluster.id, mod) %>%
    dplyr::summarise(xbar = mean(x)) -> xbar_pop
  xbar.pop.data <- rbind(data.frame(xbar.pop.data), data.frame(xbar_pop))
}  

for (sfx in c("nb", "ln")) {
  cat("\n========================================================================\n")
  cat("MODEL:", sfx, ", ", "simno:", sim, "\n")
  parnames <- c("alpha0", "gamma0", "alpha1", "gamma1",
                "sigma_beta0", "sigma_beta1", "sigma_y")
  parnames2 <- c("beta0", "beta1")
  if (sfx == "nb") {
    parnames <- c(parnames, "mu", "phi")
  } else {
    parnames <- c(parnames, "mu", "sigma")
  }
  
  # Pull out sampled data and data for sampled clusters
  sample.data %>%
    dplyr::filter(mod == sfx) -> currdat
  currdat %>%
    dplyr::distinct(cluster.id, Nj.pop, log.Nj.pop.c) -> currdat.cluster
  xbar.pop.data %>%
    dplyr::filter(mod == sfx) -> curr.pop.dat
  
  # Create stan data and run stan
  stan_data <- list(J = J, K = K, n = nrow(currdat),
                    x = currdat$x, y = currdat$y,
                    cluster_id = currdat$cluster.id,
                    Nj_sample = currdat.cluster$Nj.pop,
                    log_Nj_sample = currdat.cluster$log.Nj.pop.c)
  stan_res <- sampling(get(paste0(sfx, "_stan_", ptype)), data = stan_data,
                       iter = 2000, chains = 4,
                       control = list(adapt_delta = 0.99))
  print(stan_res, pars = parnames)
  
  # Extract samples and join with true param values
  stan_samps_wide <- data.frame(rstan::extract(stan_res, pars = parnames))
  stan_samps_long <- tidyr::gather(stan_samps_wide, key = parnames, value = value)
  stan_samps_long$parnames <- as.character(stan_samps_long$parnames)
  if (sfx == "nb") {
    stan_samps_long$parnames[stan_samps_long$parnames == "mu"] <- "mu_nb" 
  }
print("joining stan_samps and true pars")
print(str(stan_samps_long))
print(str(true.pars.df))
  stan_samps <- left_join(stan_samps_long, true.pars.df, by = "parnames")
  
  # Make density plots of posterior distributions with true parameter values
  #p <- ggplot(stan_samps, aes(x = value)) +
  #  geom_line(stat = "density") +
  #  geom_vline(aes(xintercept = true.pars)) +
  #  facet_wrap(~ parnames, scales = "free") +
  #  ggtitle(paste0("Model: ", sfx, ", ", ptype)) +
  #  theme_bw()
  #print(p)
  
  # Get beta's
  stan_samps_betas <- data.frame(rstan::extract(stan_res, pars = parnames2))
  stan_samps_betas %>%
    tidyr::gather(key = pname, value = value) %>%
    dplyr::mutate(par = substr(pname, 1, 5),
                  cluster.id = as.numeric(substr(pname, 7, nchar(pname))),
                  pname = NULL) %>%
    dplyr::group_by(par, cluster.id) %>%
    dplyr::mutate(draw.num = row_number()) %>%
    tidyr::spread(key = par, value = value) -> beta_samps
  
  # Merge beta's with xbar_pop -- will only be for *sampled* clusters
print("joining beta_samps and curr pop dat")
  beta_samps <- left_join(beta_samps, curr.pop.dat, by = "cluster.id")
  
  # Calculate generated quantities here
  # (hard to do in stan directly since get errors during warmup)
  # For sampled clusters, we can just use the posterior draws of beta0 and beta1 directly
  #beta_samps$yj_new <- beta_samps$beta0 + beta_samps$beta1 * beta_samps$xbar
  
  # For unsampled clusters, we need to first draw Nj_new, then draw beta0_new, beta1_new,
  # and *then* draw yj_new
  num_draws <- nrow(stan_samps_wide)
  ybar_new <- rep(NA, times = num_draws)
  if (sfx == "nb") {
    for (s in 1:num_draws) {
      beta_samps_sub <- dplyr::filter(beta_samps, draw.num == s)
      ybar_new[s] <- ybar_new_nb_rng(J, K, beta_samps_sub$beta0, beta_samps_sub$beta1,
                                     currdat.cluster$Nj.pop, curr.pop.dat$xbar,
                                     stan_samps_wide$mu[s], stan_samps_wide$phi[s],
                                     stan_samps_wide$alpha0[s], stan_samps_wide$gamma0[s],
                                     stan_samps_wide$sigma_beta0[s],
                                     stan_samps_wide$alpha1[s], stan_samps_wide$gamma1[s],
                                     stan_samps_wide$sigma_beta1[s], stan_samps_wide$sigma_y[s])
    }
  } else {
    for (s in 1:num_draws) {
      beta_samps_sub <- dplyr::filter(beta_samps, draw.num == s)
      ybar_new[s] <- ybar_new_ln_rng(J, K, beta_samps_sub$beta0, beta_samps_sub$beta1,
                                     currdat.cluster$Nj.pop, curr.pop.dat$xbar,
                                     stan_samps_wide$mu[s], stan_samps_wide$sigma[s],
                                     stan_samps_wide$alpha0[s], stan_samps_wide$gamma0[s],
                                     stan_samps_wide$sigma_beta0[s],
                                     stan_samps_wide$alpha1[s], stan_samps_wide$gamma1[s],
                                     stan_samps_wide$sigma_beta1[s], stan_samps_wide$sigma_y[s])
    } 
  }
  
  # Plot posterior dist of ybar_new
  ybar_new <- data.frame(ybar_new)
  # p2 <- ggplot(ybar_new, aes(x = ybar_new) ) +
  #   geom_line(stat = "density") +
  #   geom_vline(xintercept = mean(pop.data$y[pop.data$mod == sfx])) +
  #   ggtitle(paste0("Model: ", sfx, ", ", ptype)) +
  #   theme_bw()
  # print(p2)
  
  tmp <- data.frame(ybar_new, sim = rep(sim, nrow(ybar_new)),
                    model = rep(sfx, nrow(ybar_new)))
  ybar.df <- rbind(ybar.df, tmp)
  
} # end negbin/lognormal loop
  
saveRDS(ybar.df, paste0(codedir, "/tmp/ybar_df_sim_", sim, ".rds"))

