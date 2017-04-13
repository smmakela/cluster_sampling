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


# Parameters for population
J <- 1000    # number of clusters in pop
K <- 50      # number of clusters to sample

# params for negative binomial (mu_nb, phi are for the Stan parameterization)
k <- 200
p <- 0.2
mu_nb <- k * (1-p) / p
phi <- k

# params for lognormal (select so that get roughly same mean, var as above)
alpha <- mu_nb # mean of negative binomial
beta <- mu_nb + mu_nb^2 / phi # var of negative binomial
# eqns to back-solve for mu and sigma based on alpha and beta
sigma <- sqrt(log(exp(log(beta) - 2*log(alpha)) + 1))
mu <- 0.5*(2*log(alpha) - sigma^2)

# Create superpopulations
Nj.pop.nb <- rnbinom(n = J, size = k, prob = p)
Nj.pop.ln <- exp(rnorm(n = J, mean = mu, sd = sigma))

Nj.pop.df <- data.frame(Nj = c(Nj.pop.nb, Nj.pop.ln),
                        type = rep(c("negbin", "lognorm"), each = J))


# Generate population data frame
pop.data.nb <- data.frame(ids = c(1:J), Nj.pop.nb)
pop.data.ln <- data.frame(ids = c(1:J), Nj.pop.ln)

# Generate population data
alpha0 <- rnorm(1)
gamma0 <- rnorm(1)
alpha1 <- rnorm(1)
gamma1 <- rnorm(1)
sigma_beta0 <- abs(rnorm(1, 0, 0.5))
sigma_beta1 <- abs(rnorm(1, 0, 0.5))
sigma_y <- abs(rnorm(1, 0, 0.5))
true.pars <- c(mu_nb, phi, mu, sigma, alpha0, gamma0, alpha1, gamma1,
               sigma_beta0, sigma_beta1, sigma_y)
names(true.pars) <- c("mu_nb", "phi", "mu", "sigma",
                      "alpha0", "gamma0", "alpha1", "gamma1",
                      "sigma_beta0", "sigma_beta1", "sigma_y")
true.pars.df <- data.frame(true.pars)
true.pars.df$parnames <- as.character(rownames(true.pars.df))
pop.data <- data.frame()
for (sfx in c("nb", "ln")) {
  Nj.pop <- get(paste0("Nj.pop.", sfx))
  Nj.pop.rnd <- round(Nj.pop) # needed for ln; doesn't change anything for nb
  log.Nj.pop.c <- log(Nj.pop) - mean(log(Nj.pop))
  x <- sample(c(15:49), size = sum(Nj.pop.rnd), replace = TRUE)
  beta0 <- rnorm(n = J, mean = alpha0 + gamma0 * log.Nj.pop.c, sd = sigma_beta0)
  beta0_rep <- rep(beta0, Nj.pop.rnd)
  beta1 <- rnorm(n = J, mean = alpha1 + gamma1 * log.Nj.pop.c, sd = sigma_beta1)
  beta1_rep <- rep(beta1, Nj.pop.rnd)
  y <- rnorm(sum(Nj.pop.rnd), mean = beta0_rep + beta1_rep * x, sd = sigma_y)
  df <- data.frame(cluster.id = rep(c(1:J), times = Nj.pop.rnd),
                   Nj.pop = rep(Nj.pop, Nj.pop.rnd),
                   log.Nj.pop.c = rep(log.Nj.pop.c, Nj.pop.rnd),
                   unit.id = unlist(lapply(Nj.pop.rnd, seq_len)),
                   x, y)
  df$mod <- sfx
  pop.data <- rbind(pop.data, df)
}

all.data <- list(pop.data = pop.data, true.pars.df = true.pars.df)

saveRDS(all.data, paste0(codedir, "/all_data_with_rng.rds"))

