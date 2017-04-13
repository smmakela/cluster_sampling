draw_pop_cluster_sizes_for_sim <- function(J) {
  # Purpose: use the observed city sizes in the original FF sampling frame
  # to simulate J cluster sizes for the simulation by fitting a gamma
  # distribution to the log of the observed city sizes, and then drawing J
  # log cluster sizes from gamma(alpha_hat, beta_hat)
  #
  # Inputs:
  #   J -- number of cluster sizes to draw
  #
  # Outputs:
  #   pop_for_sim -- vector of J cluster sizes drawn from gamma(alpha_hat, beta_hat)

  #############################################################################
  ### Set lib paths, source files 
  #############################################################################
    libdir <- "/vega/stats/users/smm2253/rpackages"
    .libPaths(libdir)
    rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
    Sys.setenv(HOME = rootdir)
  
    # Set up libraries, directories, options, etc
    library(rstan)
    
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
    
    codedir <- paste0(rootdir, "/src/simulation/") 
    figdir <- paste0(rootdir, "/output/figures/")
  
  #############################################################################
  ### Load data and fit stan model
  #############################################################################
    pop_dat <- read.csv(paste0(codedir, "/observed_city_pops.csv"), header = TRUE)
    pops <- pop_dat$population
  
    pop_stan_dat <- list(N = length(pops), pop = pops)
  
    pop_model_code <- '
      data {
        int N; # number of cities
        vector[N] pop; # pop vector
      }
      transformed data {
        vector[N] logpop;
        logpop = log(pop);
      }
      parameters {
        real<lower=0> alpha;
        real<lower=0> beta;
      }
      model {
        logpop ~ gamma(alpha, beta);
      } '
    # ' adding this to make the highlighting nice
  
    pop_mod <- stan_model(model_code = pop_model_code)
    pop_res <- sampling(pop_mod, data = pop_stan_dat, iter = 5000, chains = 4)
  
  #############################################################################
  ### Inspect results and simulate
  #############################################################################
    # Print summary of stan results
    print(pop_res)
  
    # Extract post-warmup samples, retaining the chains
    samps <- as.array(extract(pop_res, permute = FALSE))
  
    # Plot traceplot of chains to check for lack of convergence
    p <- mcmc_trace(samps, facet_args = list(ncol = 1, strip.position = "left"),
                    pars = c("alpha", "beta"))
    ggsave(plot = p, width = 10, height = 8,
           file = paste0(rootdir, "/output/figures/obs_pop_pars_traceplot.pdf"))
  
    # Draw pop sizes to be used for simulation using posterior means of alpha
    # and beta
    alpha_hat <- mean(samps[,,"alpha"])
    beta_hat <- mean(samps[,,"beta"])
    log_pop_for_sim <- rgamma(n = J, shape = alpha_hat, scale = 1/beta_hat)
    pop_for_sim <- exp(log_pop_for_sim)
  
    # Plot densities of observed and simulated log pop sizes
    pdf(file = paste0(figdir, "/obs_pop_sim_pop_densities.pdf"),
        width = 10, height = 8)
    plot(density(log(pops), xlab = "Log pop size", ylab = "Density",
         main = paste0("Number of simulated clusters: ", J)))
    lines(density(log_pop_for_sim), lty = 2)
    legend(x = "topright", legend = c("Obs", "Sim"), lty = c(1, 2))
    dev.off()
    
    return(pop_for_sim)
}

