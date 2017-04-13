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


#############################################################################
### Load data 
#############################################################################

  # First load pop data
  all.data <- readRDS(paste0(rootdir, "/src/simulation/all_data_with_rng.rds"))
  pop.data <- all.data[["pop.data"]]
  truepars <- all.data[["true.pars.df"]]

  # Calculate true ybars
  ybar.true.nb <- mean(pop.data$y[pop.data$mod == "nb"])
  ybar.true.ln <- mean(pop.data$y[pop.data$mod == "ln"])

  # Loop through all sim files to load
  sim.files <- list.files(paste0(rootdir, "/src/simulation/tmp/"),
                          "ybar_df_sim_.*.rds")
  all.sims <- data.frame()
  for (j in 1:length(sim.files)) {
    curr.file <- readRDS(paste0(rootdir, "/src/simulation/tmp/", sim.files[j]))
    all.sims <- rbind(all.sims, curr.file)
  }

  all.sims$ybar_true <- ybar.true.nb
  all.sims$ybar_true[all.sims$model == "ln"] <- ybar.true.ln
  all.sims$model_str <- "Negative Binomial"
  all.sims$model_str[all.sims$model == "ln"] <- "Lognormal" 
  
  all.sims %>%
    dplyr::group_by(model_str, sim) %>%
    dplyr::summarise(ybar_true     = mean(ybar_true),
                     ybar_new_mean = mean(ybar_new),
                     ybar_lci_95   = quantile(ybar_new, 0.025),
                     ybar_lci_50   = quantile(ybar_new, 0.25),
                     ybar_uci_50   = quantile(ybar_new, 0.75),
                     ybar_uci_95   = quantile(ybar_new, 0.975)) -> sim.summ
    sim.summ %>% 
    dplyr::mutate(ybar_true_in_50  = ybar_true >= ybar_lci_50 &
                                       ybar_true <= ybar_uci_50,
                  ybar_true_in_95  = ybar_true >= ybar_lci_95 &
                                       ybar_true <= ybar_uci_95) %>%
    dplyr::summarise(ybar_covg_50  = mean(ybar_true_in_50),
                     ybar_covg_95  = mean(ybar_true_in_95)) -> covg
  print(covg)
  saveRDS(sim.summ, paste0(rootdir, "/src/simulation/sim_summ.rds"))

  sim.summ %>%
    dplyr::group_by(model_str) %>%
    dplyr::summarise(ybar_true = mean(ybar_true),
                     bias = mean(ybar_true - ybar_new_mean),
                     rel_bias = mean((ybar_true - ybar_new_mean) / ybar_true),
                     rmse = sqrt(mean((ybar_true - ybar_new_mean)^2)),
                     rel_rmse = sqrt(mean(((ybar_true - ybar_new_mean) / ybar_true)^2)),
                     avg_rel_cilen_95 = mean((ybar_uci_95 - ybar_lci_95) / ybar_true),
                     avg_rel_cilen_50 = mean((ybar_uci_50 - ybar_lci_50) / ybar_true)) -> sim.summ2
  print(sim.summ2)
  saveRDS(sim.summ2, paste0(rootdir, "/src/simulation/sim_summ2.rds"))

  p <- ggplot(sim.summ, aes(x = ybar_new_mean)) +
         geom_vline(aes(xintercept = ybar_true), colour = "grey50") +
         geom_line(stat = "density") +
         facet_wrap(~ model_str) +
         xlab("mean of ybar_new across sims") +
         theme_bw()
  ggsave(p, file = paste0(rootdir, "/src/simulation/ybar_density_of_means.pdf"),
         width = 10, height = 8)

  p2 <- ggplot(all.sims) +
          geom_line(aes(x = ybar_new, group = sim), stat = "density", colour = "grey50") +
          geom_vline(aes(xintercept = ybar_true)) +
          facet_wrap(~ model_str) +
          xlab("ybar_new") +
          ylab("density") +
          theme_bw()
  ggsave(p2, file = paste0(rootdir, "/src/simulation/ybar_density_of_sims.pdf"),
         width = 10, height = 8)


