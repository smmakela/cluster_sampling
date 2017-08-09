# Author: Susanna Makela
# Date: 21 Apr 2014
# Purpose: run stan

runstan <- function(num_clusters, num_units, use_sizes, outcome_type, size_model, rootdir,
                    simno, stanmod, stanmod_name, sim_data, num_iter, num_chains) {
  # num_clusters -- number of clusters to sample
  # num_units -- number of units to sample
  # use_sizes -- 0/1 for whether cluster sizes used in pop data
  # outcome_type -- whether outcome is continuous or binary
  # rootdir -- root directory where src, output, etc folders are
  # simno -- current iteration; used so that multiple instances aren't trying to write to the same file
  # stanmod -- compiled stan model
  # stanmod_name -- string for name of stan model so we know which parts of the code to run
  # sim_data -- data to use for simulation
  # num_iter -- number of iterations stan should run for
  # num_chains -- number of chains to run in stan

  ##########################################
  ### Load pop and sample data
  ##########################################
  rstan_options(auto_write = TRUE)
  options(mc_cores = parallel::detectCores())

  if (num_units <= 1) {
    nunits <- paste(100*num_units, "pct", sep = "")
  } else {
    nunits <- num_units
  }

  print("str(sim_data)")
  print(str(sim_data))

  for (j in names(sim_data)) {
    assign(j, sim_data[[j]])
  }
  rm(sim_data)
  ybar_true <- mean(pop_data$y)

  model_name <- gsub("_binary", "", stanmod_name)

  ##########################################
  ### Prep data for stan
  ##########################################
  print("making stan data")
  print(Sys.time())

  if (num_units != 999) {
    K <- num_clusters
  } else {
    K <- J
  }

  # SORT pop and sample data by cluster_id
  pop_data <- dplyr::arrange(pop_data, cluster_id)
  sample_data <- dplyr::arrange(sample_data, cluster_id)

  # Take things directly from sample_data
  x <- sample_data$x
  y <- sample_data$y
  cluster_id <- sample_data$cluster_id
  n <- nrow(sample_data)

  # Get pop data at cluster level -- calculate mean of x, grouping by cluster_id,
  # stratum_id, Mj, logMj_c (cluster_id is the most detailed)
  if (grepl("strat", stanmod_name)) {
    cluster_data_pop <- pop_data %>%
      dplyr::select(cluster_id, stratum_id, Mj, logMj_c, x) %>%
      dplyr::group_by(cluster_id, stratum_id, Mj, logMj_c) %>%
      dplyr::summarise(xbar_pop = mean(x)) %>%
      dplyr::arrange(cluster_id)
  } else {
    cluster_data_pop <- pop_data %>%
      dplyr::select(cluster_id, Mj, logMj_c, x) %>%
      dplyr::group_by(cluster_id, Mj, logMj_c) %>%
      dplyr::summarise(xbar_pop = mean(x)) %>%
      dplyr::arrange(cluster_id)
  }
print("str(cluster_data_pop)")
print(str(cluster_data_pop))
  Nj_pop <- cluster_data_pop$Mj
  N <- sum(Nj_pop)
  xbar_pop <- cluster_data_pop$xbar_pop

  # Get sample data at cluster level
  if (grepl("strat", stanmod_name)) {
    stratum_id <- sampled_cluster_data$stratum_id
  }
  sampled_cluster_data <- dplyr::arrange(sampled_cluster_data, cluster_id)
  Nj_sample <- sampled_cluster_data$Mj
  log_Nj_sample <- sampled_cluster_data$logMj_c

  # special data for bb
  n_dat <- sampled_cluster_data %>%
    dplyr::arrange(cluster_id) %>%
    dplyr::group_by(cluster_id, Mj) %>%
    dplyr::summarise(n = n_distinct(cluster_id))
  M <- nrow(n_dat)      # number of unique cluster sizes
  M_counts <- n_dat$n   # counts of unique cluster sizes
  Nj_unique <- n_dat$Mj # vector of unique cluster sizes

  ##########################################
  ### Make inputs to stan -- data list, parameters, functions
  ##########################################
  source(paste0(rootdir, "src/simulation/make_stan_data.R"))

  ##########################################
  ### Run stan
  ##########################################
  print(Sys.time())

  fit <- sampling(stanmod, data = standata,
                  iter = num_iter, chains = num_chains,
                  control = list(stepsize = 0.001, adapt_delta = 0.999))
  print("done fitting stan model")
  print(Sys.time())
  print(warnings())

  ##########################################
  ### See if there are any divergent transitions -- if so need to rerun 
  ##########################################
  source(paste0(rootdir, "src/simulation/check_for_div_trans.R"))

  # Save the current fit
  if (simno == 1) {
    saveRDS(fit, file = paste0(rootdir, "/output/simulation/stanfit_usesizes_",
                               use_sizes, "_", outcome_type, "_", size_model,
                                "_", model_name, "_nclusters_", num_clusters,
                               "_nunits_", nunits, "_sim_", simno,".rds"))
    print("done saving stanfit object")
    cat("filename:", paste0(rootdir, "/output/simulation/stanfit_usesizes_",
                               use_sizes, "_", outcome_type, "_", size_model,
                                "_", model_name, "_nclusters_", num_clusters,
                               "_nunits_", nunits, "_sim_", simno,".rds"))
    print(Sys.time())
  }

  # Print summary for the basic parameters
  print(fit, pars = parlist)

  ##########################################
  ### Extract samples, format param estimates to just keep necessary info 
  ##########################################
  source(paste0(rootdir, "src/simulation/process_fit.R"))

  ##########################################
  ### Draw new cluster sizes (if needed) and ybar_new
  ##########################################
  source(paste0(rootdir, "src/simulation/make_pp_draws.R"))

  print("done drawing Nj_new, ybar_new")
  print(Sys.time())      

  yn_list <- list(ybar_new = ybar_new, Nj_new_df = Nj_new_df)
  saveRDS(yn_list,
          paste0(rootdir, "output/simulation/yn_list_usesizes_",
                 use_sizes, "_", outcome_type, "_", size_model, "_", model_name,
                 "_nclusters_", num_clusters, "_nunits_", nunits,
                 "_sim_", simno, ".rds"))

  ##########################################
  ### Summarise Nj_new, ybar_new and add to parameter estimates
  ##########################################
  print("summarizing draws")
  if (grepl("cluster_inds_only", stanmod_name) ||
      grepl("knowsizes", stanmod_name)) {
    Nj_new_means <- NA
    tmp <- data.frame(ybar_new, draw_num = c(1:num_draws))
    tmp %>%
      dplyr::summarise(mean = mean(ybar_new),
                       sd = sd(ybar_new),
                       p025 = quantile(ybar_new, 0.025),
                       p25 = quantile(ybar_new, 0.25),
                       p50 = quantile(ybar_new, 0.50),
                       p75 = quantile(ybar_new, 0.75),
                       p975 = quantile(ybar_new, 0.975)) -> draw_summ
    draw_summ$truth <- ybar_true
  } else {
    Nj_new_df %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(Nj_new = mean(Nj_new)) -> Nj_new_means
    Nj_new_means$truth <- Nj_pop

    true_vals <- data.frame(param = c("N_new", "ybar_new"),
                            truth = c(N, ybar_true))
    tmp <- data.frame(N_new, ybar_new, draw_num = c(1:num_draws))
    tmp %>%
      tidyr::gather(key = param, value = value, -draw_num) %>%
      dplyr::group_by(param) %>%
      dplyr::summarise(mean = mean(value),
                       sd = sd(value),
                       p025 = quantile(value, 0.025),
                       p25 = quantile(value, 0.25),
                       p50 = quantile(value, 0.50),
                       p75 = quantile(value, 0.75),
                       p975 = quantile(value, 0.975)) %>%
      dplyr::left_join(., true_vals, by = "param") -> draw_summ
  }
  print(Sys.time())      
  print(draw_summ)
  print(warnings())

  ##########################################
  ### Save results
  ##########################################
  results_list <- list(par_ests = par_ests, Nj_new_means = Nj_new_means,
                       draw_summ = draw_summ, Nj_sample = Nj_sample)
  saveRDS(results_list,
          paste0(rootdir, "output/simulation/stan_results_usesizes_",
                 use_sizes, "_", outcome_type, "_", size_model, "_", model_name,
                 "_nclusters_", num_clusters, "_nunits_", nunits,
                 "_sim_", simno, ".rds"))
  cat("results saved to:", 
      paste0(rootdir, "output/simulation/stan_results_usesizes_",
             use_sizes, "_", outcome_type, "_", size_model, "_", model_name,
             "_nclusters_", num_clusters, "_nunits_", nunits,
             "_sim_", simno, ".rds"), "\n")

  return(NULL)
}

