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
  ### Source functions
  ##########################################
  source(paste0(rootdir, "src/simulation/make_stan_data.R"))
  source(paste0(rootdir, "src/simulation/check_for_div_trans.R"))
  source(paste0(rootdir, "src/simulation/process_fit.R"))
  source(paste0(rootdir, "src/simulation/make_pp_draws.R"))
  source(paste0(rootdir, "src/simulation/make_pp_draw_plots.R"))

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

  print("str(sim_data, give.attr = FALSE)")
  print(str(sim_data, give.attr = FALSE))

  # Add elements of sim_data to current enviroment
  for (j in names(sim_data)) {
    assign(j, sim_data[[j]])
  }
  rm(sim_data)
  ybar_true <- truepars$ybar_true

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

  ##########################################
  ### Make inputs to stan -- data list, parameters, functions
  ##########################################
  res <-  make_stan_data(stanmod_name, sample_data, pop_cluster_data,
                         sampled_cluster_data) 
  for (j in names(res)) {
    assign(j, res[[j]])
  }

  print("parlist post make stan data:")
  print(parlist)
  #source(paste0(rootdir, "src/simulation/make_stan_data.R"))

  ##########################################
  ### Run stan
  ##########################################
  print(Sys.time())
print(str(standata, give.attr = FALSE))
print(str(standata$Nj_sample, give.attr = FALSE))
  fit <- sampling(stanmod, data = standata,
                  iter = num_iter, chains = num_chains,
                  control = list(stepsize = 0.001, adapt_delta = 0.999))
  print("done fitting stan model")
  print(Sys.time())
  print(warnings())
  print("parlist post fitting stan model:")
  print(parlist)

  ##########################################
  ### See if there are any divergent transitions -- if so need to rerun 
  ##########################################
  fit <- check_for_div_trans(num_clusters, num_units, use_sizes, outcome_type,
                             size_model, rootdir, simno, stanmod, standata,
                             num_iter, num_chains, fit)
  print("parlist post check div trans:")
  print(parlist)

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
  print("parlist post save:")
  print(parlist)
  print(fit, pars = parlist)

  ##########################################
  ### Extract samples, format param estimates to just keep necessary info 
  ##########################################
  if (!exists("kappapars")) {
    kappapars <- NA
  }
  if (!exists("sbpars")) {
    kappapars <- NA
  }
  res <- process_fit(fit, parlist, betapars, kappapars, sbpars, stanmod_name) 
  for (j in names(res)) {
    assign(j, res[[j]])
  }
  rm(res)

  ##########################################
  ### Draw new cluster sizes (if needed) and ybar_new
  ##########################################
  if (!exists("phi_star_samps")) {
    phi_star_samps <- NA
  }
  if (!exists("sb_samps")) {
    sb_samp <- NA
  }
print("STANMOD_NAME:")
print(stanmod_name)
print("str(param_samps, give.attr = FALSE):")
print(str(param_samps, give.attr = FALSE))
  res <- make_pp_draws(stanmod, stanmod_name, param_samps, beta_samps,
                       phi_star_samps, sb_samps, pop_cluster_data,
                       sampled_cluster_data) 
  for (j in names(res)) {
    assign(j, res[[j]])
  }
  rm(res)
print("str(Mj_new_df, give.attr = FALSE):")
print(str(Mj_new_df, give.attr = FALSE))

  print("done drawing Mj_new, ybar_new")
  print(Sys.time())      

  yn_list <- list(ybar_new = ybar_new, Mj_new_df = Mj_new_df)
  saveRDS(yn_list,
          paste0(rootdir, "output/simulation/yn_list_usesizes_",
                 use_sizes, "_", outcome_type, "_", size_model, "_", model_name,
                 "_nclusters_", num_clusters, "_nunits_", nunits,
                 "_sim_", simno, ".rds"))

  ##########################################
  ### Plot Mj_new, ybar_new
  ##########################################
  if (simno == 1) {
    make_pp_draw_plots(Mj_new_df, ybar_new, ybar_true,
                       num_clusters, num_units, use_sizes, outcome_type,
                       size_model, rootdir, stanmod_name,
                       pop_cluster_data, sampled_cluster_data) 
  }

  ##########################################
  ### Summarise Mj_new, ybar_new and add to parameter estimates
  ##########################################
  print("summarizing draws")
  num_draws <- nrow(param_samps)
  if (grepl("cluster_inds_only", stanmod_name) ||
      grepl("knowsizes", stanmod_name)) {
    Mj_new_means <- NA
    tmp <- data.frame(ybar_new, draw_num = c(1:num_draws))
    draw_summ <- tmp %>%
      dplyr::summarise(mean = mean(ybar_new),
                       sd = sd(ybar_new),
                       p025 = quantile(ybar_new, 0.025),
                       p25 = quantile(ybar_new, 0.25),
                       p50 = quantile(ybar_new, 0.50),
                       p75 = quantile(ybar_new, 0.75),
                       p975 = quantile(ybar_new, 0.975)) 
    draw_summ$truth <- ybar_true
  } else {
    Mj_new_means <- Mj_new_df %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(Mj_new = mean(Mj_new)) 
    Mj_new_means$truth <- pop_cluster_data$Mj

    true_vals <- data.frame(param = c("M_tot_new", "ybar_new"),
                            truth = c(sum(pop_cluster_data$Mj), ybar_true))
    tmp <- data.frame(M_tot_new, ybar_new, draw_num = c(1:num_draws))
    draw_summ <- tmp %>%
      tidyr::gather(key = param, value = value, -draw_num) %>%
      dplyr::group_by(param) %>%
      dplyr::summarise(mean = mean(value),
                       sd = sd(value),
                       p025 = quantile(value, 0.025),
                       p25 = quantile(value, 0.25),
                       p50 = quantile(value, 0.50),
                       p75 = quantile(value, 0.75),
                       p975 = quantile(value, 0.975)) %>%
      dplyr::left_join(., true_vals, by = "param") 
  }
  print(Sys.time())      
  print(draw_summ)
  print(warnings())

  ##########################################
  ### Save results
  ##########################################
  results_list <- list(par_ests = par_ests, Mj_new_means = Mj_new_means,
                       draw_summ = draw_summ, Mj_sample = sampled_cluster_data$Mj)
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

