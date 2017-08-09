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
  ### Make data for stan
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

  #Nj_pop <- Mj
  #N <- sum(Nj_pop)
  #n <- sum(pop_data$insample)
  #tmp <- dplyr::summarise(group_by(pop_data, cluster_id),
  #                        xbar_pop = mean(x))
  #tmp <- dplyr::arrange(tmp, cluster_id)
  #xbar_pop <- tmp$xbar_pop
  #x <- sample_data$x
  #y <- sample_data$y
  #cluster_id <- sample_data$cluster_id

  # Get pop data at cluster level -- calculate mean of x, grouping by cluster_id,
  # stratum_id, Mj, logMj_c (cluster_id is the most detailed)
  #if (size_model == "ffstrat") {
print("ONE")
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
print("TWO")
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

  # delete pop and sample data to save memory
  #rm(pop_data)
  #rm(sample_data)
  ##########################################
  ### Make inputs to stan -- data list, parameters, functions
  ##########################################
  expose_stan_functions(stanmod)
  parlist <- c("alpha0", "gamma0", "alpha1", "gamma1",
               "sigma_beta0", "sigma_beta1", "sigma_y")
  betapars <- c("beta0", "beta1")
print("THREE")
  if (grepl("strat", stanmod_name)) {
    kappapars <- c("kappa0", "kappa1")
    if (grepl("lognormal", stanmod_name)) {
      sbpars <- c("mu", "phi")
    }
    if (grepl("negbin", stanmod_name)) {
      sbpars <- c("mu", "phi")
    }
  }
  if (grepl("bb", stanmod_name)) {
    standata <- list(J = J,
                     K = K,
                     n = n,
                     N = N,
                     M = M,
                     x = x,
                     y = y,
                     cluster_id = cluster_id,
                     Nj_sample = Nj_sample,
                     log_Nj_sample = log_Nj_sample,
                     Nj_unique = Nj_unique,
                     M_counts = M_counts)
  } else if (grepl("cluster_inds_only", stanmod_name)) {
    standata <- list(J = J,
                     K = K,
                     n = n,
                     x = x,
                     y = y,
                     cluster_id = cluster_id)
    parlist <- c("alpha0", "alpha1", "sigma_beta0", "sigma_beta1", "sigma_y")
  } else if (grepl("knowsizes", stanmod_name)) {
    standata <- list(J = J,
                     K = K,
                     n = n,
                     x = x,
                     y = y,
                     cluster_id = cluster_id,
                     Nj_sample = Nj_sample,
                     log_Nj_sample = log_Nj_sample)
  } else if (grepl("lognormal", stanmod_name)) {
    standata <- list(J = J,
                     K = K,
                     n = n,
                     x = x, 
                     y = y,
                     cluster_id = cluster_id,
                     Nj_sample = Nj_sample,
                     log_Nj_sample = log_Nj_sample)
    if (!grepl("strat", stanmod_name)) {
      parlist <- c(parlist, "mu_star", "sigma", "mu", "mu_star_scaled", "sigma_scaled")
    }
  } else if (grepl("negbin", stanmod_name)) {
print("FOUR")
    standata <- list(J = J,
                     K = K,
                     n = n,
                     x = x, 
                     y = y,
                     cluster_id = cluster_id,
                     Nj_sample = Nj_sample,
                     log_Nj_sample = log_Nj_sample)
    #if (size_model != "ffstrat") {
    if (!grepl("strat", stanmod_name)) {
      parlist <- c(parlist, "mu", "phi")
    }
  } else {
    stop("Invalid stan model")
  }
  #if (size_model == "ffstrat") {
  if (grepl("strat", stanmod_name)) {
print("FIVE")
    stratum_matrix <- model.matrix(~ 0 + factor(stratum_id), sampled_cluster_data)
    stratum_matrix_pop <- model.matrix(~ 0 + factor(stratum_id), cluster_data_pop)
    standata <- c(standata, list(stratum_id = stratum_id, S = 2,
                                 stratum_matrix = stratum_matrix))
    #parlist <- c(parlist, "mu", "phi")
  } #else { # don't include these for ff because they're vectors in it
    #parlist <- c(parlist, "mu_star", "phi_star", "mu", "phi")
  #}

  # don't need x in standata when outcome_type is binary (for now), and also
  # don't need alpha1, gamma1, sigma_beta1, etc in parlist
  if (outcome_type == "binary") {
    standata$x <- NULL
    betapars <- "beta0"
    kappapars <- "kappa0"
    plist <- c("alpha1", "gamma1", "sigma_beta1", "sigma_y")
    parlist <- parlist[!(parlist %in% plist)]
  }

  print(str(standata))
  rm(pop_data)
  rm(sample_data)

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

  #fit <- check_for_div_trans(num_clusters, num_units, use_sizes, outcome_type,
  #                           size_model, rootdir, simno, stanmod, standata
  #                           num_iter, num_chains, fit)

  if (simno == 1) {
    saveRDS(fit, file = paste0(rootdir, "/output/simulation/stanfit_usesizes_",
                               use_sizes, "_", outcome_type, "_", size_model,
                                "_", model_name, "_nclusters_", num_clusters,
                               "_nunits_", nunits, "_simno_", simno,".rds"))
    print("done saving stanfit object")
    cat("filename:", paste0(rootdir, "/output/simulation/stanfit_usesizes_",
                               use_sizes, "_", outcome_type, "_", size_model,
                                "_", model_name, "_nclusters_", num_clusters,
                               "_nunits_", nunits, "_simno_", simno,".rds"))
    print(Sys.time())
  }

  # Make sure things converged in terms of Rhat too
  rhat_check_df <- data.frame(summary(fit)$summary)
  rhat_check <- rhat_check_df$Rhat
  if (max(rhat_check) >= 1.1) {
    out_msg <- paste0("Rhat check failed! Max Rhat is", max(rhat_check))
    cat("Rhat check for alpha0 failed! Rhat is", max(rhat_check), "\n")
    print(rhat_check_df)
    print(summary(rhat_check_df[!grepl("\\[", rownames(rhat_check_df)), ]))
    saveRDS(out_msg,
            paste0(rootdir, "output/simulation/stan_results_usesizes_",
                   use_sizes, "_", outcome_type, "_", size_model, "_",
                   model_name, "_nclusters_", num_clusters,
                   "_nunits_", nunits, "_sim_", simno, ".rds"))
    saveRDS(fit, file = paste0(rootdir, "/output/simulation/stanfit_usesizes_",
                               use_sizes, "_", outcome_type, "_", size_model,
                                "_", model_name, "_nclusters_", num_clusters,
                               "_nunits_", nunits, "_simno_", simno,".rds"))
    return(NULL)
  }

  # Print summary for the basic parameters
  print(fit, pars = parlist)

  ##########################################
  ### Extract samples, format param estimates to just keep necessary info 
  ##########################################
  source(paste0(rootdir, "src/simulation/process_fit.R"))
  #kappapars <- ifelse(exists("kappapars"), kappapars, NA)
  #sbpars <- ifelse(exists("sbpars"), sbpars, NA)
  #res <- process_fit(fit, parlist, betapars, kappapars, sbpars, stanmod_name)

  # Pull out the objects from the list "res"
  #for (j in names(res)) {
  #  assign(j, res[[j]])
  #}

  ##########################################
  ### Draw new cluster sizes (if needed) and ybar_new
  ##########################################
  source(paste0(rootdir, "src/simulation/make_pp_draws.R"))

  print("done drawing Nj_new, ybar_new")
  print(Sys.time())      

  #yn_list <- list(ybar_new = ybar_new, Nj_new_df = Nj_new_df)
  #saveRDS(yn_list,
  #        paste0(rootdir, "output/simulation/yn_list_usesizes_",
  #               use_sizes, "_", outcome_type, "_", size_model, "_", model_name,
  #               "_nclusters_", num_clusters, "_nunits_", nunits,
  #               "_sim_", simno, ".rds"))

  ##########################################
  ### Plot draws of cluster sizes vs truth
  ##########################################
  if ((grepl("bb", stanmod_name) || grepl("lognormal", stanmod_name) ||
       grepl("negbin", stanmod_name)) && simno == 1) {
    Nj_new_df$in_sample <- Nj_new_df$cluster_id <= K
    tmpdf <- data.frame(cluster_id = c(1:J), draw_num = 9999,
                        Nj_new = Nj_pop, in_sample = FALSE)
    Nj_new_df <- rbind(Nj_new_df, tmpdf)
    Nj_new_df$draw_num <- as.integer(Nj_new_df$draw_num)
    Nj_new_df$is_truth <- ifelse(Nj_new_df$draw_num == 9999, "truth", "draws")
    Nj_sam_df <- data.frame(Nj_sample)
    Nj_pop_df <- data.frame(Nj_pop)
    plt <- ggplot(Nj_new_df, aes(x = Nj_new)) +
      geom_line(aes(group = draw_num, colour = is_truth), stat = "density") +
      geom_line(data = Nj_sam_df, stat = "density",
                aes(x = Nj_sample, colour = "sample")) +
      scale_colour_manual("", values = c("truth" = "black", "draws" = "grey50",
                                         "sample" = "grey80")) +
      scale_x_continuous(limits = c(0, max(Nj_new))) +
      xlab("Cluster size") +
      ylab("Density") +
      ggtitle(paste0("Model: ", stanmod_name, ", ", size_model)) +
      theme_bw()
    ggsave(plt, file = paste0(rootdir, "/output/figures/Nj_draws_usesizes_",
                              use_sizes, "_", outcome_type, "_", size_model, "_",
                              model_name, "_nclusters_", num_clusters,
                              "_nunits_", nunits, "_sim_", simno, ".png"),
           width = 10, height = 8)
  }
  print("done plotting Nj_new")
  print(Sys.time())      

  ##########################################
  ### Plot draws of ybar_new
  ##########################################
  if (simno == 1) {
    tmpdf <- data.frame(ybar_new, truth = ybar_true)
    plt <- ggplot(tmpdf, aes(x = ybar_new)) +
      geom_vline(aes(xintercept = truth), colour = "grey80") +
      geom_line(stat = "density") +
      xlab("Draws of ybar_new") +
      ylab("Density") +
      ggtitle(paste0("Model: ", stanmod_name, ", ", size_model)) +
      theme_bw()
    ggsave(plt, file = paste0(rootdir, "/output/figures/ybar_new_draws_usesizes_",
                              use_sizes, "_", outcome_type, "_", size_model, "_",
                              model_name, "_nclusters_", num_clusters,
                              "_nunits_", nunits, "_sim_", simno, ".png"),
           width = 10, height = 8)
  }
  print("done making plots")
  print(Sys.time())      

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
