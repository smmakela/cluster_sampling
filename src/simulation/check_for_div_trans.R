
#check_for_div_trans <- function(num_clusters, num_units, use_sizes, outcome_type,
#                                size_model, rootdir, simno, stanmod, standata
#                                num_iter, num_chains, fit) {
  # num_clusters -- number of clusters to sample
  # num_units -- number of units to sample
  # use_sizes -- 0/1 for whether cluster sizes used in pop data
  # outcome_type -- whether outcome is continuous or binary
  # size_model -- model used to create pop cluster sizes
  # rootdir -- root directory where src, output, etc folders are
  # simno -- current iteration; used so that multiple instances aren't trying to write to the same file
  # stanmod -- compiled stan model
  # standata -- data for stan model
  # num_iter -- number of iterations stan should run for
  # num_chains -- number of chains to run in stan
  # fit -- output from sampling()

  # Check if there are any divergent transitions in the current fit. if no, end.
  # Otherwise, run with a higher value of adapt_delta until we switch to the
  # ncp version of the model. If that fails, then return a message saying we
  # can't get rid of the divergent transitions

    samp_params <- get_sampler_params(fit)
    num_div_trans <- 0
    for (s in 1:length(samp_params)) {                                
      num_div_trans <- num_div_trans + sum(samp_params[[s]][((num_iter/2)+1):num_iter, "divergent__"])
    }
    ad_val <- 0.999
    if (num_div_trans > 0) {
      counter <- 1
      while (counter <= 3 & num_div_trans > 0) {
        ad_val <- ad_val + 9 * (10^(-1*(counter + 3)))
        cat("----------------------------------------------------------------\n")
        cat("Divergent transitions happened, rerunning with adapt_delta =", ad_val, "\n")
        cat("num_clusters =", num_clusters, "num_units=", num_units, "\n")
        fit <- sampling(stanmod, data = standata,
                        iter = num_iter, chains = num_chains,
                        control = list(stepsize = 0.001, adapt_delta = ad_val))
        print("done fitting stan model")
        print(Sys.time())
        print(warnings())
        samp_params <- get_sampler_params(fit)
        num_div_trans <- 0
        for (s in 1:length(samp_params)) {                                
          num_div_trans <- num_div_trans + sum(samp_params[[s]][((num_iter/2)+1):num_iter, "divergent__"])
        }
        counter <- counter + 1
      } # end while
      # if we still have divergent transitions, try the ncp version (only for continuous)
      if (num_div_trans > 0) {
        cat("----------------------------------------------------------------\n")
        cat("Divergent transitions remain, running NCP version\n")
        cat("num_clusters =", num_clusters, "num_units=", num_units, "\n")
        stanmod <- readRDS(paste0(rootdir, "/src/analysis/", stanmod_name, "_ncp.rds"))
        expose_stan_functions(stanmod)
        fit <- sampling(stanmod, data = standata,
                      iter = num_iter, chains = num_chains,
                      control = list(stepsize = 0.001, adapt_delta = 0.999))
        print("done fitting stan model")
        print(Sys.time())
        print(warnings())
        samp_params <- get_sampler_params(fit)
        num_div_trans <- 0
        for (s in 1:length(samp_params)) {                                
          num_div_trans <- num_div_trans + sum(samp_params[[s]][((num_iter/2)+1):num_iter, "divergent__"])
        }
        # try raising adapt_delta if necessary
        ad_val <- 0.999
        if (num_div_trans > 0) {
          counter <- 1
          while (counter <= 3 & num_div_trans > 0) {
            ad_val <- ad_val + 9 * (10^(-1*(counter + 3)))
            cat("----------------------------------------------------------------\n")
            cat("Divergent transitions happened for NCP, rerunning with adapt_delta =", ad_val, "\n")
            cat("num_clusters =", num_clusters, "num_units=", num_units, "\n")
            fit <- sampling(stanmod, data = standata,
                            iter = num_iter, chains = num_chains,
                            control = list(stepsize = 0.001, adapt_delta = ad_val))
            print("done fitting stan model")
            print(Sys.time())
            print(warnings())
            samp_params <- get_sampler_params(fit)
            num_div_trans <- 0
            for (s in 1:length(samp_params)) {                                
              num_div_trans <- num_div_trans + sum(samp_params[[s]][((num_iter/2)+1):num_iter, "divergent__"])
            }
            counter <- counter + 1
          } # end while
          # if we *still* have divergent transitions, give up
          if (num_div_trans > 0) {
            cat("Unable to get rid of divergent transitions :( \n")
            cat("num_clusters =", num_clusters, "num_units=", num_units, "\n")
            out_msg <- paste0("There were ", num_div_trans, " divergent transitions_")
            saveRDS(out_msg,
                    paste0(rootdir, "output/simulation/stan_results_usesizes_",
                           use_sizes, "_", outcome_type, "_", size_model, "_",
                           model_name, "_nclusters_", num_clusters,
                           "_nunits_", nunits, "_sim_", simno, ".rds"))
            #saveRDS(fit,
            #        paste0(rootdir, "/output/simulation/stanfit_usesizes_",
            #               use_sizes, "_", outcome_type, "_", size_model,
            #               "_", model_name, "_nclusters_", num_clusters,
            #               "_nunits_", nunits, "_simno_", simno,".rds"))
            return(NULL)
          } # end if for num_div_trans > 0 after ncp
        } # end if for raising adapt_delta in ncp
      } # end if for doing ncp
    } # end if for initial num_div_trans > 0

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

#  return(fit)

#}

