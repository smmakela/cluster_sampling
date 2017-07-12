#############################################################################
### Set lib paths, source files 
#############################################################################
  libdir <- "/vega/stats/users/smm2253/rpackages"
  .libPaths(libdir)
  rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
  Sys.setenv(HOME = rootdir)
  
#############################################################################
### Count the number of sim check files for the current option values
#############################################################################
  num_sims <- 100
  use_sizes <- c(0, 1)
  outcome_types <- c("binary", "continuous")
  size_models <- c("multinomial", "poisson")
  model_names <- c("bb", "cluster_inds_only", "knowsizes", "lognormal", "negbin")
  num_clusters_list <- c(5, 10, 20, 30)
  num_units_list <- c(0.05, 0.1, 0.25, 0.5, 1, 10, 30, 60)

  # Data frame to hold combos with missing sims
  miss_sims_stan <- data.frame()
  miss_sims_db <- data.frame()
  miss_sims_ff <- data.frame()

  # Results -- write out out usesizes and outcome_type bc the loops are slow
    for (mm in model_names) {
      # 0/binary
      fil_list_ff <- list.files(path = paste0(rootdir, "output/simulation"),
                                pattern = paste0("stan_results_usesizes_0_",
                                                 "binary_", mm, 
                                                 "_ff_.*sim.*.txt"))
      if (length(fil_list_ff) != num_sims) {
        tmp <- data.frame(use_sizes = 0, outcome_type = "binary",
                          model_name = mm, num_sims = length(fil_list_ff))
        miss_sims_ff <- rbind(miss_sims_ff, tmp)
      }
      # 1/binary
      fil_list_ff <- list.files(path = paste0(rootdir, "output/simulation"),
                                pattern = paste0("stan_results_usesizes_1_",
                                                 "binary_", mm, 
                                                 "_ff_.*sim.*.txt"))
      if (length(fil_list_ff) != num_sims) {
        tmp <- data.frame(use_sizes = 1, outcome_type = "binary",
                          model_name = mm, num_sims = length(fil_list_ff))
        miss_sims_ff <- rbind(miss_sims_ff, tmp)
      }
      # 0/continuous
      fil_list_ff <- list.files(path = paste0(rootdir, "output/simulation"),
                                pattern = paste0("stan_results_usesizes_0_",
                                                 "continuous_", mm, 
                                                 "_ff_.*sim.*.txt"))
      if (length(fil_list_ff) != num_sims) {
        tmp <- data.frame(use_sizes = 0, outcome_type = "continuous",
                          model_name = mm, num_sims = length(fil_list_ff))
        miss_sims_ff <- rbind(miss_sims_ff, tmp)
      }
      # 1/continuous
      fil_list_ff <- list.files(path = paste0(rootdir, "output/simulation"),
                                pattern = paste0("stan_results_usesizes_1_",
                                                 "continuous_", mm, 
                                                 "_ff_.*sim.*.txt"))
      if (length(fil_list_ff) != num_sims) {
        tmp <- data.frame(use_sizes = 1, outcome_type = "continuous",
                          model_name = mm, num_sims = length(fil_list_ff))
        miss_sims_ff <- rbind(miss_sims_ff, tmp)
      }

      for (ss in size_models) {
        for (nc in num_clusters_list) {
          for (nu in num_units_list) {
            # 0/binary
            fil_list_stan <- list.files(path = paste0(rootdir, "output/simulation"),
                                        pattern = paste0("stan_results_usesizes_0_",
                                                         "binary_", ss,
                                                         "_", mm, "_nclusters",
                                                         nc, "_nunits_", nu,
                                                         "_sim.*.txt"))
            if (length(fil_list_stan) != num_sims) {
              tmp <- data.frame(use_sizes = 0, outcome_type = "binary", 
                                size_model = ss, model_name = mm,
                                nclusters = nc, nunits = nu,
                                num_sims = length(fil_list_stan))
              miss_sims_stan <- rbind(miss_sims_stan, tmp)
            }
            # 1/binary
            fil_list_stan <- list.files(path = paste0(rootdir, "output/simulation"),
                                        pattern = paste0("stan_results_usesizes_1_",
                                                         "binary_", ss,
                                                         "_", mm, "_nclusters",
                                                         nc, "_nunits_", nu,
                                                         "_sim.*.txt"))
            if (length(fil_list_stan) != num_sims) {
              tmp <- data.frame(use_sizes = 1, outcome_type = "binary", 
                                size_model = ss, model_name = mm,
                                nclusters = nc, nunits = nu,
                                num_sims = length(fil_list_stan))
              miss_sims_stan <- rbind(miss_sims_stan, tmp)
            }
            # 0/continuous
            fil_list_stan <- list.files(path = paste0(rootdir, "output/simulation"),
                                        pattern = paste0("stan_results_usesizes_0_",
                                                         "continuous_", ss,
                                                         "_", mm, "_nclusters",
                                                         nc, "_nunits_", nu,
                                                         "_sim.*.txt"))
            if (length(fil_list_stan) != num_sims) {
              tmp <- data.frame(use_sizes = 0, outcome_type = "continuous", 
                                size_model = ss, model_name = mm,
                                nclusters = nc, nunits = nu,
                                num_sims = length(fil_list_stan))
              miss_sims_stan <- rbind(miss_sims_stan, tmp)
            }
            # 1/continuous
            fil_list_stan <- list.files(path = paste0(rootdir, "output/simulation"),
                                        pattern = paste0("stan_results_usesizes_1_",
                                                         "continuous_", ss,
                                                         "_", mm, "_nclusters",
                                                         nc, "_nunits_", nu,
                                                         "_sim.*.txt"))
            if (length(fil_list_stan) != num_sims) {
              tmp <- data.frame(use_sizes = 1, outcome_type = "continuous", 
                                size_model = ss, model_name = mm,
                                nclusters = nc, nunits = nu,
                                num_sims = length(fil_list_stan))
              miss_sims_stan <- rbind(miss_sims_stan, tmp)
            }

            # This will get repeated for every value of mm which is not
            # necessary but didn't feel like putting in an if-statement
            # 0/binary
            fil_list_db <- list.files(path = paste0(rootdir, "output/simulation"),
                                      pattern = paste0("stan_results_usesizes_0_",
                                                       "binary_", ss,
                                                       "_nclusters",
                                                       nc, "_nunits_", nu,
                                                       "_sim.*.txt"))
            if (length(fil_list_db) != num_sims) {
              tmp <- data.frame(use_sizes = 0, outcome_type = "binary", 
                                size_model = ss,
                                nclusters = nc, nunits = nu,
                                num_sims = length(fil_list_db))
              miss_sims_db <- rbind(miss_sims_db, tmp)
            }
            # 1/binary
            fil_list_db <- list.files(path = paste0(rootdir, "output/simulation"),
                                      pattern = paste0("stan_results_usesizes_1_",
                                                       "binary_", ss,
                                                       "_nclusters",
                                                       nc, "_nunits_", nu,
                                                       "_sim.*.txt"))
            if (length(fil_list_db) != num_sims) {
              tmp <- data.frame(use_sizes = 1, outcome_type = "binary", 
                                size_model = ss,
                                nclusters = nc, nunits = nu,
                                num_sims = length(fil_list_db))
              miss_sims_db <- rbind(miss_sims_db, tmp)
            }
            # 0/continuous
            fil_list_db <- list.files(path = paste0(rootdir, "output/simulation"),
                                      pattern = paste0("stan_results_usesizes_0_",
                                                       "continuous_", ss,
                                                       "_nclusters",
                                                       nc, "_nunits_", nu,
                                                       "_sim.*.txt"))
            if (length(fil_list_db) != num_sims) {
              tmp <- data.frame(use_sizes = 0, outcome_type = "continuous", 
                                size_model = ss,
                                nclusters = nc, nunits = nu,
                                num_sims = length(fil_list_db))
              miss_sims_db <- rbind(miss_sims_db, tmp)
            }
            # 1/continuous
            fil_list_db <- list.files(path = paste0(rootdir, "output/simulation"),
                                      pattern = paste0("stan_results_usesizes_1_",
                                                       "continuous_", ss,
                                                       "_nclusters",
                                                       nc, "_nunits_", nu,
                                                       "_sim.*.txt"))
            if (length(fil_list_db) != num_sims) {
              tmp <- data.frame(use_sizes = 1, outcome_type = "continuous", 
                                size_model = ss,
                                nclusters = nc, nunits = nu,
                                num_sims = length(fil_list_db))
              miss_sims_db <- rbind(miss_sims_db, tmp)
            }
          } # end nu
        } # end nc
      } # end ss
    } # end mm

  # Save
  miss_all <- list(miss_sims_ff = miss_sims_ff,
                   miss_sims_stan = miss_sims_stan,
                   miss_sims_db = miss_sims_db)
  saveRDS(miss_all, file = paste0(rootdir, "output/simulation/check_files_all.rds"))


