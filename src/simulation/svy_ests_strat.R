# Author: Susanna Makela
# Date: Jan 28, 2016
# Purpose: use the survey package to estimate ybar from the sampled data
svy_ests_strat <- function(J, num_clusters, num_units, use_sizes, outcome_type,
                           size_model, sim_data, rootdir, simno) {
  # Inputs:
  #   J -- number of clusters in pop
  #   num_clusters -- number of sampled clusters
  #   num_units -- number of sampled units
  #   use_sizes -- whether cluster sizes are related to outcome values (0/1)
  #   outcome_type -- whether outcomes are continuous or binary
  #   size_model -- model used to generate pop cluster sizes
  #   sim_data -- list output from sampledata function, has pop and sampled data
  #   rootdir -- directory where everything is stored
  #   simno -- which simulation this is

  library(survey)

  #############################################################################
  # Load pop and simno data
  #############################################################################
  # load pop data -- use this to get the total sample size so we can calculate
  # the sampling probs (here we are assuming we know M_j in all clusters, so
  # obviously we can calculate the sampling probs
    if (num_units <= 1) {
      nunits <- paste(num_units*100, "pct", sep = "")
    } else {
      nunits <- num_units
    }
    # pull everything out of sim_data and assign to its name
    for (j in names(sim_data)) {
      assign(j, sim_data[[j]])
    }
    popdata <- readRDS(paste0(rootdir, "/output/simulation/popdata_usesizes_",
                              use_sizes, "_", outcome_type, "_", size_model, ".rds"))
    sizetot <- sum(popdata[["Mj"]])
    ybar_true <- mean(popdata[["pop_data"]]$y)
    truepars <- popdata[["truepars"]]
    print(truepars)
    rm(popdata)
    rm(sim_data)

    Nj_sample <- Mj[1:num_clusters]

    # fpc is the number of clusters in the pop
    sample_data$fpc <- J 
    # prob of selecting the cluster
    sample_data$prob_cluster <- num_clusters*sample_data$Mj/sizetot 
    # prob of selecting each unit within the cluster
    if (num_units > 1) {
      sample_data$prob_unit <- num_units/sample_data$Mj 
    } else {
      sample_data$prob_unit <- num_units
    }
    sample_data$wt <- 1/sample_data$prob_cluster

  #############################################################################
  # Estimate mean + std error
  #############################################################################
  # Set up survey design
    sample_data <- sample_data %>%
      dplyr::group_by(stratum_id) %>%
      dplyr::mutate(num_clusters_per_stratum = n_distinct(cluster_id))
    strat_design <- svydesign(id = ~cluster_id + unit_id, data = sample_data,
                              strata = ~stratum_id,
                              fpc = ~num_clusters_per_stratum + tot_births)

  # Hajek (????) estimates
    ests <- svymean(~y, design = strat_design, na.rm = TRUE)
    ybar_hat_hajek <- as.numeric(ests[1])
    ybar_se_hajek <- as.numeric(sqrt(attr(ests, "var")))

  # GREG estimates
    if (outcome_type == "continuous") {
      ptot <- sum(pop_data$x)
      pop_totals <- c(`(Intercept)`=nrow(sample_data), x = ptot)
      if (num_units > 1) { # self-weighting
        strat_design2 <- calibrate(strat_design, formula = ~x, pop_totals)
      } else { # make weights be the same by cluster
        strat_design2 <- calibrate(strat_design, formula = ~x, pop_totals, aggregate = 1) 
      }
      ests2 <- svymean(~y, design = strat_design2)
      ybar_hat_greg <- as.numeric(ests2[1])
      ybar_se_greg <- as.numeric(sqrt(attr(ests2, "var")))
    } else { 
      # for binary outcome, no greg estimator
      ybar_hat_greg <- NA
      ybar_se_greg <- NA
    }

  #############################################################################
  # Print and save results
  #############################################################################
    print("**************************************************")
    print(paste0("ybar true: ", round(ybar_true, digits = 2)))
    print(paste0("Hajek estimate (se): ", round(ybar_hat_hajek, digits = 2),
                 " (", round(ybar_se_hajek, digits = 2), ")"))
    print(paste0("GREG estimate (se): ", round(ybar_hat_greg, digits = 2),
                 " (", round(ybar_se_greg, digits = 2), ")"))
    print("**************************************************")

  # store results so that statistics are wide, model names are long
    res <- data.frame(ybar_true,
                      ybar_hat_hajek, ybar_se_hajek,
                      ybar_hat_greg, ybar_se_greg)
print(res)
    res %>%
      tidyr::gather(key = tmpname, value = value, -ybar_true) %>%
      tidyr::extract(col = tmpname, into = c("stat_name", "model_name"),
                    regex = "^([^_]*_[^_]*)_(.*)$") %>%
      tidyr::spread(stat_name, value) -> res
    res$model_name <- gsub("\\_", "", res$model_name)
    res$ybar_hat_lci50 <- res$ybar_hat - res$ybar_se
    res$ybar_hat_uci50 <- res$ybar_hat + res$ybar_se
    res$ybar_hat_lci95 <- res$ybar_hat - 1.96*res$ybar_se
    res$ybar_hat_uci95 <- res$ybar_hat + 1.96*res$ybar_se

    saveRDS(res,
            paste0(rootdir, "output/simulation/svy_ests_usesizes_",
                   use_sizes, "_", outcome_type, "_", size_model, "_nclusters_",
                   num_clusters, "_nunits_", nunits, "_simno_", simno, ".rds"))
    return(NULL)
}

