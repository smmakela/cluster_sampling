# Author: Susanna Makela
# Date: 21 Apr 2014
# Purpose: sample data from population

sampledata <- function(num_clusters, num_units, use_sizes, outcome_type,
                       size_model) {
  # num_clusters -- number of clusters to sample
  # num_units -- number of units to sample
  # use_sizes -- 0/1 for whether cluster sizes were used to generate pop data
  # outcome_type -- whether outcome is normal or binary
  # size_model -- model used to generate pop cluster sizes

  #############################################################################
  ### Useful functions -- sample nunits units in each clsuter
  #############################################################################
    unit_sampler <- function(cluster_size, nunits) {
      selection <- sort(sample.int(cluster_size, nunits, replace = FALSE))
      return(selection)
    }

  #############################################################################
  ### Load the pop data
  #############################################################################
    popdata <- readRDS(file = paste0(rootdir,
                                     "/output/simulation/popdata_usesizes_",
                                     use_sizes, "_", outcome_type, "_",
                                     size_model, ".rds"))

    # Unlist the items in popdata into their own objects
    nms <- names(popdata)
    for (j in nms) {
      assign(j, popdata[[j]])
    }

  #############################################################################
  ### Do the sampling -- FF case
  #############################################################################
  if (size_model == "ff") {
    # Sample clusters, stratified by stratum_id -- here we won't get
    # inclusion probabilities but it doesn't matter since we're not going to
    # use the Sarndal formulas in svy_ests.R for size_model = "ff", which is
    # the only place we'd need them for anyway
    strat_dat <- data.frame(stratum_id, Mj, cluster_id = c(1:J), tot_births)
    strat_dat <- dplyr::arrange(strat_dat, stratum_id)
    strat_dat$num_units_to_sample <- ifelse(strat_dat$stratum_id == 9, 100, 325)
    sampled_cluster_list_orig <- ppssstrat(strat_dat$Mj, strat_dat$stratum_id,
                                           c(rep(1, times = 8), 8))

    # ppsstrat returns cluster ids in the sorted order (from arrange()), so
    # the cluster id's it returns don't correspond to our actual cluster ids.
    # use the output of ppsstrat to index into the actual cluster id variable
    # to get the correct sampled cluster ids
    sampled_cluster_list <- strat_dat$cluster_id[sampled_cluster_list_orig]

    # Also pull out the correct number of units to be sampled and the total
    # number of units in each cluster
    num_units_vec <- strat_dat$num_units_to_sample[sampled_cluster_list_orig]
    tot_births <- strat_dat$tot_births[sampled_cluster_list_orig]

    # RENUMBER cluster ids so that sampled clusters run from 1:num_clusters
    #   and the rest from (num_clusters + 1):J
    all_cluster_ids <- sort(unique(pop_data$cluster_id))
    nonsampled_cluster_ids <- sort(setdiff(all_cluster_ids, sampled_cluster_list))
    idmap <- data_frame(cluster_id = c(sampled_cluster_list, nonsampled_cluster_ids), new_cluster_id = c(1:J))
    pop_data <- merge(pop_data, idmap, by = "cluster_id")
    pop_data$orig_cluster_id <- pop_data$cluster_id
    pop_data$cluster_id <- pop_data$new_cluster_id
    pop_data$new_cluster_id <- NULL
    # pull out Mj again
    tmp <- dplyr::distinct(pop_data, cluster_id, Mj)
    tmp <- dplyr::arrange(tmp, cluster_id)
    Mj <- tmp$Mj

    # use unit_sampler to sample the right number of units in each sampled
    # cluster because we 1) renumbered clusters in pop_data and 2) remade Mj
    # from the updated pop_data with the renumbered clusters, we can just
    # pull out the first num_clusters obs in both Mj and num_units_vec
    tt <- cbind(cluster_size = tot_births, nunits = num_units_vec)
    sampled_unit_list <- plyr::mlply(tt, unit_sampler)

    pop_data$insample <- 0
    # pull out sampled data
    for (j in 1:num_clusters) {
      curr_units <- sampled_unit_list[[j]]
      pop_data <- dplyr::mutate(pop_data,
                                insample = replace(insample,
                                                   cluster_id == j &
                                                    unit_id %in% curr_units, 1))
    }
    sample_data <- dplyr::filter(pop_data, insample == 1)
    sample_data$insample <- NULL

  #############################################################################
  ### Do the sampling -- non-FF case
  #############################################################################
  } else { # end if size_model == "ff"
    # Sample clusters
    pik <- inclusionprobabilities(Mj, num_clusters)
    PI_full <- UPtillepi2(pik)
    sampled_cluster_list <- c(1:J)[UPtille(pik) == 1]
    rm(Mj) # remove Mj b/c will overwrite it after clusters are renumbered

    # RENUMBER cluster ids so that sampled clusters run from 1:num_clusters
    #   and the rest from (num_clusters + 1):J
    all_cluster_ids <- sort(unique(pop_data$cluster_id))
    nonsampled_cluster_ids <- sort(setdiff(all_cluster_ids, sampled_cluster_list))
    idmap <- data_frame(cluster_id = c(sampled_cluster_list, nonsampled_cluster_ids), new_cluster_id = c(1:J))
    pop_data <- merge(pop_data, idmap, by = "cluster_id")
    pop_data$orig_cluster_id <- pop_data$cluster_id
    pop_data$cluster_id <- pop_data$new_cluster_id
    pop_data$new_cluster_id <- NULL
    # pull out Mj again
    tmp <- dplyr::distinct(pop_data, cluster_id, Mj)
    tmp <- dplyr::arrange(tmp, cluster_id)
    Mj <- tmp$Mj
    # if num_units <= 1, then it's actually a proportion, so we need to convert
    # it to a vector of integers by multiplying by Mj
    if (num_units <= 1) {
      num_units_vec <- round(Mj*num_units)
      # need at least 2 samples per cluster so that we can estimate a variance,
      # so if any elements of num_units_vec are < 2, set them to at least 2
      bad_inds <- which(num_units_vec < 2)
      num_units_vec[bad_inds] <- num_units_vec[bad_inds] + 1
    } else {
      num_units_vec <- rep(num_units, times = J)
    }

    # reorder PI_full so that the sampled clusters are in the first 1:num_clusters
    # rows/columns
    new_inds <- c(sampled_cluster_list, nonsampled_cluster_ids)
    PI_ij <- PI_full[new_inds, new_inds]
    PI_i <- pik[new_inds]

    # use unit_sampler to sample the right number of units in each sampled
    # cluster because we 1) renumbered clusters in pop_data and 2) remade Mj
    # from the updated pop_data with the renumbered clusters, we can just
    # pull out the first num_clusters obs in both Mj and num_units_vec
    tt <- cbind(Mj[1:num_clusters], num_units_vec[1:num_clusters])
    sampled_unit_list <- plyr::mlply(tt, unit_sampler)

    pop_data$insample <- 0
    # pull out sampled data
    for (j in 1:num_clusters) {
      curr_units <- sampled_unit_list[[j]]
      pop_data <- dplyr::mutate(pop_data,
                                insample = replace(insample,
                                                   cluster_id == j &
                                                    unit_id %in% curr_units, 1))
    }
    sample_data <- dplyr::filter(pop_data, insample == 1)
    sample_data$insample <- NULL
  }

  ##########################################
  ### Save
  ##########################################
    if (size_model == "ff") {
      simdata <- list(pop_data = pop_data, sample_data = sample_data, J = J,
                      Mj = Mj, logMj_c = logMj_c)
    } else {
      simdata <- list(pop_data = pop_data, sample_data = sample_data, J = J,
                      Mj = Mj, logMj_c = logMj_c, PI_i = PI_i, PI_ij = PI_ij)
    }

    return(simdata)
}
