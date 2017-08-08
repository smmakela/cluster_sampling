# Author: Susanna Makela
# Date: 21 Apr 2014
# Purpose: sample data from population

sampledata <- function(J, num_clusters, num_units, use_sizes, outcome_type,
                       size_model) {
  # J -- number of clusters in population
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
    if (grepl("ff", size_model)) {
      popdata <- readRDS(file = paste0(rootdir,
                                       "/output/simulation/popdata_usesizes_",
                                       use_sizes, "_", outcome_type, "_ff.rds"))
    } else {
      popdata <- readRDS(file = paste0(rootdir,
                                       "/output/simulation/popdata_usesizes_",
                                       use_sizes, "_", outcome_type, "_",
                                       size_model, ".rds"))
    }
print("str(popdata)")
print(str(popdata))
    # Unlist the items in popdata into their own objects
    nms <- names(popdata)
    for (j in nms) {
      assign(j, popdata[[j]])
    }

  #############################################################################
  ### Do the sampling -- FF case
  #############################################################################
  if (grepl("ff", size_model)) {
    cluster_level_data <- cluster_level_data %>%
      dplyr::mutate(cluster_id = c(1:J))
    #cluster_level_data$num_units_to_sample <- ifelse(cluster_level_data$stratum_id == 9, 100, 325)
    cluster_level_data$num_units_to_sample <- ifelse(cluster_level_data$stratum_id == 2, 100, 325)
    # stratified sampling
    if (size_model == "ffstrat") {
      cluster_level_data_withsamp <- data.frame()
      for (s in 1:2) {
        curr_cluster_level_data <- cluster_level_data[cluster_level_data$stratum_id == s, ]
        curr_cluster_level_data$pik <- inclusionprobabilities(curr_cluster_level_data$Mj, num_clusters/2)
        curr_cluster_level_data$is_sampled_cluster <- UPtille(curr_cluster_level_data$pik) == 1
        curr_cluster_level_data$is_certainty <- curr_cluster_level_data$pik >= 1
        cluster_level_data_withsamp <- rbind(cluster_level_data_withsamp, curr_cluster_level_data)
      }
      cluster_level_data <- cluster_level_data_withsamp
    } else { # non-stratified
      cluster_level_data$pik <- inclusionprobabilities(cluster_level_data$Mj, num_clusters)
      cluster_level_data$is_sampled_cluster <- UPtille(cluster_level_data$pik) == 1
    }
print("str(cluster_level_data):")
print(str(cluster_level_data))
    # Create map for renumber clusters so that sampled clusters run from 
    # 1:num_clusters and the rest from (num_clusters + 1):J
    sampled_cluster_list <- cluster_level_data$cluster_id[cluster_level_data$is_sampled_cluster]
print("sampled_cluster_list")
print(sampled_cluster_list)
    all_cluster_ids <- c(1:J)
    nonsampled_cluster_ids <- sort(setdiff(all_cluster_ids, sampled_cluster_list))
    idmap <- data.frame(cluster_id = c(sampled_cluster_list,
                                       nonsampled_cluster_ids),
                        new_cluster_id = c(1:J))

print("idmap")
print(idmap)
    # Join map wither cluster_level_data
    cluster_level_data <- cluster_level_data %>%
      dplyr::left_join(., idmap, by = "cluster_id") %>%
      dplyr::mutate(orig_cluster_id = cluster_id,
                    cluster_id = new_cluster_id) %>%
      dplyr::select(-new_cluster_id) %>%
      dplyr::arrange(cluster_id)
print("str(cluster_level_data):")
print(str(cluster_level_data))

    # Create data frame of just sampled clusters
    sampled_cluster_data <- cluster_level_data %>%
      dplyr::filter(is_sampled_cluster == 1)
print("str(sampled_cluster_data):")
print(str(sampled_cluster_data))

    # Merge map to pop_data also
    pop_data <- pop_data %>%
     dplyr::rename(orig_cluster_id = cluster_id) %>%
     dplyr::left_join(., cluster_level_data, by = "orig_cluster_id") %>%
     dplyr::select(-contains(".y")) %>%
     dplyr::rename_(.dots = setNames(names(.), gsub("\\.x", "", names(.))))
print("str(pop_data):")
print(str(pop_data))

    # Calculate PI_i and PI_ij
    PI_i <- cluster_level_data$pik
    if (size_model == "ff") {
      PI_ij <- UPtillepi2(cluster_level_data$pik) # pairwise selection probabilities
    }

    # use unit_sampler to sample the right number of units in each sampled
    # cluster because we 1) renumbered clusters in pop_data and 2) remade Mj
    # from the updated pop_data with the renumbered clusters, we can just
    # pull out the first num_clusters obs in both Mj and num_units_vec
    tt <- cbind(cluster_size = sampled_cluster_data$tot_births,
                nunits = sampled_cluster_data$num_units_to_sample)
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
    sample_data <- pop_data %>%
      dplyr::filter(insample == 1) %>%
      dplyr::select(-insample)
print("str(sample_data):")
print(str(sample_data))

  #############################################################################
  ### Do the sampling -- non-FF case
  #############################################################################
  } else { # end if size_model == "ff"
    # Sample clusters
    cluster_level_data$pik <- inclusionprobabilities(cluster_level_data$Mj, num_clusters)
    PI_full <- UPtillepi2(cluster_level_data$pik)
    cluster_level_data$is_sampled_cluster <- UPtille(cluster_level_data$pik) == 1

    # Create map for renumber clusters so that sampled clusters run from 
    # 1:num_clusters and the rest from (num_clusters + 1):J
    sampled_cluster_list <- cluster_level_data$cluster_id[cluster_level_data$is_sampled_cluster]
    all_cluster_ids <- c(1:J)
    nonsampled_cluster_ids <- sort(setdiff(all_cluster_ids, sampled_cluster_list))
    idmap <- data.frame(cluster_id = c(sampled_cluster_list,
                                       nonsampled_cluster_ids),
                        new_cluster_id = c(1:J))

    # Join map wither cluster_level_data
    cluster_level_data <- cluster_level_data %>%
      dplyr::left_join(., idmap, by = "cluster_id") %>%
      dplyr::mutate(orig_cluster_id = cluster_id,
                    cluster_id = new_cluster_id) %>%
      dplyr::select(-new_cluster_id) %>%
      dplyr::arrange(cluster_id)

    # Create data frame of just sampled clusters
    sampled_cluster_data <- cluster_level_data %>%
      dplyr::filter(is_sampled_cluster == 1)

    # Merge map to pop_data also
    pop_data <- pop_data %>%
     dplyr::rename(orig_cluster_id = cluster_id) %>%
     dplyr::left_join(., cluster_level_data, by = "orig_cluster_id") %>%
     dplyr::select(-contains(".y")) %>%
     dplyr::rename_(.dots = setNames(names(.), gsub("\\.x", "", names(.))))

    # Calculate PI_i and PI_ij
    PI_i <- cluster_level_data$pik
    PI_ij <- UPtillepi2(cluster_level_data$pik) # pairwise selection probabilities

    # if num_units <= 1, then it's actually a proportion, so we need to convert
    # it to a vector of integers by multiplying by Mj
    if (num_units <= 1) {
      sampled_cluster_data$num_units_vec <- round(sampled_cluster_data$Mj * num_units)
      # need at least 2 samples per cluster so that we can estimate a variance,
      # so if any elements of num_units_vec are < 2, set them to at least 2
      bad_inds <- which(sampled_cluster_data$num_units_vec < 2)
      sampled_cluster_data$num_units_vec[bad_inds] <- sampled_cluster_data$num_units_vec[bad_inds] + 1
    } else {
      sampled_cluster_data$num_units_vec <- num_units
    }

    # use unit_sampler to sample the right number of units in each sampled
    # cluster because we 1) renumbered clusters in pop_data and 2) remade Mj
    # from the updated pop_data with the renumbered clusters, we can just
    # pull out the first num_clusters obs in both Mj and num_units_vec
    tt <- cbind(cluster_size = sampled_cluster_data$Mj,
                nunits = sampled_cluster_data$num_units_vec)
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
    if (size_model == "ffstrat") {
      simdata <- list(pop_data = pop_data, sample_data = sample_data,
                      J = as.numeric(J),
                      cluster_level_data = cluster_level_data,
                      sampled_cluster_data = sampled_cluster_data)
    } else {
      simdata <- list(pop_data = pop_data, sample_data = sample_data,
                      J = as.numeric(J),
                      cluster_level_data = cluster_level_data,
                      sampled_cluster_data = sampled_cluster_data,
                      PI_i = PI_i, PI_ij = PI_ij)
    }

    return(simdata)
}
