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
print("str(popdata, give.attr = FALSE)")
print(str(popdata, give.attr = FALSE))
    # Unlist the items in popdata into their own objects
    nms <- names(popdata)
    for (j in nms) {
      assign(j, popdata[[j]])
    }

  #############################################################################
  ### Do the sampling -- FF case
  #############################################################################
  if (grepl("ff", size_model)) {
    pop_cluster_data <- pop_cluster_data %>%
      dplyr::mutate(cluster_id = c(1:J))
    #pop_cluster_data$num_units_to_sample <- ifelse(pop_cluster_data$stratum_id == 9, 100, 325)
    pop_cluster_data$num_units_to_sample <- ifelse(pop_cluster_data$stratum_id == 2, 100, 325)
    # stratified sampling
    if (size_model == "ffstrat") {
      pop_cluster_data_withsamp <- data.frame()
      for (s in 1:2) {
        curr_pop_cluster_data <- pop_cluster_data[pop_cluster_data$stratum_id == s, ]
        curr_pop_cluster_data$pik <- inclusionprobabilities(curr_pop_cluster_data$Mj, num_clusters/2)
        curr_pop_cluster_data$is_sampled_cluster <- UPtille(curr_pop_cluster_data$pik) == 1
        curr_pop_cluster_data$is_certainty <- curr_pop_cluster_data$pik >= 1
        pop_cluster_data_withsamp <- rbind(pop_cluster_data_withsamp, curr_pop_cluster_data)
      }
      pop_cluster_data <- pop_cluster_data_withsamp
    } else { # non-stratified
      pop_cluster_data$pik <- inclusionprobabilities(pop_cluster_data$Mj, num_clusters)
      pop_cluster_data$is_sampled_cluster <- UPtille(pop_cluster_data$pik) == 1
    }
print("str(pop_cluster_data, give.attr = FALSE):")
print(str(pop_cluster_data, give.attr = FALSE))
    # Create map for renumber clusters so that sampled clusters run from 
    # 1:num_clusters and the rest from (num_clusters + 1):J
    sampled_cluster_list <- pop_cluster_data$cluster_id[pop_cluster_data$is_sampled_cluster]
print("sampled_cluster_list")
print(sampled_cluster_list)
    all_cluster_ids <- c(1:J)
    nonsampled_cluster_ids <- sort(setdiff(all_cluster_ids, sampled_cluster_list))
    idmap <- data.frame(cluster_id = c(sampled_cluster_list,
                                       nonsampled_cluster_ids),
                        new_cluster_id = c(1:J))

print("idmap")
print(idmap)
    # Join map wither pop_cluster_data
    pop_cluster_data <- pop_cluster_data %>%
      dplyr::left_join(., idmap, by = "cluster_id") %>%
      dplyr::mutate(orig_cluster_id = cluster_id,
                    cluster_id = new_cluster_id) %>%
      dplyr::select(-new_cluster_id) %>%
      dplyr::arrange(cluster_id)
print("str(pop_cluster_data, give.attr = FALSE):")
print(str(pop_cluster_data, give.attr = FALSE))

    # Create data frame of just sampled clusters
    sampled_cluster_data <- pop_cluster_data %>%
      dplyr::filter(is_sampled_cluster == 1)
print("str(sampled_cluster_data, give.attr = FALSE):")
print(str(sampled_cluster_data, give.attr = FALSE))

    # Merge map to pop_data also
    pop_data <- pop_data %>%
     dplyr::rename(orig_cluster_id = cluster_id) %>%
     dplyr::left_join(., pop_cluster_data, by = "orig_cluster_id") %>%
     dplyr::select(-contains(".y")) %>%
     dplyr::rename_(.dots = setNames(names(.), gsub("\\.x", "", names(.))))
print("str(pop_data, give.attr = FALSE):")
print(str(pop_data, give.attr = FALSE))

    # Calculate PI_i and PI_ij
    PI_i <- pop_cluster_data$pik
    if (size_model == "ff") {
      PI_ij <- UPtillepi2(pop_cluster_data$pik) # pairwise selection probabilities
    }

    # UNITS ALREADY SAMPLED (see makepopdata.R), so we only need to keep the
    # clusters we just picked
    pop_data$insample <- 0
    pop_data$insample[pop_data$cluster_id %in% c(1:num_clusters)] <- 1
    sample_data <- pop_data %>%
      dplyr::filter(insample == 1) %>%
      dplyr::select(-insample)
print("str(sample_data, give.attr = FALSE):")
print(str(sample_data, give.attr = FALSE))

  #############################################################################
  ### Do the sampling -- non-FF case
  #############################################################################
  } else { # end if size_model == "ff"
    # Sample clusters
    pop_cluster_data$pik <- inclusionprobabilities(pop_cluster_data$Mj, num_clusters)
    PI_full <- UPtillepi2(pop_cluster_data$pik)
    pop_cluster_data$is_sampled_cluster <- UPtille(pop_cluster_data$pik) == 1

    # if num_units <= 1, then it's actually a proportion, so we need to convert
    # it to a vector of integers by multiplying by Mj
    if (num_units <= 1) {
      pop_cluster_data$num_units_to_sample <- round(pop_cluster_data$Mj * num_units)
      # need at least 2 samples per cluster so that we can estimate a variance,
      # so if any elements of num_units_to_sample are < 2, set them to at least 2
      bad_inds <- which(pop_cluster_data$num_units_to_sample < 2)
      pop_cluster_data$num_units_to_sample[bad_inds] <- pop_cluster_data$num_units_to_sample[bad_inds] + 1
    } else {
      pop_cluster_data$num_units_to_sample <- num_units
    }

    # Create map for renumber clusters so that sampled clusters run from 
    # 1:num_clusters and the rest from (num_clusters + 1):J
    sampled_cluster_list <- pop_cluster_data$cluster_id[pop_cluster_data$is_sampled_cluster]
    all_cluster_ids <- c(1:J)
    nonsampled_cluster_ids <- sort(setdiff(all_cluster_ids, sampled_cluster_list))
    idmap <- data.frame(cluster_id = c(sampled_cluster_list,
                                       nonsampled_cluster_ids),
                        new_cluster_id = c(1:J))

    # Join map wither pop_cluster_data
    pop_cluster_data <- pop_cluster_data %>%
      dplyr::left_join(., idmap, by = "cluster_id") %>%
      dplyr::mutate(orig_cluster_id = cluster_id,
                    cluster_id = new_cluster_id) %>%
      dplyr::select(-new_cluster_id) %>%
      dplyr::arrange(cluster_id)

    # Create data frame of just sampled clusters
    sampled_cluster_data <- pop_cluster_data %>%
      dplyr::filter(is_sampled_cluster == 1)
print("str(sampled_cluster_data, give.attr = FALSE)")
print(str(sampled_cluster_data, give.attr = FALSE))

    # Merge map to pop_data also
    pop_data <- pop_data %>%
     dplyr::rename(orig_cluster_id = cluster_id) %>%
     dplyr::left_join(., pop_cluster_data, by = "orig_cluster_id") %>%
     dplyr::select(-contains(".y")) %>%
     dplyr::rename_(.dots = setNames(names(.), gsub("\\.x", "", names(.))))
print("str(pop_data, give.attr = FALSE)")
print(str(pop_data, give.attr = FALSE))

    # Calculate PI_i and PI_ij
    PI_i <- pop_cluster_data$pik
    PI_ij <- UPtillepi2(pop_cluster_data$pik) # pairwise selection probabilities

    # Sample units within selected clusters by first inner_join-ing with
    # sampled_cluster_data, then by randomly sampling the units in each cluster
    # using a data.table trick (makes it go way faster)
    sample_data <- pop_data %>%
      dplyr::filter(is_sampled_cluster)
print("str(sample_data, give.attr = FALSE)")
print(str(sample_data, give.attr = FALSE))

    sample_data <- setDT(sample_data)
    sample_data <- sample_data[, .SD[sample(.N, num_units_to_sample)], by = cluster_id]
  }

  ##########################################
  ### Save
  ##########################################
    if (size_model == "ffstrat") {
      simdata <- list(pop_data = pop_data, sample_data = sample_data,
                      J = as.numeric(J),
                      pop_cluster_data = pop_cluster_data,
                      sampled_cluster_data = sampled_cluster_data,
                      truepars = truepars)
    } else {
      simdata <- list(pop_data = pop_data, sample_data = sample_data,
                      J = as.numeric(J),
                      pop_cluster_data = pop_cluster_data,
                      sampled_cluster_data = sampled_cluster_data,
                      PI_i = PI_i, PI_ij = PI_ij,
                      truepars = truepars)
    }

    return(simdata)
}
