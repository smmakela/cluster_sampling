# Author: Susanna Makela
# Date: 21 Apr 2014
# Purpose: sample data from population

sampledata <- function(num.clusters, num.units, use.sizes,
                       outcome.type, rootdir, sim) {
  # num.clusters -- number of clusters to sample
  # num.units -- number of units to sample
  # use.sizes -- 0/1 for whether cluster sizes were used to generate pop data
  # outcome.type -- whether outcome is normal or binary
  # rootdir -- root directory where Code, Data folders are
  # sim -- current simulation; used so that multiple instances aren't trying to write to the same file

  #############################################################################
  ### Load libraries and functions
  #############################################################################
    library(dplyr)
    source(paste(rootdir, "/src/simulation/rspps.r", sep = ""))

  #############################################################################
  ### Load the data we made
  #############################################################################
    load(file = paste0(rootdir, "//popdata_usesizes_", use.sizes, "_",
                       outcome.type, ".RData"))
    nms <- names(popdata)
    for (j in nms) {
      assign(j, popdata[[j]])
    }
    tmp <- dplyr::distinct(pop.data, cluster.id, Mj)
    Mj <- tmp$Mj

  #############################################################################
  ### Do the sampling
  #############################################################################
    # sample clusters
      sampled.cluster.list <- sort(rspps(Mj, c(1:J), num.clusters)) # randomized systematic PPS sampling
      #print("-----------`---------------------------------------")
      #print("sampled clusters:")   
      #print(sampled.cluster.list) 
      #print("---------------------------------------------------")
      rm(Mj) # remove this b/c will overwrite it once clusters are sampled and renumbered

    # RENUMBER cluster ids so that sampled clusters run from 1:num.clusters
    #   and the rest from (num.clusters + 1):J
      all.cluster.ids <- sort(unique(pop.data$cluster.id))
      nonsampled.cluster.ids <- sort(setdiff(all.cluster.ids, sampled.cluster.list))
      idmap <- data.frame(cluster.id = c(sampled.cluster.list, nonsampled.cluster.ids), new.cluster.id = c(1:J))
      pop.data <- merge(pop.data, idmap, by = "cluster.id")
      pop.data$orig.cluster.id <- pop.data$cluster.id
      pop.data$cluster.id <- pop.data$new.cluster.id
      pop.data$new.cluster.id <- NULL

    # pull out Mj again
      pop.data <- tbl_df(pop.data)
      tmp <- distinct(pop.data, cluster.id, Mj)
      tmp <- arrange(tmp, cluster.id)
      Mj <- tmp$Mj

    # if num.units<=1, then it's actually a proportion, so we need to convert it to a vector of integers
      if (num.units <= 1) {
        num.units.vec <- round(Mj*num.units)
        bad.inds <- which(num.units.vec < 2) # need at least 2 samples per cluster so that we can estimate a variance
        num.units.vec[bad.inds] <- num.units.vec[bad.inds] + 1
      } else {
        num.units.vec <- rep(num.units, times = J)
      }

    # function to sample num.units units in each sampled cluster
      unit.sampler <- function(unit.size, num.units) {
        selection <- sort(sample.int(unit.size, num.units, replace = FALSE))
        return(selection)
      }

    # use unit.sampler to sample the right number of units in each sampled cluster
    # because we 1) renumbered clusters in pop.data and 2) remade Mj from
    # the updated pop.data with the renumbered clusters, we can just
    # pull out the first num.clusters obs in both Mj and num.units.vec
      tt <- cbind(Mj[1:num.clusters], num.units.vec[1:num.clusters])
      sampled.unit.list <- mlply(tt, unit.sampler)

    # pull out sampled data
      pop.data$insample <- 0
      for (j in 1:num.clusters) {
        curr.units <- sampled.unit.list[[j]]
        pop.data$insample[pop.data$cluster.id == j & pop.data$unit.id %in% curr.units] <- 1
      }

      #pop.data$new.cluster.id <- NA
      #for (j in 1:length(sampled.cluster.list)) {
      #  curr.cluster <- sampled.cluster.list[j]
      #  curr.units <- sampled.unit.df[[j]]
      #  pop.data$insample[pop.data$cluster.id == curr.cluster & pop.data$unit.id %in% curr.units] <- 1
      #  #pop.data$new.cluster.id[pop.data$cluster.id == curr.cluster] <- j
      #}
      sample.data <- pop.data[pop.data$insample == 1, ]
      sample.data$insample <- NULL

  ##########################################
  ### Save
  ##########################################
    if (num.units <= 1) {
      nunits <- paste(100*num.units, "pct", sep = "")
    } else {
      nunits <- num.units
    }
    save(sample.data, file = paste(rootdir, "/Data/Simplify/vary_K/sampledata_usesizes_", use.sizes, 
                                   "_nclusters_", num.clusters, "_nunits_", nunits, "_sim_", sim, ".RData", sep = ""))

    # also save pop.data with insample indicator
    pop.data.inds <- pop.data
    save(pop.data.inds, file = paste(rootdir, "/Data/Simplify/vary_K/popdata_usesizes_", use.sizes,
                                     "_inds_nclusters_", num.clusters, "_nunits_", nunits,
                                     "_sim_", sim, ".RData", sep = ""))
                                     
  
}
