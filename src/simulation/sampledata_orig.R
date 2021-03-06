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
  ### Source functions
  #############################################################################
    source(paste(rootdir, "/src/simulation/rspps.r", sep = ""))

  #############################################################################
  ### Load the data we made
  #############################################################################
    popdata <- readRDS(file = paste0(rootdir,
                                     "/output/simulation/popdata_usesizes_",
                                     use.sizes, "_", outcome.type, ".rds"))

    # Unlist the items in popdata into their own objects
    nms <- names(popdata)
    for (j in nms) {
      assign(j, popdata[[j]])
    }
if (num.units == 999) {
  sample.data <- pop.data
  pop.data$insample <- 1
} else {
  #############################################################################
  ### Do the sampling
  #############################################################################
    # sample clusters
    #sampled.cluster.list <- sort(rspps(Mj, c(1:J), num.clusters)) # randomized systematic PPS sampling
    #sampled.cluster.list <- sort(sample(c(1:J), num.clusters, prob = Mj, replace = FALSE))
    pik <- inclusionprobabilities(Mj, num.clusters)
    PI.full <- UPtillepi2(pik)
    sampled.cluster.list <- c(1:J)[UPtille(pik) == 1]
    rm(Mj) # remove Mj b/c will overwrite it after clusters are renumbered

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
    tmp <- dplyr::distinct(pop.data, cluster.id, Mj)
    tmp <- dplyr::arrange(tmp, cluster.id)
    Mj <- tmp$Mj
    # if num.units <= 1, then it's actually a proportion, so we need to convert
    # it to a vector of integers by multiplying by Mj
    if (num.units <= 1) {
      num.units.vec <- round(Mj*num.units)
      # need at least 2 samples per cluster so that we can estimate a variance,
      # so if any elements of num.units.vec are < 2, set them to at least 2
      bad.inds <- which(num.units.vec < 2)
      num.units.vec[bad.inds] <- num.units.vec[bad.inds] + 1
    } else {
      num.units.vec <- rep(num.units, times = J)
    }

    # reorder PI.full so that the sampled clusters are in the first 1:num.clusters
    # rows/columns
    new.inds <- c(sampled.cluster.list, nonsampled.cluster.ids)
    PI_ij <- PI.full[new.inds, new.inds]
    PI_i <- pik[new.inds]

    # function to sample num.units units in each sampled cluster
    unit.sampler <- function(cluster.size, num.units) {
      selection <- sort(sample.int(cluster.size, num.units, replace = FALSE))
      return(selection)
    }

    # use unit.sampler to sample the right number of units in each sampled
    # cluster because we 1) renumbered clusters in pop.data and 2) remade Mj
    # from the updated pop.data with the renumbered clusters, we can just
    # pull out the first num.clusters obs in both Mj and num.units.vec
    tt <- cbind(Mj[1:num.clusters], num.units.vec[1:num.clusters])
    sampled.unit.list <- plyr::mlply(tt, unit.sampler)

    pop.data$insample <- 0
    # pull out sampled data
    for (j in 1:num.clusters) {
      curr.units <- sampled.unit.list[[j]]
      pop.data <- dplyr::mutate(pop.data,
                                insample = replace(insample,
                                                   cluster.id == j &
                                                    unit.id %in% curr.units, 1))
    }
    sample.data <- dplyr::filter(pop.data, insample == 1)
    sample.data$insample <- NULL
}
  ##########################################
  ### Save
  ##########################################
    if (num.units <= 1) {
      nunits <- paste(100*num.units, "pct", sep = "")
    } else {
      nunits <- num.units
    }
    simdata <- list(pop.data = pop.data, sample.data = sample.data, J = J,
                    Mj = Mj, logMj_c = logMj_c, PI_i = PI_i, PI_ij = PI_ij)
    saveRDS(simdata,
            file = paste0(rootdir, "output/simulation/simdata_usesizes_",
                          use.sizes, "_", outcome.type,
                          "_nclusters_", num.clusters, "_nunits_", nunits,
                          "_simno_", sim, ".rds"))
}
