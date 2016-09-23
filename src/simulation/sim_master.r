# Author: Susanna Makela
# Date: 13 Nov 2014
# Purpose: master file for binomial model

sim_master <- function(sim, seed = NULL, numclusters, clustersize.range, unitcovar.range) {

  ##########################################
  ### Setup of directories
  ##########################################
    which.comp <- 1 # set to 0 for laptop, 1 for cluster
    if (which.comp == 1) {
      libdir <- "/vega/stats/users/smm2253/rpackages"
      .libPaths(libdir)
      rootdir <- "/vega/stats/users/smm2253/Projects/Cluster_Sampling/"
      Sys.setenv(HOME = rootdir)
    } else {
      rootdir <- "/Users/susanna/Google Drive/Survey Sampling/FF Modeling/"
    }

  ##########################################
  ### Source files
  ##########################################
    source(paste(rootdir, "/Code/Simplify/vary_K/sampledata.r", sep = ""))
    source(paste(rootdir, "/Code/Simplify/vary_K/runstan.r", sep = ""))
    source(paste(rootdir, "/Code/Simplify/vary_K/lmer_compare.r", sep = ""))
    source(paste(rootdir, "/Code/Simplify/vary_K/svy_ests.r", sep = ""))

  ##########################################
  ### Load libraries, set up parallel stuff
  ##########################################
    library(lme4)
    library(rstan)
    #library(foreach)
    #library(doParallel) 
    #library(iterators)
    #cl <- makeCluster(min(detectCores(), 10), outfile = paste(rootdir, "/Code/Simplify/outfiles/parallel_output.txt", sep = ""))
    #registerDoParallel(cl)
    #clusterEvalQ(cl, .libPaths( "/vega/stats/users/smm2253/rpackages"))
    #clusterEvalQ(cl, library(iterators, lib.loc = "/vega/stats/users/smm2253/rpackages"))
    #clusterEvalQ(cl, library(foreach, lib.loc = "/vega/stats/users/smm2253/rpackages"))
    #clusterEvalQ(cl, library(doParallel, lib.loc = "/vega/stats/users/smm2253/rpackages"))

  ##########################################
  ### Compile stan models
  ##########################################
    #print("About to compile stan models")
    #print(Sys.time())
    #if (stanmod_name == "bb") {
    #  stanmod <- stan_model(file = paste(rootdir, "/Code/Simplify/vary_K/bb.stan", sep = ""))
    #} else if (stanmod_name == "knowsizes") {
    #  stanmod <- stan_model(file = paste(rootdir, "/Code/Simplify/vary_K/knowsizes.stan", sep = ""))
    #} else if (stanmod_name == "negbin") {
    #  stanmod <- stan_model(file = paste(rootdir, "/Code/Simplify/vary_K/negbin.stan", sep = ""))
    #} else {
    #  stanmod <- stan_model(file = paste(rootdir, "/Code/Simplify/vary_K/cluster_inds_only.stan", sep = ""))
    #}
    
  ##########################################
  ### Loop through cluster/unit lists
  ##########################################
    stanmod_list <- c("cluster_inds_only", "knowsizes", "bb", "negbin")
    #stanmod_list <- c("negbin")
    popseed <- 1
    #num.cluster.list <- c(15, 30, 60, 100)
    num.cluster.list <- c(15, 30, 60)
    #num.cluster.list <- c(15)
    #loopres <- foreach (k = 1:length(num.cluster.list), .packages=c("rstan", "plyr"),
    #                    .export=ls(envir=globalenv()), .verbose=TRUE) %dopar% {
    for (u in c(0,1)) {
      use.sizes <- u
      for (k in 1:length(num.cluster.list)) {    
        num.clusters <- num.cluster.list[k]
        #num.units.list <- c(.05, .1, .25, .5, 1, 10, 50, 100)
        num.units.list <- c(.05, .1, 10, 50)
        #num.units.list <- c(.1)
        for (j in 1:length(num.units.list)) {
          num.units <- num.units.list[j]
          # skip the ones we already did
          if (num.clusters == 15 & num.units == 0.1) {
            next
          }
          # sample data
          print("Sampling data")
          print(Sys.time())
          sampledata(num.clusters, num.units, rootdir, sim, popseed, use.sizes)
          for (p in 1:length(stanmod_list)) {
            stanmod_name <- stanmod_list[p]
            ##### TEMP -- don't redo stuff we don't need to
            if (num.units <= 1) {
              nunits <- paste(100*num.units, "pct", sep = "")
            } else {
              nunits <- num.units
            }
            outname <- paste(rootdir, "/Results/Simplify/vary_K/ybar_stan_usesizes_", use.sizes, "_nclusters_", num.clusters,
                             "_nunits_", nunits, "_sim_", sim, "_", stanmod_name, ".txt", sep = "")
            tt <- file.info(outname)
            fdate <- tt$mtime
            fdate <- substr(fdate, 1, 10)
            if (file.exists(outname) & (fdate == "2016-07-29" | fdate == "2016-07-30")) {
              next
            }
            ##### TEMP
            cat("##################################################################################\n")
            cat("##################################################################################\n")
            print(paste("Starting to run use.sizes = ", use.sizes, ", stanmod = ", stanmod_name,
                        ", num.clusters = ", num.clusters, ", num.units = ", num.units, sep = ""))
            print(Sys.time())
            # run stan
            print("Compiling stan model")
            print(Sys.time())
            stanmod <- stan_model(file = paste(rootdir, "/Code/Simplify/vary_K/", stanmod_name, ".stan", sep = ""))
            print("Running stan")
            print(Sys.time())
            stanout <- runstan(num.clusters, num.units, rootdir, sim, popseed, stanmod, stanmod_name, use.sizes)
            cat("##################################################################################\n")
          } # end stanmod loop
          # compare using lmer
          print("Running lmer_compare")
          print(Sys.time())
          lmer_compare(num.clusters, num.units, rootdir, sim, popseed, use.sizes)
          # estimate ybar using survey package
          print("Running svy_ests")
          print(Sys.time())
          J <- numclusters # number of clusters in the population
          svy_ests(J, num.clusters, num.units, rootdir, sim, popseed, use.sizes)
          print("##################################################################################")
          print("##################################################################################")
          # delete files
          if (num.units <= 1) {
            nunits <- paste(100*num.units, "pct", sep = "")
          } else {
            nunits <- num.units
          }
          fil1 <- paste(rootdir, "/Data/Simplify/vary_K/sampledata_usesizes_", use.sizes, "_nclusters_", num.clusters,
                        "_nunits_", nunits, "_sim_", sim, ".RData", sep = "")
          if (file.exists(fil1)) {
            file.remove(fil1)
          }
        } # end unit loop
      } # end cluster loop
    } # end use sizes loop

    #stopCluster(cl) # close workers

    # remove the pop data file with the indices
      file.remove(paste(rootdir, "/Data/Simplify/vary_K/popdata_usesizes_", use.sizes,
				 "_inds_nclusters_", num.clusters, "_nunits_", nunits,
				 "_sim_", sim, ".RData", sep = ""))
} # end function
