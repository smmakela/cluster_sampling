# Author: Susanna Makela
# Date: 21 Apr 2014
# Purpose: run stan

runstan <- function(num.clusters, num.units, rootdir, sim, popseed, stanmod, stanmod_name, use.sizes) {
  # num.clusters -- number of clusters to sample
  # num.units -- number of units to sample
  # rootdir -- root directory where Code, Data folders are
  # sim -- current iteration; used so that multiple instances aren't trying to write to the same file
  # popseed -- which seed was used to make the pop data
  # stanmod -- compiled stan model
  # stanmod_name -- string for name of stan model so we know which parts of the code to run
  # use.sizes -- 0/1 for whether cluster sizes used in pop data

  ##########################################
  ### Setup of libraries and directories
  ##########################################
    library(Rcpp)
    library(inline)
    library(rstan)
    library(plyr)    
    library(dplyr)
    if (num.units <= 1) {
      nunits <- paste(100*num.units, "pct", sep = "")
    } else {
      nunits <- num.units
    }
    load(paste(rootdir, "/Data/Simplify/vary_K/popdata_usesizes_", use.sizes, "_inds_nclusters_", num.clusters,
               "_nunits_", nunits, "_sim_", sim, ".RData", sep = ""))
    load(paste(rootdir, "/Data/Simplify/vary_K/popdata_usesizes_", use.sizes, "_seed_", popseed, ".RData", sep = ""))
    load(paste(rootdir, "/Data/Simplify/vary_K/sampledata_usesizes_", use.sizes, "_nclusters_", num.clusters,
               "_nunits_", nunits, "_sim_", sim, ".RData", sep = ""))
    #nms <- c("pop.data", "Mj", "J")
    #for (j in nms) {
    #  assign(j, popdata[[j]])
    #}
    J <- popdata[["J"]]
    rm(popdata)
    #rm(pop.data)
    pop.data <- pop.data.inds
 
  ##########################################
  ### Renumber cluster ids in pop data so that 1:J_sam are the sampled clusters and J_sam+1:J_pop are the unsampled ones
  ##########################################
    #sampled.cluster.ids <- sort(unique(sample.data$cluster.id))
    #J_sam <- length(sampled.cluster.ids)
    #J_pop <- J
    #all.cluster.ids <- sort(unique(pop.data$cluster.id))
    #nonsampled.cluster.ids <- sort(setdiff(all.cluster.ids, sampled.cluster.ids))
    #idmap <- data.frame(cluster.id = c(sampled.cluster.ids, nonsampled.cluster.ids), new.cluster.id2 = c(1:J))
    #pop.data <- merge(pop.data, idmap, by = "cluster.id")
    #pop.data$new.cluster.id <- pop.data$new.cluster.id2
    #pop.data$new.cluster.id2 <- NULL

  ##########################################
  ### Sort pop and sample data by new.cluster.id
  ##########################################
    #pop.data <- tbl_df(pop.data)
    #sample.data <- tbl_df(sample.data)
    #pop.data <- arrange(pop.data, new.cluster.id)
    #sample.data <- arrange(sample.data, new.cluster.id)

   
  ##########################################
  ### Make data for stan
  ##########################################
    print("making stan data")
    print(Sys.time())

    # SORT pop and sample data by cluster.id
    pop.data <- tbl_df(pop.data)
    pop.data <- arrange(pop.data, cluster.id)
    sample.data <- tbl_df(sample.data)
    sample.data <- arrange(sample.data, cluster.id)

    J_pop <- J
    J_sam <- num.clusters

    # pop data at cluster level
print(names(pop.data))
    pop.data.clev <- distinct(pop.data, cluster.id, Mj, logMj_c)
    Mj_pop <- pop.data.clev$Mj
    N_pop <- sum(Mj_pop)
    Tx <- N_pop
print("ONE")
    # sample data at cluster level
    sam.dat <- filter(pop.data, insample == 1)
    sample.data.clev <- distinct(sam.dat, cluster.id, Mj, logMj_c)
print("TWO")
    Mj_sam <- sample.data.clev$Mj
    N_sam <- sum(pop.data$insample)
    logMj_sam <- sample.data.clev$logMj_c
    sample.data.clev2 <- summarise(group_by(sample.data, cluster.id),
                                   ybar_sam = mean(yi))
    ybar_sam <- sample.data.clev2$ybar_sam
    # special data for bb:
    #  counts of unique values of Mj in sample data
    #  number of unique values of Mj in sample data
    n.dat <- summarise(group_by(sample.data.clev, Mj),
                       n = n_distinct(cluster.id))
print("THREE")
print(names(sample.data))
    n <- n.dat$n
    M <- n_distinct(sample.data$Mj)
print("FOUR")

    # sample data
    xi_sam <- sample.data$xi
    yi_sam <- sample.data$yi
    nj.dat <- summarise(group_by(sample.data, cluster.id),
                        nj = n())
    nj <- nj.dat$nj

    # Mj in nonsampled clusters
    nonsamp.dat <- filter(pop.data, insample == 0)
    nonsamp.dat <- arrange(nonsamp.dat, cluster.id)
    mis.dat <- distinct(nonsamp.dat, cluster.id, Mj)
    Mj_mis <- mis.dat$Mj

    # xbar_pop: in unsampled clusters, mean(xi) for all units
    #   in the cluster; for sampled clusters, mean(xi) for UNSAMPLED
    #   units in the cluster
    xbar.dat <- summarise(group_by(nonsamp.dat, cluster.id),
                          xbar_pop = mean(xi))
    xbar_pop <- xbar.dat$xbar_pop

    # cluster ids, of length N_sam
print(dim(sample.data))
print(N_sam)
    cluster_id_long <- sample.data$cluster.id

    #new.cluster.id.rep <- sample.data$new.cluster.id # length = sample size
    #new.cluster.id <- sort(unique(sample.data$new.cluster.id)) # length = number of sampled clusters
    #N_pop <- sum(Mj)
    #N_sam <- nrow(sample.data)
    #xi_sam <- sample.data$xi
    #yi_sam <- sample.data$yi
    #xi_pop <- pop.data$xi
    #logMj_pop_c <- log(Mj) - mean(log(Mj))
    #logMj_pop_c <- log(Mj)
    
    # Mj, logMj in sample data
    #tmp <- unique(sample.data[, c("new.cluster.id", "logMj_c", "Mj")])
    #tmp2 <- tmp[with(tmp, order(new.cluster.id)),]
    #logMj_sam_c <- tmp2$logMj_c
    #Mj_sam <- tmp2$Mj

    # xbar_pop -- in sampled clusters, this is actually xbar for the UNsampled units
    #tmp <- ddply(pop.data[,c("new.cluster.id", "xi")], .(new.cluster.id), summarise, mean = mean(xi))
    #tmp <- tmp[with(tmp, order(new.cluster.id)),]
    #xbar_pop <- tmp$mean # elements 1:J_sam are the sampled clusters, J_sam+1:J_pop the unsampled ones
    #if (num.units != 1) { # when we've sampled all the units (num.units = 1 meaning 100%), don't need to do this
    #  tmp <- ddply(pop.data[pop.data$cluster.id %in% sampled.cluster.ids & pop.data$insample == 0,],
    #               .(new.cluster.id), xbar = mean(xi), summarise) # mean of xi for unsampled units in sampled clusters
    #  xbar_pop[1:J_sam] <- tmp$xbar
    #}

    # Mj_mis -- number of nonsampled units in each cluster
    #tmp <- ddply(sample.data, .(new.cluster.id), nj = length(new.cluster.id), summarise)
    #nj <- tmp$nj
    #Mj_mis <- Mj
    #Mj_mis[1:J] <- Mj - tmp$nj
 
    # xbar_sam, ybar_sam, Tx
    #tmp <- ddply(sample.data, .(new.cluster.id), summarise, mean_xi = mean(xi), mean_yi = mean(yi))
    #tmp <- tmp[with(tmp, order(new.cluster.id)),]
    #xbar_sam <- tmp$mean_xi
    #ybar_sam <- tmp$mean_yi
    #Tx <- sum(exp(logMj_pop_c)) # because everything is with centered LOG sizes
    
    # M (unique Mj's in sampled data) and n (vector of counts for unique Mj's)
    #tmp <- ddply(sample.data, .(Mj), summarise, ct = length(unique(new.cluster.id)))
    #tmp <- tmp[!is.na(tmp$Mj), ]
    #M <- nrow(tmp) # number of unique Mj's
    #n <- tmp$ct
    #if (sum(n) != num.clusters) {
    #  stop(paste("sum(n) != ", num.clusters, "!", sep = ""))
    #}
    if (stanmod_name == "bb") {
      standata <- list(J_pop = J_pop,
                       J_sam = J_sam,
                       N_sam = N_sam,
                       M = M,
                       cluster_id_long = cluster_id_long,
                       xi_sam = xi_sam,
                       logMj_sam = logMj_sam,
                       yi_sam = yi_sam,
                       ybar_sam = ybar_sam,
                       xbar_pop = xbar_pop,
                       Mj_mis = Mj_mis,
                       Mj_sam = Mj_sam,
                       Tx = Tx,
                       n = n)
    } else if (stanmod_name == "knowsizes") {
      standata <- list(J_pop = J_pop,
                       J_sam = J_sam,
                       N_pop = N_pop,
                       N_sam = N_sam,
                       cluster_id_long = cluster_id_long,
                       xi_sam = xi_sam,
                       yi_sam = yi_sam,
                       Mj_pop = Mj_pop,
                       logMj_sam = logMj_sam,
                       Mj_sam = Mj_sam,
                       Mj_mis = Mj_mis,
                       ybar_sam = ybar_sam,
                       xbar_pop = xbar_pop)
    } else if (stanmod_name == "negbin") {
      standata <- list(J_pop = J_pop,
                       J_sam = J_sam,
                       N_sam = N_sam,
                       cluster_id_long = cluster_id_long,
                       xi_sam = xi_sam, 
                       logMj_sam = logMj_sam,
                       yi_sam = yi_sam,
                       ybar_sam = ybar_sam,
                       xbar_pop = xbar_pop,
                       Mj_mis = Mj_mis,
                       nj = nj,
                       Mj_sam = Mj_sam)
    } else { # stanmod_name == "cluster_inds_only"
      standata <- list(J_pop = J_pop,
                       J_sam = J_sam,
                       N_pop = N_pop,
                       N_sam = N_sam,
                       cluster_id_long = cluster_id_long,
                       Mj_pop = Mj_pop,
                       Mj_sam = Mj_sam,
                       Mj_mis = Mj_mis,
                       xi_sam = xi_sam,
                       yi_sam = yi_sam,
                       ybar_sam = ybar_sam,
                       xbar_pop = xbar_pop)
    }

#    if (stanmod_name == "bb") {
#      standata <- list(J_pop = J_pop, J_sam = J_sam,
#                       N_sam = N_sam, M = M,
#                       new_cluster_id_rep = new.cluster.id.rep,
#                       xi_sam = xi_sam, logMj_sam = logMj_sam_c,
#                       yi_sam = yi_sam, ybar_sam = ybar_sam,
#                       xbar_pop = xbar_pop, Mj_mis = Mj_mis,
#                       Mj_sam = Mj_sam, nj = nj, Tx = Tx, n = n)
#    } else if (stanmod_name == "knowsizes") {
#      standata <- list(J_pop = J_pop, J_sam = J_sam,
#                       N_pop = N_pop, N_sam = N_sam,
#                       new_cluster_id_rep = new.cluster.id.rep,
#                       new_cluster_id = new.cluster.id,
#                       xi_sam = xi_sam, yi_sam = yi_sam,
#                       Mj_pop = Mj, logMj_sam = logMj_sam_c,
#                       Mj_sam = Mj_sam, Mj_mis = Mj_mis,
#                       ybar_sam = ybar_sam, xbar_pop = xbar_pop)
#    } else if (stanmod_name == "negbin") {
#      standata <- list(J_pop = J_pop, J_sam = J_sam,
#                       N_pop = N_pop, N_sam = N_sam,
#                       new_cluster_id_rep = new.cluster.id.rep,
#                       new_cluster_id = new.cluster.id,
#                       xi_sam = xi_sam, yi_sam = yi_sam,
#                       Mj_pop = Mj, logMj_sam = logMj_sam_c,
#                       Mj_sam = Mj_sam, Mj_mis = Mj_mis, nj = nj,
#                       ybar_sam = ybar_sam, xbar_pop = xbar_pop)
#    } else { # stanmod_name == "cluster_inds_only"
#      standata <- list(J_pop = J_pop, J_sam = J_sam,
#                       N_pop = N_pop, N_sam = N_sam,
#                       new_cluster_id_rep = new.cluster.id.rep,
#                       new_cluster_id = new.cluster.id,
#                       Mj_pop = Mj,
#                       Mj_sam = Mj_sam, Mj_mis = Mj_mis,
#                       xi_sam = xi_sam, yi_sam = yi_sam,
#                       ybar_sam = ybar_sam, xbar_pop = xbar_pop)
#    }
print(str(standata))
print(summary(Mj_sam))
  ##########################################
  ### Run stan
  ##########################################
    print("***********************************************************************")
    print(paste("********** about to run stan for num.clusters = ", num.clusters, 
                ", num.units = ", nunits, ", model ", stanmod_name, sep = ""))
    print(Sys.time())
    fit <- sampling(stanmod, data = standata, iter = 5000, chains = 4)
    if (sim %in% c(1, 100, 200, 300, 400, 500)) {
      save(fit, file = paste(rootdir, "/Results/Simplify/vary_K/stanfit_usesizes_", use.sizes, "_nclusters_", num.clusters,
                               "_nunits_", nunits, "_sim_", sim, "_", stanmod_name, ".RData", sep = ""))
      print("done saving stanfit object")
      print(Sys.time())
    }
    print("done fitting stan model")
    print(Sys.time())

    if (stanmod_name == "bb") {
      parlist1 <- c("gamma0", "gamma1", "alpha0", "alpha1", "sigma_beta0", "sigma_beta1", "sigma_y",
                 "ybar_hat")
      #parlist2 <- c("n_star", "phi", "phi_star", "Mj_unsamp", "Mj_all")
      parlist2 <- c("n_star", "phi", "phi_star", "Mj_tot_est")
      samps <- extract(fit, permute = FALSE, pars = c(parlist1, parlist2))
    } else if (stanmod_name == "knowsizes") {
      parlist <- c("gamma0", "gamma1", "alpha0", "alpha1", "sigma_beta0", "sigma_beta1", "sigma_y",
                 "ybar_hat")
      samps <- extract(fit, permute = FALSE, pars = parlist)
    } else if (stanmod_name == "negbin") {
      parlist1 <- c("gamma0", "gamma1", "alpha0", "alpha1", "sigma_beta0", "sigma_beta1", "sigma_y",
                 "ybar_hat")
      parlist2 <- c("mu_star", "phi_star", "mu", "phi", "Mj_tot_est")
      samps <- extract(fit, permute = FALSE, pars = c(parlist1, parlist2))
    } else {
      parlist <- c("gamma0", "alpha0", "sigma_beta0", "sigma_beta1", "sigma_y",
                 "ybar_hat")
      samps <- extract(fit, permute = FALSE, pars = parlist)
    }
    print("done making samps")
    print(Sys.time())

    tt <- summary(fit)$summary

    if (stanmod_name == "bb" | stanmod_name == "negbin") {
      par.ests1 <- tt[parlist1, ]
      inds <- NULL
      for (j in 1:length(parlist2)) {
        inds <- c(inds, grep(parlist2[j], rownames(tt)))
      }
      inds <- unique(inds) 
      par.ests2 <- tt[inds, ]
      par.ests <- rbind(par.ests1, par.ests2)
    } else {
      par.ests <- tt[parlist, ]
    }
    rm(fit)
    print("done making par.ests")
    print(Sys.time())

    if (sim %in% c(100)) {
      save(samps, file = paste(rootdir, "/Results/Simplify/vary_K/samples_usesizes_", use.sizes, "_nclusters_", num.clusters,
                               "_nunits_", nunits, "_sim_", sim, "_", stanmod_name, ".RData", sep = ""))
      print("done saving samps")
      print(Sys.time())
    }

    print("making par.ests")
    print(Sys.time())
    par.ests.rownames <- attr(par.ests, "dimnames")[[1]]
    par.ests.rownames <- gsub("\\[", "", par.ests.rownames) 
    par.ests.rownames <- gsub("\\]", "", par.ests.rownames) 
    par.ests.colnames <- attr(par.ests, "dimnames")[[2]]
    par.ests <- data.frame(par.ests, row.names = par.ests.rownames)
    colnames(par.ests) <- par.ests.colnames

    # save par.ests
    write.table(par.ests,
                file = paste(rootdir, "/Results/Simplify/vary_K/parests_stan_usesizes_", use.sizes, "_nclusters_", num.clusters,
                             "_nunits_", nunits, "_sim_", sim, "_", stanmod_name, ".txt", sep = ""), sep = ",")
    
    print("printing results")
    print(Sys.time())      
    ybar.true <- mean(pop.data$yi)
    ybar.hat <- par.ests["ybar_hat",1]
    ybar.hat.lci50 <- par.ests["ybar_hat", 5] # 25% quantile
    ybar.hat.uci50 <- par.ests["ybar_hat", 7] # 75% quantile
    ybar.hat.lci95 <- par.ests["ybar_hat", 4] # 2.5% quantile
    ybar.hat.uci95 <- par.ests["ybar_hat", 8] # 97.5% quantile
    rnames <- c("gamma0", "gamma1", "alpha0", "alpha1", "sigma_beta0", "sigma_beta1", "sigma_y")
    cnames <- c("mean", "2.5%", "97.5%", "Rhat", "n_eff")
    cat("########################################################\n")
    cat("-----------------  use.sizes = ", use.sizes, ", stanmod = ", stanmod_name,
        ", num.clusters = ", num.clusters, ", num.units = ", num.units, "---------------------\n")
    cat("true: ", ybar.true, "\n")
    cat("est: ", ybar.hat, "\n")
    cat("est 95% CI: ", ybar.hat.lci95, ", ", ybar.hat.uci95, "\n")
    cat("param posterior means: \n")
    print(par.ests[rnames, cnames])
    cat("summary of Mj_sam: \n")
    print(summary(Mj_sam))
    cat("summary of Mj_unsamp posterior means:\n")
    print(summary(par.ests[grep("Mj_unsamp", rownames(par.ests)), "mean"]))
    cat("########################################################\n")


   
    print("saving results")
    print(Sys.time())
    res <- c(ybar.true, ybar.hat, ybar.hat.lci50, ybar.hat.uci50, ybar.hat.lci95, ybar.hat.uci95) 
    write.table(res,
                file = paste(rootdir, "/Results/Simplify/vary_K/ybar_stan_usesizes_", use.sizes, "_nclusters_", num.clusters,
                             "_nunits_", nunits, "_sim_", sim, "_", stanmod_name, ".txt", sep = ""), sep = ",")


  return(res)
}
