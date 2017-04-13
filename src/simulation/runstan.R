# Author: Susanna Makela
# Date: 21 Apr 2014
# Purpose: run stan

runstan <- function(num.clusters, num.units, use.sizes, rootdir, sim,
                    stanmod, stanmod_name, num.iter, num.chains) {
  # num.clusters -- number of clusters to sample
  # num.units -- number of units to sample
  # use.sizes -- 0/1 for whether cluster sizes used in pop data
  # rootdir -- root directory where Code, Data folders are
  # sim -- current iteration; used so that multiple instances aren't trying to write to the same file
  # stanmod -- compiled stan model
  # stanmod_name -- string for name of stan model so we know which parts of the code to run
  # num.iter -- number of iterations stan should run for
  # num.chains -- number of chains to run in stan

  ##########################################
  ### Load pop and sample data
  ##########################################
    if (num.units <= 1) {
      nunits <- paste(100*num.units, "pct", sep = "")
    } else {
      nunits <- num.units
    }
    simdata <- readRDS(paste0(rootdir, "output/simulation/simdata_usesizes_",
                              use.sizes, "_", outcome.type, "_nclusters_",
                              num.clusters, "_nunits_", nunits, "_simno_", sim,
                              ".rds"))
    for (j in names(simdata)) {
      assign(j, simdata[[j]])
    }
    rm(simdata)
    ybar.true <- mean(pop.data$y)
   
  ##########################################
  ### Make data for stan
  ##########################################
    print("making stan data")
    print(Sys.time())

    # SORT pop and sample data by cluster.id
    pop.data <- dplyr::arrange(pop.data, cluster.id)
    sample.data <- dplyr::arrange(sample.data, cluster.id)

    J_pop <- J
    if (num.units != 999) {
      J_sam <- num.clusters
    } else {
      J_sam <- J_pop
    }
    Mj_pop <- Mj
    N_pop <- sum(Mj_pop)
    Tx <- N_pop
    N_sam <- sum(pop.data$insample)
    tmp <- dplyr::summarise(group_by(pop.data, cluster.id),
                            xbar_pop = mean(x))
    tmp <- dplyr::arrange(tmp, cluster.id)
    xbar_pop <- tmp$xbar_pop
    x <- sample.data$x
    y <- sample.data$y
    cluster_id_long <- sample.data$cluster.id

    # get sample data at cluster level
    sam.dat <- dplyr::filter(pop.data, insample == 1)
    sam.dat <- dplyr::distinct(sam.dat, cluster.id, Mj, logMj_c)
    Mj_sam <- sam.dat$Mj
    logMj_sam <- sam.dat$logMj_c
    # special data for bb:
    #  counts of unique values of Mj in sample data
    #  number of unique values of Mj in sample data
    n.dat <- summarise(group_by(sam.dat, Mj),
                       n = n_distinct(cluster.id))
    n <- n.dat$n
    M <- n_distinct(sample.data$Mj)

    # delete pop and sample data to save memory
    rm(pop.data)
    rm(sample.data)

    # make standata lists
    if (stanmod_name == "bb") {
      standata <- list(J_pop = J_pop,
                       J_sam = J_sam,
                       N_sam = N_sam,
                       M = M,
                       cluster_id_long = cluster_id_long,
                       x = x,
                       y = y,
                       Mj_sam = Mj_sam,
                       logMj_sam = logMj_sam,
                       xbar_pop = xbar_pop,
                       Tx = Tx,
                       n = n)
    } else if (grepl("knowsizes", stanmod_name) | stanmod_name == "cluster_inds_only") {
      standata <- list(J_pop = J_pop,
                       J_sam = J_sam,
                       N_sam = N_sam,
                       cluster_id_long = cluster_id_long,
                       x = x,
                       y = y,
                       Mj_sam = Mj_sam,
                       logMj_sam = logMj_sam,
                       Mj_pop = Mj_pop,
                       xbar_pop = xbar_pop)
    } else { #(stanmod_name == "negbin")
      standata <- list(J_pop = J_pop,
                       J_sam = J_sam,
                       N_sam = N_sam,
                       cluster_id_long = cluster_id_long,
                       x = x, 
                       y = y,
                       Mj_sam = Mj_sam,
                       logMj_sam = logMj_sam,
                       xbar_pop = xbar_pop)
    } 
print(str(standata))
  ##########################################
  ### Run stan
  ##########################################
    print("***********************************************************************")
    print(paste("********** about to run stan for num.clusters = ", num.clusters, 
                ", num.units = ", nunits, ", model ", stanmod_name, sep = ""))
    print(Sys.time())
    fit <- sampling(stanmod, data = standata, iter = num.iter, chains = num.chains,
                    control = list(stepsize = 0.001, adapt_delta = 0.999))
    if (sim %in% c(1, 200, 300, 400, 500, 999)) {
      save(fit, file = paste0(rootdir, "/output/simulation/stanfit_usesizes_",
                              use.sizes, "_nclusters_", num.clusters,
                              "_nunits_", nunits, "_simno_", sim, "_",
                              stanmod_name, ".RData"))
      print("done saving stanfit object")
      print(Sys.time())
    }
    print("done fitting stan model")
    print(Sys.time())

    if (stanmod_name == "bb") {
      parlist1 <- c("alpha0", "gamma0", "alpha1", "gamma1",
                    "sigma_beta0", "sigma_beta1", "sigma_y", "ybar_new")
      parlist2 <- c("n_star", "phi", "phi_star", "M_tot_est")
      samps <- extract(fit, permute = FALSE, pars = c(parlist1, parlist2))
    } else if (grepl("knowsizes", stanmod_name) & stanmod_name != "knowsizes_ccc") {
      parlist <- c("alpha0", "gamma0", "alpha1", "gamma1",
                   "sigma_beta0", "sigma_beta1", "sigma_y", "ybar_new")
    } else if (stanmod_name == "knowsizes_ccc") {
      parlist <- c("alpha0", "gamma0",
                   "sigma_beta0", "sigma_y", "ybar_new")
      samps <- extract(fit, permute = FALSE, pars = parlist)
    } else if (stanmod_name == "negbin") {
      parlist1 <- c("alpha0", "gamma0", "alpha1", "gamma1",
                    "sigma_beta0", "sigma_beta1", "sigma_y", "ybar_new")
      parlist2 <- c("mu_star", "phi_star", "mu", "phi", "M_tot_est")
      samps <- extract(fit, permute = FALSE, pars = c(parlist1, parlist2))
    } else if (stanmod_name == "lognormal") {
      parlist1 <- c("alpha0", "gamma0", "alpha1", "gamma1",
                    "sigma_beta0", "sigma_beta1", "sigma_y", "ybar_new")
      parlist2 <- c("mu_sb", "sigma_M", "mu", "M_tot_est")
      samps <- extract(fit, permute = FALSE, pars = c(parlist1, parlist2))
    } else {
      parlist <- c("alpha0", "alpha1",
                   "sigma_beta0", "sigma_beta1", "sigma_y", "ybar_new")
      samps <- extract(fit, permute = FALSE, pars = parlist)
    }
    print("done making samps")
    print(Sys.time())

#pdf(paste0(rootdir, "/output/simulation/pairs_usesizes_", use.sizes,
#           "_", stanmod_name, ".pdf"), width = 10, height = 6, onefile = TRUE)
#pairs(fit, pars = c(parlist1, parlist2))
#dev.off()

    tt <- summary(fit)$summary

    if (stanmod_name == "bb" | stanmod_name == "negbin" | stanmod_name == "lognormal") {
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

    #if (sim %in% c(100)) {
    #  save(samps, file = paste(rootdir, "/Results/Simplify/vary_K/samples_usesizes_", use.sizes, "_nclusters_", num.clusters,
    #                           "_nunits_", nunits, "_sim_", sim, "_", stanmod_name, ".RData", sep = ""))
    #  print("done saving samps")
    #  print(Sys.time())
    #}

    print("making par.ests")
    print(Sys.time())
    par.ests.rownames <- attr(par.ests, "dimnames")[[1]]
    par.ests.rownames <- gsub("\\[", "", par.ests.rownames) 
    par.ests.rownames <- gsub("\\]", "", par.ests.rownames) 
    par.ests.colnames <- attr(par.ests, "dimnames")[[2]]
    par.ests <- data.frame(par.ests, row.names = par.ests.rownames)
    colnames(par.ests) <- par.ests.colnames

    # save par.ests
    #write.table(par.ests,
    #            file = paste0(rootdir,
    #                          "/output/simulation/parests_stan_usesizes_",
    #                          use.sizes, "_nclusters_", num.clusters,
    #                          "_nunits_", nunits, "_sim_", sim, "_",
    #                          stanmod_name, ".txt"),
    #            sep = ",")
    
    print("printing results")
    print(Sys.time())      
print(str(par.ests))
    ybar.hat <- par.ests["ybar_new", "mean"]
    ybar.hat.lci50 <- par.ests["ybar_new", "25%"] # 25% quantile
    ybar.hat.uci50 <- par.ests["ybar_new", "75%"] # 75% quantile
    ybar.hat.lci95 <- par.ests["ybar_new", "2.5%"] # 2.5% quantile
    ybar.hat.uci95 <- par.ests["ybar_new", "97.5%"] # 97.5% quantile
    #rnames <- c("alpha0", "gamma0", "alpha1", "gamma1", "sigma_beta0",
    #            "sigma_beta1", "sigma_y")
    rnames <- c("alpha0", "gamma0", "sigma_beta0", "sigma_y")
    cnames <- c("mean", "2.5%", "97.5%", "Rhat", "n_eff")
    cat("########################################################\n")
    cat("-----------------  use.sizes = ", use.sizes, ", stanmod = ",
        stanmod_name, ", num.clusters = ", num.clusters, ", num.units = ",
        num.units, "---------------------\n")
    cat("true: ", ybar.true, "\n")
    cat("est: ", ybar.hat, "\n")
    cat("est 95% CI: ", ybar.hat.lci95, ", ", ybar.hat.uci95, "\n")
    cat("param posterior means: \n")
    print(par.ests[rnames, cnames])
    cat("summary of Mj_sam: \n")
    print(summary(Mj_sam))
    cat("summary of Mj_unsamp posterior means:\n")
    print(summary(par.ests[grep("Mj_mis", rownames(par.ests)), "mean"]))
    cat("########################################################\n")

    print("saving results")
    print(Sys.time())
    ybar.ests <- data.frame(ybar.true, ybar.hat, ybar.hat.lci50, ybar.hat.uci50,
                            ybar.hat.lci95, ybar.hat.uci95) 
    #write.table(res,
    #            file = paste0(rootdir, "/output/simulation/ybar_stan_usesizes_",
    #                          use.sizes, "_nclusters_", num.clusters,
    #                          "_nunits_", nunits, "_sim_", sim, "_",
    #                          stanmod_name, ".txt"),
    #            sep = ",")

  toreturn <- list(par.ests = par.ests, ybar.ests = ybar.ests)
  return(toreturn)
}
