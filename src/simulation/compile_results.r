# Author: Susanna Makela
# Date: 13 Jan 2016
# Purpose: process results from simulation

##########################################
### Setup of directories and libraries
##########################################
  libdir <- "/vega/stats/users/smm2253/rpackages"
  .libPaths(libdir)
  rootdir <- "/vega/stats/users/smm2253/cluster_sampling/"
  Sys.setenv(HOME = rootdir)

  # set working directory, figure directory
  resdir <- paste0(rootdir, "output/simulation/", sep = "")
  setwd(resdir)

  # libraries
  library(plyr)

  # print time
  print(Sys.time())

  # calculate true Mj total
    load(file = paste(rootdir, "/Data/Simplify/vary_K/popdata_usesizes_0_seed_1.RData", sep = ""))
    Mj_tot_true_0 <- sum(popdata[["Mj"]])
    load(file = paste(rootdir, "/Data/Simplify/vary_K/popdata_usesizes_1_seed_1.RData", sep = ""))
    Mj_tot_true_1 <- sum(popdata[["Mj"]])

  # calculate true value of ybar
    ybar_true <- mean(popdata[["pop.data"]]$yi)

  # list of unit and cluster numbers used
    #unitlist <- c(10, 50, 100)
    #num.units.list <- c(.05, .1, .25, .5, 1, 10, 50, 100)
    #num.units.list2 <- c("5pct", "10pct", "25pct", "50pct", "100pct", 10, 50, 100)
    num.units.list <- c(.05, .1, 10, 50)
    num.units.list2 <- c("5pct", "10pct", "10", "50")
    #cluster.list <- c(15, 30, 60, 100)
    cluster.list <- c(15, 30, 60)
    #num.units.list <- c(.1, 10)
    #num.units.list2 <- c("10pct", "10")
    #cluster.list <- c(15, 30)
    #num.units.list <- .1
    #num.units.list2 <- "10pct"
    #cluster.list <- 15
    u.list <- c(0, 1)

  # which model to pull results for
    modelname.list <- c("knowsizes", "cluster_inds_only", "negbin", "bb")
    #modelname.list <- c("negbin", "bb")

##########################################
### Loop through res files -- these are the stan results
##########################################
cat("--------------------------- ybar_stan -------------------------\n")
  fil.list <- NULL
  for (i in 1:length(cluster.list)) {
    for (j in 1:length(num.units.list)) {
      for (k in 1:length(modelname.list)) {
        for (u in u.list) {
          fil.list <- c(fil.list,
                        list.files(resdir,
                                   paste("ybar_stan_usesizes_", u, "_nclusters_", cluster.list[i],
                                         "_nunits_", num.units.list2[j],
                                         "_sim_.*_", modelname.list[k], ".txt", sep = "")))
        }
      }
    }
  }
  nres <- length(fil.list)
  ybar_stan <- data.frame()
  for (f in 1:nres) {
    filnam <- fil.list[f]
    tt <- file.info(paste(resdir, "/", filnam, sep = ""))
    fdate <- tt$mtime
    fdate <- substr(fdate, 1, 10)
    fsize <- tt$size
    if (!(fdate %in% resdate) | (fsize == 0 & !is.na(fsize))) {
      if (fsize == 0 & !is.na(fsize)) {
        print(paste("********** Output file of size 0 for ", filnam, " ****************", sep = ""))
      }
      next
    }
    fil <- read.table(filnam, sep = ",")
    ybar.true <- fil[1,1]
    ybar.hat <- fil[2,1]
    ybar.lci50 <- fil[3,1]
    ybar.uci50 <- fil[4,1]
    ybar.lci95 <- fil[5,1]
    ybar.uci95 <- fil[6,1]
    ff <- unlist(strsplit(filnam, "[[:punct:]]"))
    nclusters <- as.numeric(ff[6])
    nunits <- ff[8]
    simno <- as.numeric(ff[10])
    modelname <- ff[11]
    use.sizes <- as.numeric(ff[4])
    covers.50 <- ybar.lci50 <= ybar.true & ybar.true <= ybar.uci50
    covers.95 <- ybar.lci95 <= ybar.true & ybar.true <= ybar.uci95
    df <- data.frame(ybar.true, ybar.hat, ybar.lci50, ybar.uci50, ybar.lci95, ybar.uci95,
                     covers.50, covers.95, simno, nclusters, nunits, modelname, use.sizes)
    ybar_stan <- rbind(ybar_stan, df)
  }
  rm(df)

##########################################
### Loop through svy files
##########################################
cat("--------------------------- ybar_svy -------------------------\n")
  fil.list <- NULL
  for (i in 1:length(cluster.list)) {
    for (j in 1:length(num.units.list)) {
      for (u in u.list) {
        fil.list <- c(fil.list,
                      list.files(resdir,
                                 paste("ybar_svy_usesizes_", u, "_nclusters_", cluster.list[i],
                                       "_nunits_", num.units.list2[j],
                                       "_sim_.*.txt", sep = "")))
      }
    }
  }
  nres <- length(fil.list)
  ybar_svy <- data.frame()
  for (f in 1:nres) {
    filnam <- fil.list[f]
    tt <- file.info(paste(resdir, "/", filnam, sep = ""))
    fdate <- tt$mtime
    fdate <- substr(fdate, 1, 10)
    fsize <- tt$size
    if (!(fdate %in% resdate) | (fsize == 0 & !is.na(fsize))) {
      if (fsize == 0 & !is.na(fsize)) {
        print(paste("********** Output file of size 0 for ", filnam, " ****************", sep = ""))
      }
      next
    }
    fil <- read.csv(filnam)
    names(fil) <- c("simno", "ybar_true", "ybar_hat", "ybar_se", "ybar_hat2", "ybar_se2")
    fil$parname <- rownames(fil)
    rownames(fil) <- NULL
    ff <- unlist(strsplit(filnam, "[[:punct:]]"))
    nclusters <- as.numeric(ff[6])
    nunits <- ff[8]
    simno <- as.numeric(ff[10])
    fil$nclusters <- nclusters
    fil$nunits <- nunits
    fil$simno <- simno
    fil$use.sizes <- as.numeric(ff[4])
    fil$lci.50.ht <- fil$ybar_hat - fil$ybar_se
    fil$uci.50.ht <- fil$ybar_hat + fil$ybar_se
    fil$lci.95.ht <- fil$ybar_hat - 1.96*fil$ybar_se
    fil$uci.95.ht <- fil$ybar_hat + 1.96*fil$ybar_se
    fil$covers.50.ht <- fil$lci.50.ht <= fil$ybar_true & fil$ybar_true <= fil$uci.50.ht
    fil$covers.95.ht <- fil$lci.95.ht <= fil$ybar_true & fil$ybar_true <= fil$uci.95.ht
    fil$lci.50.gr <- fil$ybar_hat2 - fil$ybar_se2
    fil$uci.50.gr <- fil$ybar_hat2 + fil$ybar_se2
    fil$lci.95.gr <- fil$ybar_hat2 - 1.96*fil$ybar_se2
    fil$uci.95.gr <- fil$ybar_hat2 + 1.96*fil$ybar_se2
    fil$covers.50.gr <- fil$lci.50.gr <= fil$ybar_true & fil$ybar_true <= fil$uci.50.gr
    fil$covers.95.gr <- fil$lci.95.gr <= fil$ybar_true & fil$ybar_true <= fil$uci.95.gr
    ybar_svy <- rbind(ybar_svy, fil)
  }
  rm(fil)

print(str(ybar_stan))
 print(str(ybar_svy)) 

##########################################
### Combine stan and svy ests of ybar
##########################################
  namlist.ht <- c("ybar_true", "ybar_hat", "lci.50.ht", "uci.50.ht", "lci.95.ht", "uci.95.ht",
                  "covers.50.ht", "covers.95.ht", "nunits", "nclusters", "simno", "use.sizes")
  namlist.gr <- c("ybar_true", "ybar_hat2", "lci.50.gr", "uci.50.gr", "lci.95.gr", "uci.95.gr",
                  "covers.50.gr", "covers.95.gr", "nunits", "nclusters", "simno", "use.sizes")
  ybar_svy_ht <- ybar_svy[, namlist.ht]
  names(ybar_svy_ht) <- c("true", "hat", "lci50", "uci50", "lci95", "uci95",
                          "covers50", "covers95", "nunits", "nclusters", "simno", "use.sizes")
  ybar_svy_ht$modelname <- "svy_ht"
  ybar_svy_gr <- ybar_svy[, namlist.gr]
  names(ybar_svy_gr) <- c("true", "hat", "lci50", "uci50", "lci95", "uci95",
                          "covers50", "covers95", "nunits", "nclusters", "simno", "use.sizes")
  ybar_svy_gr$modelname <- "svy_gr"
  ybar_stan2 <- ybar_stan[, c("ybar.true", "ybar.hat", "ybar.lci50", "ybar.uci50", "ybar.lci95", "ybar.uci95",
                              "covers.50", "covers.95", "nunits", "nclusters", "simno", "use.sizes", "modelname")]
  names(ybar_stan2) <- c("true", "hat", "lci50", "uci50", "lci95", "uci95", "covers50", "covers95",
                         "nunits", "nclusters", "simno", "use.sizes", "modelname")
#  ybar_stan2$model <- "stan"
print(str(ybar_stan2))
print(str(ybar_svy_ht))
print(str(ybar_svy_gr))

  ybar_all <- rbind(ybar_stan2, ybar_svy_ht)
  ybar_all <- rbind(ybar_all, ybar_svy_gr)

##########################################
### Meta
##########################################
  ybar_all$metadata <- "J = 300, Nj ~ (100, 1000)"

##########################################
### Save
##########################################
  today <- Sys.Date()
  today <- gsub("-", "_", today)
  write.csv(ybar_all, file = paste(resdir, "/simulation_results_ybar_", today, ".csv", sep = ""))
  rm(ybar_stan, ybar_stan2, ybar_svy_ht, ybar_svy_gr, ybar_all)

##########################################
### Loop through parest files
##########################################
cat("--------------------------- parests_stan -------------------------\n")
  fil.list <- NULL
  for (i in 1:length(cluster.list)) {
    for (j in 1:length(num.units.list)) {
      for (k in 1:length(modelname.list)) {
        for (u in u.list) {
          fil.list <- c(fil.list,
                        list.files(resdir,
                                   paste("parests_stan_usesizes_", u, "_nclusters_", cluster.list[i],
                                         "_nunits_", num.units.list2[j],
                                         "_sim_.*_", modelname.list[k], ".txt", sep = "")))
        }
      }
    }
  }
  nres <- length(fil.list)
  parests_stan <- data.frame()
  for (f in 1:nres) {
    filnam <- fil.list[f]
    tt <- file.info(paste(resdir, "/", filnam, sep = ""))
    fdate <- tt$mtime
    fdate <- substr(fdate, 1, 10)
    fsize <- tt$size
    if (!(fdate %in% resdate) | (fsize == 0 & !is.na(fsize))) {
      if (fsize == 0 & !is.na(fsize)) {
        print(paste("********** Output file of size 0 for ", filnam, " ****************", sep = ""))
      }
      next
    }
    fil <- read.csv(filnam)
    fil$parname <- rownames(fil)
    rownames(fil) <- NULL
    ff <- unlist(strsplit(filnam, "[[:punct:]]"))
    nclusters <- as.numeric(ff[6])
    nunits <- ff[8]
    simno <- as.numeric(ff[10])
    fil$nclusters <- nclusters
    fil$nunits <- nunits
    fil$simno <- simno
    fil$use.sizes <- as.numeric(ff[4])
    fil$modelname <- ff[11]
    #fil$stan_covers_50 <- fil[, "X25."] <= fil$truth & fil$truth <= fil[, "X75."]
    #fil$stan_covers_95 <- fil[, "X2.5."] <= fil$truth & fil$truth <= fil[, "X97.5."]
    parests_stan <- rbind(parests_stan, fil)
  }
  
##########################################
### Loop through lmer files
##########################################
cat("--------------------------- parests_lmer -------------------------\n")
  fil.list <- NULL
  for (i in 1:length(cluster.list)) {
    for (j in 1:length(num.units.list)) {
      for (u in u.list) {
        fil.list <- c(fil.list,
                      list.files(resdir,
                                 paste("parests_lmer_usesizes_", u, "_nclusters_", cluster.list[i],
                                       "_nunits_", num.units.list2[j],
                                       "_sim_.*.txt", sep = "")))
      }
    }
  }
print("****************************************")
print("Lenght of fil.list:")
print(length(fil.list))
print("Lneght of unique fil.list:")
print(length(unique(fil.list)))
print(sort(fil.list))
  nres <- length(fil.list)
  parests_lmer <- data.frame()
  for (f in 1:nres) {
    filnam <- fil.list[f]
    tt <- file.info(paste(resdir, "/", filnam, sep = ""))
    fdate <- tt$mtime
    fdate <- substr(fdate, 1, 10)
    fsize <- tt$size
    if (!(fdate %in% resdate) | (fsize == 0 & !is.na(fsize))) {
      if (fsize == 0 & !is.na(fsize)) {
        print(paste("********** Output file of size 0 for ", filnam, " ****************", sep = ""))
      }
      next
    }
    fil <- read.csv(filnam)
    names(fil) <- c("lmer_est", "lmer_stderr", "truth", "whichmodel")
    # whichmodel = 1 -> HT; whichmodel = 2 -> GR
    #fil$parname <- rownames(fil)
    fil$parname <- c("gamma0", "alpha0", "gamma1", "alpha1",
                     "sigma_beta0", "sigma_beta1", "sigma_y",
                     "gamma0", "alpha0", "sigma_beta0", "sigma_beta1", "sigma_y")
    rownames(fil) <- NULL
    ff <- unlist(strsplit(filnam, "[[:punct:]]"))
    nclusters <- as.numeric(ff[6])
    nunits <- ff[8]
    simno <- as.numeric(ff[10])
    fil$nclusters <- nclusters
    fil$nunits <- nunits
    fil$simno <- simno
    fil$use.sizes <- as.numeric(ff[4])
    #fil$lmer_covers_50 <- (fil$lmer_est - fil$lmer_stderr) <= fil$truth &
    #                       fil$truth <= (fil$lmer_est + fil$lmer_stderr)
    #fil$lmer_covers_95 <- (fil$lmer_est - 2*fil$lmer_stderr) <= fil$truth &
    #                       fil$truth <= (fil$lmer_est + 2*fil$lmer_stderr)
    f1 <- fil[fil$whichmodel == 1, ]
    f1$whichmodel <- NULL
    names(f1)[names(f1) == "lmer_est"] <- "lmer_est_m1"
    names(f1)[names(f1) == "lmer_stderr"] <- "lmer_stderr_m1"
    f2 <- fil[fil$whichmodel == 2, c("parname", "lmer_est", "lmer_stderr")]
    f2$whichmodel <- NULL
    names(f2)[names(f2) == "lmer_est"] <- "lmer_est_m2"
    names(f2)[names(f2) == "lmer_stderr"] <- "lmer_stderr_m2"
    f3 <- merge(f1, f2, by = c("parname"), all.x = TRUE)
    parests_lmer <- rbind(parests_lmer, f3)
  }
  rm(fil)
print("done with parests_lmer")
  
##########################################
### Combine stan and lmer parests
##########################################
parests_stan <- parests_stan[with(parests_stan, order(parname, simno, nunits, nclusters, use.sizes)), ]
parests_lmer <- parests_lmer[with(parests_lmer, order(parname, simno, nunits, nclusters, use.sizes)), ]
print(str(parests_stan))
print(parests_stan[1:4, ])
print(str(parests_lmer))
print(parests_lmer[1:4, ])
print(max(parests_lmer$simno))
  parests_all <- merge(parests_stan, parests_lmer,
                       by = c("parname", "simno", "nunits", "nclusters", "use.sizes"), all.x = TRUE)
  # added all.x = TRUE to keep mu, phi, Mj_unsamp, etc 
print(str(parests_all))
print(parests_all[1:4, ])
  rm(parests_stan, parests_lmer)
  parests_all$truth[parests_all$parname == "ybar_hat"] <- ybar_true
  parests_all$truth[parests_all$use.sizes == 0 & parests_all$parname == "Mj_tot_est"] <- Mj_tot_true_0
  parests_all$truth[parests_all$use.sizes == 1 & parests_all$parname == "Mj_tot_est"] <- Mj_tot_true_1
  parests_all$stan_covers_50 <- parests_all[, "X25."] <= parests_all$truth &
                                parests_all$truth <= parests_all[, "X75."]
  parests_all$stan_covers_95 <- parests_all[, "X2.5."] <= parests_all$truth &
                                parests_all$truth <= parests_all[, "X97.5."]
  parests_all$lmer_m1_covers_50 <- (parests_all$lmer_est_m1 - parests_all$lmer_stderr_m1) <= parests_all$truth &
                                 parests_all$truth <= (parests_all$lmer_est_m1 + parests_all$lmer_stderr_m1)
  parests_all$lmer_m1_covers_95 <- (parests_all$lmer_est_m1 - 2*parests_all$lmer_stderr_m1) <= parests_all$truth &
                                 parests_all$truth <= (parests_all$lmer_est_m1 + 2*parests_all$lmer_stderr_m1)
  parests_all$lmer_m2_covers_50 <- (parests_all$lmer_est_m2 - parests_all$lmer_stderr_m2) <= parests_all$truth &
                                 parests_all$truth <= (parests_all$lmer_est_m2 + parests_all$lmer_stderr_m2)
  parests_all$lmer_m2_covers_95 <- (parests_all$lmer_est_m2 - 2*parests_all$lmer_stderr_m2) <= parests_all$truth &
                                 parests_all$truth <= (parests_all$lmer_est_m2 + 2*parests_all$lmer_stderr_m2)

##########################################
### Meta
##########################################
  parests_all$metadata <- "J = 300, Nj ~ (100, 1000)"

##########################################
### save
##########################################
  today <- Sys.Date()
  today <- gsub("-", "_", today)
  write.csv(parests_all, file = paste(resdir, "/simulation_results_parests_", today, ".csv", sep = ""))


