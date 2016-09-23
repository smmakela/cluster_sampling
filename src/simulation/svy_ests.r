# Author: Susanna Makela
# Date: Jan 28, 2016
# Purpose: use the survey package to estimate ybar from the sampled data
svy_ests <- function(J, num.clusters, num.units, rootdir, sim, popseed, use.sizes) {
  # Inputs:
  #   J -- number of clusters in pop
  #   num.clusters -- number of sampled clusters
  #   num.units -- number of sampled units
  #   rootdir -- directory where everything is stored
  #   sim -- which simulation this is
  #   popseed -- which population to use

  # load survey package
    library(survey)

  # load pop data -- use this to get the total sample size so we can calculate the sampling probs (here we are assuming we know M_j in all clusters, so obviously we can calculate the sampling probs
    load(paste(rootdir, "/Data/Simplify/vary_K/popdata_usesizes_", use.sizes, "_seed_", popseed, ".RData", sep = ""))
    pop.data <- popdata[["pop.data"]]
    sizetot <- sum(pop.data$Mj)
    ybar_true <- mean(pop.data$yi)

  # load sample data
    if (num.units <= 1) {
      nunits <- paste(100*num.units, "pct", sep = "")
    } else {
      nunits <- num.units
    }
    load(paste(rootdir, "/Data/Simplify/vary_K/sampledata_usesizes_", use.sizes, "_nclusters_", num.clusters,
               "_nunits_", nunits, "_sim_", sim, ".RData", sep = ""))
    sample.data$fpc <- J # fpc is the number of clusters in the pop
    sample.data$prob <- num.clusters*sample.data$Mj/sizetot # prob of selecting the cluster
    if (num.units > 1) {
      sample.data$prob2 <- num.units/sample.data$Mj # prob of selecting each unit within the cluster
    } else {
      sample.data$prob2 <- num.units
    }
    sample.data$wt <- 1/sample.data$prob

  # HT ESTIMATE
  # describe survey design
    #des <- svydesign(id = ~new.cluster.id, fpc = ~fpc, weights = ~wt, data = sample.data)
    des <- svydesign(id = ~cluster.id+unit.id, fpc = ~prob+prob2, data = sample.data, pps = "brewer")

  # estimate pop mean, pull out std err
    tt <- svymean(~yi, des)
    ybar_hat <- as.numeric(tt[1])
    ybar_se <- as.numeric(sqrt(attr(tt,"var")))

  # if we got an estimate that's way off from the truth, record the info here
    if (abs((ybar_hat - ybar_true)/ybar_true) > 2) {
      write.table(sample, file = paste(rootdir, "/Results/Simplify/vary_K/svydesign_usesizes_", use.sizes, "_nclusters_", num.clusters,
                             "_nunits_", nunits, "_sim_", sim, ".txt", sep = "")) 
      save(des, file = paste(rootdir, "/Results/Simplify/vary_K/svydesign_usesizes_", use.sizes, "_nclusters_", num.clusters,
                             "_nunits_", nunits, "_sim_", sim, ".txt", sep = "")) 
    }

  # GREG ESTIMATE
    ptot <- sum(pop.data$xi)
    pop.totals <- c(`(Intercept)`=nrow(sample.data), xi = ptot)
    if (num.units > 1) { # self-weighting
      des2 <- calibrate(des, formula = ~xi, pop.totals)
    } else { # make weights be the same by cluster
      des2 <- calibrate(des, formula = ~xi, pop.totals, aggregate = 1) 
    }
    tt2 <- svymean(~yi, des2)
    ybar_hat2 <- as.numeric(tt2[1])
    ybar_se2 <- as.numeric(sqrt(attr(tt2,"var")))

  # if we got an estimate that's way off from the truth, record the info here
    d1 <- abs((ybar_hat - ybar_true)/ybar_true)
    d2 <- abs((ybar_hat2 - ybar_true)/ybar_true)
    if ((d1 > 2) | (d2 > 2)) {
      write.csv(sample.data, file = paste(rootdir, "/Results/Simplify/vary_K/sampdat_weird_svy_usesizes_", use.sizes,
                                          "_nclusters_", num.clusters, "_nunits_", nunits, "_sim_", sim, ".txt", sep = ""))
                              
    }

  # print
    print("**************************************************")
    print(paste("ybar true: ", round(ybar_true, digits=2), sep = ""))
    print(paste("HT estimate (se): ", round(ybar_hat, digits=2), " (", round(ybar_se, digits=2), ")", sep = ""))
    print(paste("GREG estimate (se): ", round(ybar_hat2, digits=2), " (", round(ybar_se2, digits=2), ")", sep = ""))
    print("**************************************************")

  # save results
    res <- data.frame(sim, ybar_true, ybar_hat, ybar_se, ybar_hat2, ybar_se2)
    write.table(res, file = paste(rootdir, "/Results/Simplify/vary_K/ybar_svy_usesizes_", use.sizes, 
                                  "_nclusters_", num.clusters, "_nunits_", nunits, "_sim_", sim, ".txt", sep = ""), sep = ",")

}
