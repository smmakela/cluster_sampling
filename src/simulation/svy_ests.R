# Author: Susanna Makela
# Date: Jan 28, 2016
# Purpose: use the survey package to estimate ybar from the sampled data
svy_ests <- function(J, num.clusters, num.units, use.sizes, outcome.type, rootdir, sim) {
  # Inputs:
  #   J -- number of clusters in pop
  #   num.clusters -- number of sampled clusters
  #   num.units -- number of sampled units
  #   use.sizes -- whether cluster sizes are related to outcome values (0/1)
  #   outcome.type -- whether outcomes are continuous or binary
  #   rootdir -- directory where everything is stored
  #   sim -- which simulation this is

  # load pop data -- use this to get the total sample size so we can calculate the sampling probs (here we are assuming we know M_j in all clusters, so obviously we can calculate the sampling probs
    if (num.units <= 1) {
      nunits <- paste(num.units*100, "pct", sep = "")
    } else {
      nunits <- num.units
    }
    popdata <- readRDS(paste0(rootdir, "/output/simulation/popdata_usesizes_",
                       use.sizes, "_", outcome.type, ".rds"))
    pop.data <- popdata[["pop.data"]]
    truepars <- popdata[["truepars"]]
    print(truepars)
    rm(popdata)

    sizetot <- sum(pop.data$Mj)
    ybar.true <- mean(pop.data$y)

    simdata <- readRDS(paste0(rootdir, "output/simulation/simdata_usesizes_",
                              use.sizes, "_", outcome.type, "_nclusters_",
                              num.clusters, "_nunits_", nunits, "_simno_", sim,
                              ".rds"))
print(str(simdata))
print(names(simdata))
    for (j in names(simdata)) {
      assign(j, simdata[[j]])
    }
    sample.data$fpc <- J # fpc is the number of clusters in the pop
    sample.data$prob <- num.clusters*sample.data$Mj/sizetot # prob of selecting the cluster
    if (num.units > 1) {
      sample.data$prob2 <- num.units/sample.data$Mj # prob of selecting each unit within the cluster
    } else {
      sample.data$prob2 <- num.units
    }
    sample.data$wt <- 1/sample.data$prob
    rm(simdata)


  # HT ESTIMATE
  # describe survey design
    #des <- svydesign(id = ~new.cluster.id, fpc = ~fpc, weights = ~wt, data = sample.data)
    des <- svydesign(id = ~cluster.id+unit.id, fpc = ~prob+prob2,
                     data = sample.data, pps = "brewer")

  # estimate pop mean, pull out std err
    tt <- svymean(~y, des)
    ybar.hat.hajek <- as.numeric(tt[1])
    ybar.se.hajek <- as.numeric(sqrt(attr(tt,"var")))

  # if we got an estimate that's way off from the truth, record the info here
    #if (abs((ybar_hat - ybar_true)/ybar_true) > 2) {
    #  write.table(sample, file = paste(rootdir, "/Results/Simplify/vary_K/svydesign_usesizes_", use.sizes, "_nclusters_", num.clusters,
    #                         "_nunits_", nunits, "_sim_", sim, ".txt", sep = "")) 
    #  save(des, file = paste(rootdir, "/Results/Simplify/vary_K/svydesign_usesizes_", use.sizes, "_nclusters_", num.clusters,
    #                         "_nunits_", nunits, "_sim_", sim, ".txt", sep = "")) 
    #}

  # GREG ESTIMATE
    ptot <- sum(pop.data$x)
    pop.totals <- c(`(Intercept)`=nrow(sample.data), x = ptot)
    if (num.units > 1) { # self-weighting
      des2 <- calibrate(des, formula = ~x, pop.totals)
    } else { # make weights be the same by cluster
      des2 <- calibrate(des, formula = ~x, pop.totals, aggregate = 1) 
    }
    tt2 <- svymean(~y, des2)
    ybar.hat.greg <- as.numeric(tt2[1])
    ybar.se.greg <- as.numeric(sqrt(attr(tt2,"var")))

  # if we got an estimate that's way off from the truth, record the info here
    #d1 <- abs((ybar_hat - ybar_true)/ybar_true)
    #d2 <- abs((ybar_hat2 - ybar_true)/ybar_true)
    #if ((d1 > 2) | (d2 > 2)) {
    #  write.csv(sample.data, file = paste(rootdir, "/Results/Simplify/vary_K/sampdat_weird_svy_usesizes_", use.sizes,
    #                                      "_nclusters_", num.clusters, "_nunits_", nunits, "_sim_", sim, ".txt", sep = ""))
    #}

  # print
    print("**************************************************")
    print(paste0("ybar true: ", round(ybar.true, digits = 2)))
    print(paste0("HT estimate (se): ", round(ybar.hat.hajek, digits = 2),
                 " (", round(ybar.se.hajek, digits = 2), ")"))
    print(paste0("GREG estimate (se): ", round(ybar.hat.greg, digits = 2),
                 " (", round(ybar.se.greg, digits = 2), ")"))
    print("**************************************************")

  # store results so that statistics are wide, model names are long
    res <- data.frame(ybar.true, ybar.hat.hajek, ybar.se.hajek,
                      ybar.hat.greg, ybar.se.greg)
    res %>%
      tidyr::gather(key = tmpname, value = value, -ybar.true) %>%
      tidyr::extract(col = tmpname, into = c("stat.name", "model.name"),
                    regex = "(^[^.]+[.][^.]+)(.+$)") %>%
      tidyr::spread(stat.name, value) -> res
    res$model.name <- gsub("\\.", "", res$model.name)
    res$ybar.hat.lci50 <- res$ybar.hat - res$ybar.se
    res$ybar.hat.uci50 <- res$ybar.hat + res$ybar.se
    res$ybar.hat.lci95 <- res$ybar.hat - 1.96*res$ybar.se
    res$ybar.hat.uci95 <- res$ybar.hat + 1.96*res$ybar.se
    res$ybar.se <- NULL
    return(res)
}
