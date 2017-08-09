#process_fit <- function(fit, parlist, betapars, kappapars = NA, sbpars = NA,
#                        stanmod_name) {

  # fit -- output of sampling()
  # parlist -- vector of regression parameter names
  # betapars -- vector of beta parameter names
  # kappapars -- vector of kappa parameter names (only used for strat models)
  # sbpars -- vector of size biased parameter names (only used for lognormal and negbin)
  # stanmod_name -- string for name of stan model so we know which parts of the code to run

  ##########################################
  ### Extract samples
  ##########################################
  param_samps <- data.frame(rstan::extract(fit, pars = parlist))

  # since the betas have to be passed in as a vector, deal with them separately
  beta_samps_orig <- data.frame(rstan::extract(fit, pars = betapars))
  nc <- nchar("beta0")
  beta_samps <- beta_samps_orig %>%
      tidyr::gather(key = pname, value = value) %>%
      dplyr::mutate(par = substr(pname, 1, nc),
                    cluster_id = as.numeric(substr(pname, nc+2, nchar(pname))),
                    pname = NULL) %>%
      dplyr::group_by(par, cluster_id) %>%
      dplyr::mutate(draw_num = row_number()) %>%
      tidyr::spread(key = par, value = value)

  # same for the kappas and sb pars, if we're doing ff
  #if (size_model == "ffstrat") {
  if (grepl("strat", stanmod_name)) {
    kappa_samps_orig <- data.frame(rstan::extract(fit, pars = kappapars))
    # this is for when we only have two strata
    kappa_samps <- kappa_samps_orig %>%
      dplyr::mutate(draw_num = row_number(),
                    stratum_id = 2)
    kappa_tmp <- kappa_samps
    kappa_tmp$stratum_id <- 1
    kappa_tmp$kappa0 <- 0.0
    kappa_tmp$kappa1 <- 0.0
    kappa_samps <- rbind(kappa_samps, kappa_tmp)
    #nck <- nchar("kappa0")
    #kappa_samps <- kappa_samps_orig %>%
    #    tidyr::gather(key = pname, value = value) %>%
    #    dplyr::mutate(par = substr(pname, 1, nck),
    #                  cluster_id = as.numeric(substr(pname, nck+2, nchar(pname))),
    #                  pname = NULL) %>%
    #    dplyr::group_by(par, cluster_id) %>%
    #    dplyr::mutate(draw_num = row_number()) %>%
    #    tidyr::spread(key = par, value = value) 
    if (grepl("lognormal", stanmod_name) || grepl("negbin", stanmod_name)) {
      sb_samps_orig <- data.frame(rstan::extract(fit, pars = sbpars))
      sb_samps <- sb_samps_orig %>%
          dplyr::mutate(draw_num = row_number()) %>%
          tidyr::gather(key = pname, value = value, -draw_num) %>%
          tidyr::separate(col = pname, into = c("par", "stratum_id"), sep = "\\.") %>%
          dplyr::mutate(stratum_id = as.numeric(stratum_id)) 
      sb_samps <- sb_samps %>%
          tidyr::spread(key = par, value = value)
    }
  }
  # same for phi_star, if we're doing bb
  if (grepl("bb", stanmod_name)) {
    nc2 <- nchar("phi_star")
    phi_star_samps_orig <- data.frame(rstan::extract(fit, pars = "phi_star"))
    print(str(phi_star_samps_orig))
    phi_star_samps_orig %>%
      tidyr::gather(key = pname, value = value) %>%
      dplyr::mutate(par = substr(pname, 1, nc2),
                    cluster_id = as.numeric(substr(pname, nc2+2, nchar(pname))),
                    pname = NULL) %>%
      dplyr::group_by(par, cluster_id) %>%
      dplyr::mutate(draw_num = row_number()) %>%
      tidyr::spread(key = par, value = value) -> phi_star_samps
    print(str(phi_star_samps))
  }

  print("done making samps")
  print(Sys.time())

  # Now make par_ests, the summary of the important parameters
  tt <- summary(fit)$summary
  if (!grepl("strat", stanmod_name)) {
    par_ests <- tt[parlist, ]
  } else {
    ii <- which(grepl(sbpars[1], rownames(tt)) | grepl(sbpars[2], rownames(tt)))
    jj <- which(grepl(kappapars[1], rownames(tt)) | grepl(kappapars[2], rownames(tt)))
    par_ests <- tt[c(parlist, rownames(tt)[c(ii, jj)]), ]
  }
  rm(fit)
  print("done making par_ests")
  print(Sys.time())

  print("making par_ests")
  print(Sys.time())
  par_ests_rownames <- attr(par_ests, "dimnames")[[1]]
  par_ests_rownames <- gsub("\\[", "", par_ests_rownames) 
  par_ests_rownames <- gsub("\\]", "", par_ests_rownames) 
  par_ests_colnames <- attr(par_ests, "dimnames")[[2]]
  par_ests <- data.frame(par_ests, row_names = par_ests_rownames)
  colnames(par_ests) <- par_ests_colnames

  print("printing par_ests")
  print(Sys.time())      
  print(str(par_ests))

#  to_return <- list(par_ests = par_ests, beta_samps = beta_samps)
#  if (grepl("bb", stanmod_name)) {
#    to_return <- c(to_return, list(phi_star_samps = phi_star_samps))
#  }
#  if (grepl("strat", stanmod_name)) {
#    to_return <- c(to_return, list(kappa_samps = kappa_samps))
#    if (grepl("lognormal", stanmod_name) || grepl("negbin", stanmod_name)) {
#      to_return <- c(to_return, list(sb_samps = sb_samps))
#    }
#  }
#  return(to_return)
#
#} # end function

