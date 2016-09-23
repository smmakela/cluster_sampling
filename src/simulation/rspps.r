rspps <- function(sizes, ids, n) {
  # Purpose: do random systematic pps sampling
  # Inputs:
  #   sizes -- vector of PSU sizes
  #   ids -- vector of corresponding PSU ids
  #   n -- target sample size
  # Outputs:
  #   sel.ids -- vector of selected PSU ids
  
  # input checks
  if (length(sizes) != length(ids)) {
    stop("sizes and ids must be the same length")
  }
  if (length(unique(ids)) != length(ids)) {
    stop("ids must not contain duplicates")
  }
  
  # randomly order the sizes by permuting the ids
  set.seed(Sys.time())
  ids.orig <- ids
  sizes.orig <- sizes
  ids <- sample(ids.orig, length(ids.orig), replace = FALSE)
  sizes <- sizes.orig[ids]
  
  # make the z vector
  tot <- sum(sizes)
  z <- sizes/tot
  if (sum(n*z > 1) != 0) {
    stop("xi/X < 1/n must hold for all i!")
  }
  
  # make A and (u+k)
  A <- c(0, cumsum(n*z))
  u <- runif(1)
  uk <- u + c(0:(n-1))
  
  # figure out which indices we're sampling
  inds <- rep(NA, times = n)
  for (i in 1:n) {
    l <- min(which(A >= uk[i]))
    u <- max(which(A < uk[i])) + 1
    if (l == u) {
      inds[i] <- l-1 # subtract one because indexing into A conceptually starts at 0
    }
  }
  
  # to get the cluster *IDS*, we need to index into ids with the selected *INDICES*
  sel.ids <- ids[inds]
  
  return(sel.ids)
}
