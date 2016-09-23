data {
  int<lower=0> J_pop;                       // number of clusters
  int<lower=0> J_sam;                       // number of clusters
  int<lower=0> N_sam;                       // population size
  int cluster_id_long[N_sam];            // renumbered cluster ids, length = N_sam
  vector[N_sam] xi_sam;                     // individual-level covar ("age")
  vector[J_sam] logMj_sam;                  // log of cluster sizes
  vector[N_sam] yi_sam;                     // outcomes
  vector[J_sam] ybar_sam;                   // sample means by cluster
  vector[J_pop] xbar_pop;                   // cluster avgs of xij
  vector[J_pop] Mj_mis;                     // number of nonsampled units
  vector[J_sam] nj;                         // number of sampled units per sampled cluster
  int Mj_sam[J_sam];                        // vector of cluster sizes for sampled clusters
}
transformed data {
  int JminusK;

  JminusK <- J_pop - J_sam; // number of unsampled clusters
}
parameters {
  real<lower=0> sigma_beta0;
  real<lower=0> sigma_beta1;
  real<lower=0> sigma_y;
  real gamma0;
  real gamma1;
  real alpha0;
  real alpha1;
  vector[J_sam] beta0j;
  vector[J_sam] beta1j;
  vector[J_sam] eta0;
  vector[J_sam] eta1;
  real<lower=0> phi;
  real<lower=0> mu;
}
transformed parameters {
  vector[N_sam] yhat;
  real<lower=0> mu_star;
  real<lower=0> phi_star;
  real<lower=0> k;
  real<lower=0> p;
  
  k <- phi;
  p <- phi / (phi + mu);
  
  // size-biased params
  phi_star <- k + 1;
  mu_star <- (k + 1) * (1 - p) / p;

  for (i in 1:N_sam) {
    yhat[i] <- beta0j[cluster_id_long[i]] + xi_sam[i] * beta1j[cluster_id_long[i]];
  }
}
model {
  mu ~ normal(1000,5);
  phi ~ normal(1,.1);
  sigma_beta0 ~ cauchy(0, 2.5);
  sigma_beta1 ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);
  gamma0 ~ normal(0, 1);
  gamma1 ~ normal(0, 1);
  alpha0 ~ normal(0, 1);
  alpha1 ~ normal(0, 1);
  eta0 ~ normal(0, 1);
  eta1 ~ normal(0, 1);
  beta0j ~ normal(gamma0 + gamma1*logMj_sam, sigma_beta0);
  beta1j ~ normal(alpha0 + alpha1*logMj_sam, sigma_beta1);
  yi_sam ~ normal(yhat, sigma_y);
  Mj_sam ~ neg_binomial_2(mu_star, phi_star);
}
generated quantities {
  real beta0j_draw;
  real beta1j_draw;
  vector[J_pop] y_j_hat;
  real ybar_hat;
  real ybar_samp_mis;
  //real epsilon0;
  real epsilon1;
  real epsilon2;
  //real epsilon3;
  int Mj_unsamp[JminusK];
  int kk;
  int Mj_tot_est;
  vector[J_pop] Mj_all;
  vector[J_pop] logMj_all;

  for (j in 1:JminusK) {
    Mj_unsamp[j] <- neg_binomial_2_rng(mu, phi);
  }
  Mj_tot_est <- sum(Mj_sam) + sum(Mj_unsamp);
  Mj_all <- append_row(to_vector(Mj_sam), to_vector(Mj_unsamp));
  // if we drew a cluster size of zero, set it to 1 here so we can take logs
  for (j in 1:JminusK) {
    if (Mj_unsamp[j] == 0) {
      Mj_all[j + J_sam] <- 1.0;
    }
  }
  logMj_all <- log(Mj_all) - mean(log(Mj_all));

  kk <- 1;
  ybar_hat <- 0.0;
  for (j in 1:J_pop) {
    // if the jth cluster was sampled, then we can use beta0j, beta1j
    if (j <= J_sam) {
      //epsilon0 <- normal_rng(0,1);
      ybar_samp_mis <- beta0j[j] + beta1j[j]*xbar_pop[j];
      y_j_hat[j] <- nj[j]*ybar_sam[j] + Mj_mis[j]*ybar_samp_mis ;
    }
    else {
      epsilon1 <- normal_rng(0,1);
      epsilon2 <- normal_rng(0,1);
      //epsilon3 <- normal_rng(0,1);
      //beta0j_draw <- gamma0 + (gamma1 * log(Mj_unsamp[kk])) + (epsilon1 * sigma_beta0);
      //beta1j_draw <- alpha0 + (alpha1 * log(Mj_unsamp[kk])) + (epsilon2 * sigma_beta1);
      beta0j_draw <- gamma0 + (gamma1 * logMj_all[j]) + (epsilon1 * sigma_beta0);
      beta1j_draw <- alpha0 + (alpha1 * logMj_all[j]) + (epsilon2 * sigma_beta1);
      if (Mj_unsamp[kk] == 0) {
        y_j_hat[j] <- 0;
      } else {
        y_j_hat[j] <- Mj_unsamp[kk]*(beta0j_draw + (beta1j_draw .* xbar_pop[j]));
      }
      kk <- kk + 1;
    }
    ybar_hat <- ybar_hat + y_j_hat[j];
  }
  ybar_hat <- ybar_hat / Mj_tot_est;

}
