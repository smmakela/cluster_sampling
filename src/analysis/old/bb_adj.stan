data {
  int<lower=0> J_pop;                       // number of clusters
  int<lower=0> J_sam;                       // number of clusters
  int<lower=0> N_sam;                       // population size
  int<lower=0> M;			    // number of unique cluster sizes
  int new_cluster_id_rep[N_sam];            // renumbered cluster ids, length = N_sam
  vector[N_sam] xi_sam;                     // individual-level covar ("age")
  vector[J_sam] logMj_sam;                  // log of cluster sizes
  vector[N_sam] yi_sam;                     // outcomes
  vector[J_sam] ybar_sam;                   // sample means by cluster
  vector[J_pop] xbar_pop;                   // cluster avgs of xij
  vector[J_pop] Mj_mis;                     // number of nonsampled units
  vector[J_sam] Mj_sam;                     // vector of cluster sizes for sampled clusters
  vector[J_sam] nj;                         // vector of sample sizes for sampled clusters
  real Tx;                                  // total pop size (of exp(log(Mj)-mean(log(Mj))))
  int n[M];    			            // counts of sampled clusters
}
transformed data {
  //vector[J_sam] Mj_sam;
  int JminusK;
  vector[M] alpha_phi; // for prior on phi

  //Mj_sam <- exp(logMj_sam);
  JminusK <- J_pop - J_sam; // number of unsampled clusters
  alpha_phi <- rep_vector(1, M); // create vector of M zeros
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
  simplex[M] phi;
}
transformed parameters {
  vector[N_sam] yhat;
  vector[M] phi_star;
  real cee; // normalizer for phi_star
  real pi;

  for (i in 1:N_sam) {
    yhat[i] <- beta0j[new_cluster_id_rep[i]] + xi_sam[i] * beta1j[new_cluster_id_rep[i]];
  }
  for (m in 1:M) {
    pi <- M*Mj_sam[m]/Tx;
    phi_star[m] <- phi[m]*(1-pi)/pi;
  }
  cee <- sum(phi_star); // calculate normalizing constant for phi_star
  phi_star <- phi_star/cee; // normalize phi_star
}
model {
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
  phi ~ dirichlet(alpha_phi);
  n ~ multinomial(phi);
}
generated quantities {
  int count;
  int start;
  int stop;
  vector[JminusK] Mj_unsamp;
  real beta0j_draw;
  real beta1j_draw;
  vector[J_pop] y_j_hat;
  real ybar_hat;
  real ybar_samp_mis;
  real epsilon0;
  real epsilon1;
  real epsilon2;
  real epsilon3;
  real logMj_star;
  int n_star[M]; // counts of sizes for nonsampled clusters
  vector[J_pop] Mj_all;
  int kk;
  real Mj_tot_est;

  n_star <- multinomial_rng(phi_star, JminusK);

  // populate the vector of unsampled sizes using the appropriate number of times from n_star
  count <- 1;
  for (m in 1:M) {
    start <- count;
    stop <- start + n_star[m] - 1; // so that stop - start has n_star[m] elements inclusive
    for (k in start:stop) {
      Mj_unsamp[k] <- Mj_sam[m]; // have to do this manually since stan doesn't allow segments on the lhs
    }
    count <- start + n_star[m];
  }
  Mj_all <- append_row(Mj_sam, Mj_unsamp);
  Mj_tot_est <- sum(Mj_all);

  kk <- 1;
  for (j in 1:J_pop) {
    // if the jth cluster was sampled, then we can use beta0j, beta1j
    if (j <= J_sam) {
      epsilon0 <- normal_rng(0,1);
      ybar_samp_mis <- beta0j[j] + beta1j[j]*xbar_pop[j] + epsilon0*sigma_y/Mj_mis[j];
      y_j_hat[j] <- (Mj_sam[j]*ybar_sam[j] + Mj_mis[j]*ybar_samp_mis)/Mj_all[j];
    }
    else {
      epsilon1 <- normal_rng(0,1);
      epsilon2 <- normal_rng(0,1);
      epsilon3 <- normal_rng(0,1);
      logMj_star <- log(Mj_unsamp[kk]);
      kk <- kk + 1;
      beta0j_draw <- gamma0 + (gamma1 * logMj_star) + (epsilon1 * sigma_beta0);
      beta1j_draw <- alpha0 + (alpha1 * logMj_star) + (epsilon2 * sigma_beta1);
      y_j_hat[j] <- beta0j_draw + (beta1j_draw .* xbar_pop[j]) + ((epsilon3 * sigma_y) / (Mj_all[j]));
    }
  }
  ybar_hat <- sum(y_j_hat .* Mj_all) / sum(Mj_all);
}
