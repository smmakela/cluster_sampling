data {
  int<lower=0> J_pop;                       // number of clusters
  int<lower=0> J_sam;                       // number of clusters
  int<lower=0> N_sam;                       // population size
  int<lower=0> M;			    // number of unique cluster sizes
  int cluster_id_long[N_sam];               // renumbered cluster ids, length = N_sam
  vector[N_sam] xi_sam;                     // individual-level covar ("age")
  vector[N_sam] yi_sam;                     // outcomes
  vector[J_sam] Mj_sam;                     // vector of cluster sizes for sampled clusters
  vector[J_sam] logMj_sam;                  // log of cluster sizes
  vector[J_pop] xbar_pop;                   // cluster avgs of xij
  real Tx;                                  // total pop size (of exp(log(Mj)-mean(log(Mj))))
  int n[M];    			            // counts of sampled clusters
}
transformed data {
  int JminusK;
  vector[M] alpha_phi; // for prior on phi

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
  vector[J_sam] eta0;
  vector[J_sam] eta1;
  simplex[M] phi;
}
transformed parameters {
  vector[J_sam] beta0jhat;
  vector[J_sam] beta1jhat;
  vector[N_sam] yhat;
  vector[M] phi_star_unnorm;
  vector[M] phi_star;
  real cee; // normalizer for phi_star
  real pi;

  beta0jhat <- gamma0 + (gamma1 * logMj_sam) + eta0 * sigma_beta0;
  beta1jhat <- alpha0 + (alpha1 * logMj_sam) + eta1 * sigma_beta1;

  for (i in 1:N_sam) {
    yhat[i] <- beta0jhat[cluster_id_long[i]] + xi_sam[i] * beta1jhat[cluster_id_long[i]];
  }
  for (m in 1:M) {
    pi <- M*Mj_sam[m]/Tx;
    phi_star_unnorm[m] <- phi[m]*(1-pi)/pi;
  }
  cee <- sum(phi_star_unnorm); // calculate normalizing constant for phi_star
  phi_star <- phi_star_unnorm/cee; // normalize phi_star
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
  yi_sam ~ normal(yhat, sigma_y);
  phi ~ dirichlet(alpha_phi);
  n ~ multinomial(phi);
}
generated quantities {
  int n_star[M]; // counts of sizes for nonsampled clusters
  int count;
  int start;
  int stop;
  vector[JminusK] Mj_unsamp;
  vector[J_pop] Mj_star;
  real Mj_tot_est;
  vector[J_pop] logMj_star;
  vector[J_pop] y_j_hat;
  real beta0jhat_draw;
  real beta1jhat_draw;
  real epsilon1;
  real epsilon2;
  real ybar_hat;

  // draw counts of how many times each distinct cluster size is repeated
  n_star <- multinomial_rng(phi_star, JminusK);

  // populate the vector of unsampled sizes using the appropriate number of times from n_star
  count <- 1;
  for (m in 1:M) {
    start <- count;
    stop <- start + n_star[m] - 1; // so that stop - start has n_star[m] elements inclusive
    for (p in start:stop) {
      Mj_unsamp[p] <- Mj_sam[m];
    }
    count <- start + n_star[m];
  }
  Mj_star <- append_row(Mj_sam, Mj_unsamp);
  Mj_tot_est <- sum(Mj_star);
  logMj_star <- log(Mj_star) - mean(log(Mj_star));

  for (j in 1:J_pop) {
    // if the jth cluster was sampled, then we can use beta0j, beta1j
    if (j <= J_sam) {
      y_j_hat[j] <- beta0jhat[j] + beta1jhat[j] * xbar_pop[j];
    }
    else {
      epsilon1 <- normal_rng(0,1);
      epsilon2 <- normal_rng(0,1);
      beta0jhat_draw <- gamma0 + (gamma1 * logMj_star[j]) + (epsilon1 * sigma_beta0);
      beta1jhat_draw <- alpha0 + (alpha1 * logMj_star[j]) + (epsilon2 * sigma_beta1);
      y_j_hat[j] <- beta0jhat_draw + beta1jhat_draw * xbar_pop[j];
    }
  }
  ybar_hat <- sum(y_j_hat .* Mj_star) / sum(Mj_star);
}
