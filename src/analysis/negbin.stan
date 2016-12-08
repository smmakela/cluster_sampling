data {
  int<lower=0> J_pop;                       // number of clusters
  int<lower=0> J_sam;                       // number of clusters
  int<lower=0> N_sam;                       // population size
  int cluster_id_long[N_sam];            // renumbered cluster ids, length = N_sam
  vector[N_sam] xi_sam;                     // individual-level covar ("age")
  vector[N_sam] yi_sam;                     // outcomes
  int Mj_sam[J_sam];                        // vector of cluster sizes for sampled clusters
  vector[J_sam] logMj_sam;                  // log of cluster sizes
  vector[J_pop] xbar_pop;                   // cluster avgs of xij
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
  vector[J_sam] eta0;
  vector[J_sam] eta1;
  real<lower=0> phi;
  real<lower=0> mu;
}
transformed parameters {
  vector[J_sam] beta0jhat;
  vector[J_sam] beta1jhat;
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

  beta0jhat <- gamma0 + gamma1*logMj_sam + eta0*sigma_beta0;
  beta1jhat <- alpha0 + alpha1*logMj_sam + eta1*sigma_beta1;
  for (i in 1:N_sam) {
    yhat[i] <- beta0jhat[cluster_id_long[i]] + xi_sam[i] * beta1jhat[cluster_id_long[i]];
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
  yi_sam ~ normal(yhat, sigma_y);
  Mj_sam ~ neg_binomial_2(mu_star, phi_star);
}
generated quantities {
  int<lower=0> Mj_mis[JminusK];
  vector[J_pop] Mj_pop_est;
  vector[J_pop] logMj_pop_est;
  int ll;
  real beta0jhat_draw;
  real beta1jhat_draw;
  vector[J_pop] y_j_hat;
  real epsilon;
  real ybar_hat;

  for (j in 1:JminusK) {
    ll <- J_sam + j;
    Mj_mis[ll] <- neg_binomial_2_rng(mu, phi);
    while (Mj_mis[ll] == 0) {
      Mj_mis[ll] <- neg_binomial_2_rng(mu, phi);
    }
  }
  Mj_pop_est <- append_row(to_vector(Mj_sam), to_vector(Mj_mis));
  logMj_pop_est <- log(Mj_pop_est) - mean(log(Mj_pop_est));

  for (j in 1:J_pop) {
    if (j <= J_sam) {
      y_j_hat[j] <- beta0jhat[j] + beta1jhat[j] * xbar_pop[j];
    }
    else {
      epsilon <- normal_rng(0, 1);
      beta0jhat_draw <- gamma0 + gamma1 * logMj_pop_est[j] + sigma_beta0 * epsilon;
      epsilon <- normal_rng(0, 1);
      beta1jhat_draw <- alpha0 + alpha1 * logMj_pop_est[j] + sigma_beta1 * epsilon;
      y_j_hat[j] <- beta0jhat_draw + beta1jhat_draw * xbar_pop[j];
    }
  }
  ybar_hat <- sum(y_j_hat .* Mj_pop_est) / sum(Mj_pop_est);
}

