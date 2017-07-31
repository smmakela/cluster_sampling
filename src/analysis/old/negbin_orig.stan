data {
  int<lower=0> J_pop;                // number of clusters
  int<lower=0> J_sam;                // number of clusters
  int<lower=0> N_sam;                // population size
  int cluster_id_long[N_sam];        // renumbered cluster ids, length = N_sam
  vector[N_sam] x;                   // individual-level covar ("age")
  vector[N_sam] y;                   // outcomes
  int Mj_sam[J_sam];                 // vector of sizes for sampled clusters
  vector[J_sam] logMj_sam;           // log of cluster sizes
  vector[J_pop] xbar_pop;            // cluster avgs of xij
}
transformed data {
  int J_mis;
  int Mj_sam_minus1[J_sam]; // vector of observed cluster sizes minus 1

  J_mis = J_pop - J_sam; // number of unsampled clusters
  
  for (j in 1:J_sam) {
    Mj_sam_minus1[j] = Mj_sam[j] - 1;
  }
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
  vector[J_sam] beta0;
  vector[J_sam] beta1;
  vector[N_sam] yhat;
  real<lower=0> mu_star;
  real<lower=0> phi_star;
  
  beta0 = alpha0 + gamma0 * logMj_sam + eta0 * sigma_beta0;
  beta1 = alpha1 + gamma1 * logMj_sam + eta1 * sigma_beta1;

  // size-biased params
  phi_star = phi + 1;
  mu_star = mu + (mu / phi);

  for (i in 1:N_sam) {
    yhat[i] = beta0[cluster_id_long[i]] + x[i] * beta1[cluster_id_long[i]];
  }
}
model {
  //mu ~ normal(1000,5);
  //phi ~ normal(1,.1);
  sigma_beta0 ~ cauchy(0, 2.5);
  sigma_beta1 ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 1);
  gamma0 ~ normal(0, 1);
  alpha1 ~ normal(0, 1);
  gamma1 ~ normal(0, 1);
  eta0 ~ normal(0, 1);
  eta1 ~ normal(0, 1);
  y ~ normal(yhat, sigma_y);
  Mj_sam_minus1 ~ neg_binomial_2(mu_star, phi_star);
}
generated quantities {
  int<lower=0> Mj_mis[J_mis];
  vector[J_pop] Mj_pop_est;
  vector[J_pop] logMj_pop_est;
  real M_tot_est;
  vector[J_mis] beta0_new;
  vector[J_mis] beta1_new;
  vector[J_pop] y_new;
  real ybar_new;

  for (j in 1:J_mis) {
    Mj_mis[j] = neg_binomial_2_rng(mu, phi);
    while (Mj_mis[j] == 0) {
      Mj_mis[j] = neg_binomial_2_rng(mu, phi);
    }
  }
  Mj_pop_est = append_row(to_vector(Mj_sam), to_vector(Mj_mis));
  M_tot_est = sum(Mj_pop_est);
  logMj_pop_est = log(Mj_pop_est) - mean(log(Mj_pop_est));

  // for the sample clusters, use posterior means of beta0, beta1
  y_new[1:J_sam] = beta0 + beta1 .* xbar_pop[1:J_sam];

  // for unsampled clusters, need to first draw new beta0, beta1 from their posteriors
  for (j in 1:J_mis) {
    beta0_new[j] = normal_rng(alpha0 + gamma0 * logMj_pop_est[J_sam + j], sigma_beta0);
    beta1_new[j] = normal_rng(alpha1 + gamma1 * logMj_pop_est[J_sam + j], sigma_beta1);
  }
  y_new[(J_sam + 1):J_pop] = beta0_new + beta1_new .* xbar_pop[(J_sam + 1):J_pop];

  ybar_new = sum(y_new .* Mj_pop_est) / sum(Mj_pop_est);
}

