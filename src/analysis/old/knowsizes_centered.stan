data {
  int<lower=0> J_pop;         // number of clusters in population
  int<lower=0> J_sam;         // number of clusters in sample
  int<lower=0> N_sam;         // sample size
  int cluster_id_long[N_sam]; // cluster ids for sample units
  vector[N_sam] x;            // individual-level covar ("age")
  vector[N_sam] y;            // outcomes
  vector[J_sam] Mj_sam;       // vector of cluster sizes for sampled clusters
  vector[J_sam] logMj_sam;    // log of cluster sizes for sampled clusters
  vector[J_pop] Mj_pop;       // cluster sizes
  vector[J_pop] xbar_pop;     // cluster avgs of xij
}
transformed data{
  vector[J_pop] logMj_pop;
  int J_mis;

  logMj_pop = log(Mj_pop) - mean(log(Mj_pop));
  J_mis = J_pop - J_sam;
}
parameters {
  real<lower=0> sigma_beta0;
  real<lower=0> sigma_beta1;
  real<lower=0> sigma_y;
  real alpha0;
  real gamma0;
  real alpha1;
  real gamma1;
  vector[J_sam] beta0;
  vector[J_sam] beta1;
}
transformed parameters {
  vector[N_sam] ymean;

  for (i in 1:N_sam) {
    ymean[i] = beta0[cluster_id_long[i]] + beta1[cluster_id_long[i]]*x[i];
  }
}
model {
  sigma_beta0 ~ cauchy(0, 2.5);
  sigma_beta1 ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 1);
  gamma0 ~ normal(0, 1);
  alpha1 ~ normal(0, 1);
  gamma1 ~ normal(0, 1);
  beta0 ~ normal(alpha0 + gamma0 * logMj_sam, sigma_beta0);
  beta1 ~ normal(alpha1 + gamma1 * logMj_sam, sigma_beta1);
  y ~ normal(ymean, sigma_y);
}
generated quantities {
  vector[J_mis] beta0_new;
  vector[J_mis] beta1_new;
  vector[J_pop] y_new;
  real ybar_new;

  // for the sample clusters, use posterior means of beta0, beta1
  y_new[1:J_sam] = beta0 + beta1 .* xbar_pop[1:J_sam];

  // for unsampled clusters, need to first draw new beta0, beta1 from their posteriors
  for (j in 1:J_mis) {
    beta0_new[j] = normal_rng(alpha0 + gamma0 * logMj_pop[J_sam + j], sigma_beta0);
    beta1_new[j] = normal_rng(alpha1 + gamma1 * logMj_pop[J_sam + j], sigma_beta1);
  }
  y_new[(J_sam + 1):J_pop] = beta0_new + beta1_new .* xbar_pop[(J_sam + 1):J_pop];

  /* commented out for now -- delete if above works
  for (j in 1:J_pop) {
    if (j <= J_sam) {
      y_new[j] = beta0[j] + beta1[j] * xbar_pop[j];
    }
    else {
      beta0_new = normal_rng(alpha0 + gamma0 * logMj_pop[j], sigma_beta0);
      beta1_new = normal_rng(alpha1 + gamma1 * logMj_pop[j], sigma_beta1);
      y_new[j] = beta0_draw + beta1_draw * xbar_pop[j];
    }
  }
  */
  ybar_new = sum(y_new .* Mj_pop) / sum(Mj_pop);
}

