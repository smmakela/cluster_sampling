data {
  int<lower=0> J_pop;         // number of clusters in population
  int<lower=0> J_sam;         // number of clusters in sample
  int<lower=0> N_sam;         // sample size
  int cluster_id_long[N_sam]; // cluster ids for sample units
  vector[N_sam] x;            // individual-level covar ("age")
  vector[N_sam] y;            // outcomes
  vector[J_pop] Mj_pop;       // cluster sizes
  vector[J_pop] xbar_pop;     // cluster avgs of xij
}
transformed data {
  int<lower=0> J_mis;
  J_mis = J_pop - J_sam;
}
parameters {
  real<lower=0> sigma_beta0;
  real<lower=0> sigma_beta1;
  real<lower=0> sigma_y;
  real alpha0;
  real alpha1;
  vector[J_sam] beta0;
  vector[J_sam] beta1;
}
model {
  sigma_beta0 ~ cauchy(0, 2.5);
  sigma_beta1 ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 1);
  alpha1 ~ normal(0, 1);
  beta0 ~ normal(alpha0, sigma_beta0);
  beta1 ~ normal(alpha1, sigma_beta1);
  for (i in 1:N_sam) {
    y[i] ~ normal(beta0[cluster_id_long[i]] + beta1[cluster_id_long[i]]*x[i], sigma_y);
  }
}
generated quantities {
  real beta0_new[J_mis];
  real beta1_new[J_mis];
  vector[J_pop] y_new;
  real ybar_new;

  // for the sample clusters, use posterior means of beta0, beta1
  y_new[1:J_sam] = beta0 + beta1 .* xbar_pop[1:J_sam];

  // for unsampled clusters, need to first draw new beta0, beta1 from their posteriors
  beta0_new = normal_rng(alpha0, sigma_beta0);
  beta1_new = normal_rng(alpha1, sigma_beta1);
  y_new[(J_sam + 1):J_pop] = beta0_new + beta1_new * xbar_pop[(J_sam + 1):J_pop];

  ybar_new = sum(y_new .* Mj_pop) / sum(Mj_pop);
}

