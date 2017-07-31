data {
  int<lower=0> J_pop;         // number of clusters in population
  int<lower=0> J_sam;         // number of clusters in sample
  int<lower=0> N_sam;         // sample size
  int cluster_id_long[N_sam]; // cluster ids for sample units
  vector[N_sam] x;       // individual-level covar ("age")
  vector[N_sam] y;       // outcomes
  vector[J_sam] Mj_sam;       // vector of cluster sizes for sampled clusters
  vector[J_sam] logMj_sam;    // log of cluster sizes for sampled clusters
  vector[J_pop] Mj_pop;       // cluster sizes
  vector[J_pop] xbar_pop;     // cluster avgs of xij
}
transformed data{
  vector[J_pop] logMj_pop;
  logMj_pop = log(Mj_pop) - mean(log(Mj_pop));
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
  for (i in 1:N_sam) {
    y[i] ~ normal(beta0[cluster_id_long[i]] + beta1[cluster_id_long[i]] * x[i], sigma_y);
  }
}
generated quantities {
  vector[J_pop] y_new;
  real ybar_new;

  y_new = beta0 + beta1 .* xbar_pop;
  ybar_new = sum(y_new .* Mj_pop) / sum(Mj_pop);
}

