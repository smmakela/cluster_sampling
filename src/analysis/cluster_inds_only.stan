data {
  int<lower=0> J_pop;         // number of clusters in population
  int<lower=0> J_sam;         // number of clusters in sample
  int<lower=0> N_sam;         // sample size
  int cluster_id_long[N_sam]; // cluster ids for sample units
  vector[N_sam] xi_sam;       // individual-level covar ("age")
  vector[N_sam] yi_sam;       // outcomes
  vector[J_sam] Mj_sam;       // vector of cluster sizes for sampled clusters
  vector[J_sam] logMj_sam;    // log of cluster sizes for sampled clusters
  vector[J_pop] Mj_pop;       // cluster sizes
  vector[J_pop] xbar_pop;     // cluster avgs of xij
}
transformed data{
  vector[J_pop] logMj_pop;
  logMj_pop <- log(Mj_pop) - mean(log(Mj_pop));
}
parameters {
  real<lower=0> sigma_beta0;
  real<lower=0> sigma_beta1;
  real<lower=0> sigma_y;
  real gamma0;
  real alpha0;
  vector[J_sam] eta0;
  vector[J_sam] eta1;
}
transformed parameters {
  vector[J_sam] beta0jhat;
  vector[J_sam] beta1jhat;
  vector[N_sam] yhat;
  
  beta0jhat <- gamma0 + eta0 * sigma_beta0;
  beta1jhat <- alpha0 + eta1 * sigma_beta1;

  for (i in 1:N_sam) {
    yhat[i] <- beta0jhat[cluster_id_long[i]] + xi_sam[i] * beta1jhat[cluster_id_long[i]];
  }
}
model {
  sigma_beta0 ~ cauchy(0, 2.5);
  sigma_beta1 ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);
  gamma0 ~ normal(0, 1);
  alpha0 ~ normal(0, 1);
  eta0 ~ normal(0, 1);
  eta1 ~ normal(0, 1);
  yi_sam ~ normal(yhat, sigma_y);
}
generated quantities {
  real beta0jhat_draw;
  real beta1jhat_draw;
  vector[J_pop] y_j_hat;
  real epsilon;
  real ybar_hat;

  for (j in 1:J_pop) {
    if (j <= J_sam) {
      y_j_hat[j] <- beta0jhat[j] + beta1jhat[j] * xbar_pop[j];
    }
    else {
      epsilon <- normal_rng(0, 1);
      beta0jhat_draw <- gamma0 + sigma_beta0 * epsilon;
      epsilon <- normal_rng(0, 1);
      beta1jhat_draw <- alpha0 + sigma_beta1 * epsilon;
      y_j_hat[j] <- beta0jhat_draw + beta1jhat_draw * xbar_pop[j];
    }
  }
  ybar_hat <- sum(y_j_hat .* Mj_pop) / sum(Mj_pop);
}

