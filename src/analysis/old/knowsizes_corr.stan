data {
  int<lower=0> J_pop;         // number of clusters in population
  int<lower=0> J_sam;         // number of clusters in sample
  int<lower=0> N_sam;         // sample size
  int cluster_id_long[N_sam]; // cluster ids for sample units
  vector[N_sam] x;            // individual-level covar ("age")
  vector[N_sam] y;            // outcomes
  vector[J_sam] u;            // matrix with cluster-level predictors
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
  real<lower=0> sigma_y;
  corr_matrix[2] Omega;
  vector<lower=0>[2] sigma_betas;
  matrix[2, 2] mu_betas;
  vector[2] betas[J_sam];
}
transformed parameters {
  vector[N_sam] ymean;
  cov_matrix[2] Sigma_Betas;

  Sigma_Betas = quad_form_diag(Omega, sigma_betas);

  for (i in 1:N_sam) {
    ymean[i] = betai[cluster_id_long[i]] + beta1[cluster_id_long[i]]*x[i];
  }
}
model {
  sigma_betas ~ cauchy(0, 2.5);
  Omega ~ ljk_corr(1);
  {
    row_vector[2] u_times_mu_betas[J_sam];
    for (j in 1:J_sam) {
      u_times_mu_betas[j] = u[j] * mu_betas;
    }
    betas ~ multi_normal(u_times_mu_beta, Sigma_Betas);
  }
  y ~ normal(ymean, sigma_y);
}
/*
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
*/
