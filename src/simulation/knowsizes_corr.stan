data {
  int<lower=0> J_pop;         // number of clusters in population
  int<lower=0> J_sam;         // number of clusters in sample
  int<lower=0> N_sam;         // sample size
  int cluster_id_long[N_sam]; // cluster ids for sample units
  matrix[N_sam, 2] x;            // individual-level covar ("age")
  vector[N_sam] y;            // outcomes
  row_vector[2] u[J_sam];
}
parameters {
  real<lower=0> sigma_y;
  corr_matrix[2] Omega;
  vector<lower=0>[2] sigma_betas;
  matrix[2, 2] mu_betas;
  vector[2] betas[J_sam];
}
model {
  sigma_betas ~ cauchy(0, 2.5);
  Omega ~ lkj_corr(1);
  to_vector(mu_betas) ~ normal(0, 5);
  {
    row_vector[2] u_times_mu_beta[J_sam];
    for (j in 1:J_sam) {
      u_times_mu_beta[j] = u[j] * mu_betas;
    }
    betas ~ multi_normal(u_times_mu_beta, quad_form_diag(Omega, sigma_betas));
  }
  sigma_y ~ cauchy(0, 2.5);
  for (i in 1:N_sam) {
    y[i] ~ normal(x[i] * betas[cluster_id_long[i]], sigma_y);
  }
}
