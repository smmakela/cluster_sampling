data {
  int<lower=0> J_pop;                       // number of clusters
  int<lower=0> J_sam;                       // number of clusters
  int<lower=0> N_pop;                       // population size
  int<lower=0> N_sam;                       // population size
  int cluster_id_long[N_sam];            // renumbered cluster ids, length = N_sam
  vector[J_pop] Mj_pop;                     // cluster sizes
  vector[J_sam] Mj_sam;                     // vector of cluster sizes for sampled clusters
  vector[J_pop] Mj_mis;                     // number of missing sizes for each clusters;
  vector[N_sam] xi_sam;                     // individual-level covar ("age")
  vector[N_sam] yi_sam;                     // outcomes
  vector[J_sam] ybar_sam;                   // sample means by cluster
  vector[J_pop] xbar_pop;                   // cluster avgs of xij for unsampled units
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
  real ybar_samp_mis;
  real beta0j_draw;
  real beta1j_draw;
  vector[J_pop] y_j_hat;
  real ybar_hat;
  //real epsilon0;
  real epsilon1;
  real epsilon2;
  //real epsilon3;

  for (j in 1:J_pop) {
    // if the jth cluster was sampled, then we can use beta0j, beta1j
    if (j <= J_sam) {
      //epsilon0 <- normal_rng(0,1);
      ybar_samp_mis <- beta0jhat[j] + beta1jhat[j]*xbar_pop[j];
      y_j_hat[j] <- (Mj_sam[j]*ybar_sam[j] + Mj_mis[j]*ybar_samp_mis)/Mj_pop[j];
    }
    else {
      epsilon1 <- normal_rng(0,1);
      epsilon2 <- normal_rng(0,1);
      beta0j_draw <- gamma0 + (epsilon1 * sigma_beta0);
      beta1j_draw <- alpha0 + (epsilon2 * sigma_beta1);
      //epsilon3 <- normal_rng(0,1);
      y_j_hat[j] <- beta0j_draw + (beta1j_draw .* xbar_pop[j]);
    }
  }
  
  ybar_hat <- sum(y_j_hat .* Mj_pop) / sum(Mj_pop);
}
