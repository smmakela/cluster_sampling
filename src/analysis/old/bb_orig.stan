data {
  int<lower=0> J_pop;                       // number of clusters
  int<lower=0> J_sam;                       // number of clusters
  int<lower=0> N_sam;                       // population size
  int<lower=0> M;			    // number of unique cluster sizes
  int cluster_id_long[N_sam];               // renumbered cluster ids, length = N_sam
  vector[N_sam] x;                          // individual-level covar ("age")
  vector[N_sam] y;                          // outcomes
  vector[J_sam] Mj_sam;                     // vector of cluster sizes for sampled clusters
  vector[J_sam] logMj_sam;                  // log of cluster sizes
  vector[J_pop] xbar_pop;                   // cluster avgs of xij
  real Tx;                                  // total pop size (of exp(log(Mj)-mean(log(Mj))))
  int n[M];    			            // counts of sampled clusters
}
transformed data {
  int J_mis;
  vector[M] alpha_phi; // for prior on phi

  J_mis = J_pop - J_sam; // number of unsampled clusters
  alpha_phi = rep_vector(1, M); // create vector of M zeros
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
  vector[J_sam] beta0;
  vector[J_sam] beta1;
  vector[N_sam] yhat;
  vector[M] phi_star_unnorm;
  vector[M] phi_star;
  real cee; // normalizer for phi_star
  real pi;

  beta0 <- alpha0 + gamma0 * logMj_sam + eta0 * sigma_beta0;
  beta1 <- alpha1 + gamma1 * logMj_sam + eta1 * sigma_beta1;

  for (i in 1:N_sam) {
    yhat[i] = beta0[cluster_id_long[i]] + x[i] * beta1[cluster_id_long[i]];
  }
  for (m in 1:M) {
    pi = M*Mj_sam[m]/Tx;
    phi_star_unnorm[m] = phi[m]*(1-pi)/pi;
  }
  cee = sum(phi_star_unnorm); // calculate normalizing constant for phi_star
  phi_star = phi_star_unnorm/cee; // normalize phi_star
}
model {
  sigma_beta0 ~ cauchy(0, 2.5);
  sigma_beta1 ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 1);
  gamma0 ~ normal(0, 1);
  alpha1 ~ normal(0, 1);
  gamma1 ~ normal(0, 1);
  eta0 ~ normal(0, 1);
  eta1 ~ normal(0, 1);
  //beta0 ~ normal(alpha0 + gamma0 * logMj_sam, sigma_beta0);
  //beta1 ~ normal(alpha1 + gamma1 * logMj_sam, sigma_beta1);
  y ~ normal(yhat, sigma_y);
  phi ~ dirichlet(alpha_phi);
  n ~ multinomial(phi);
}
generated quantities {
  int n_star[M]; // counts of sizes for nonsampled clusters
  int count;
  int start;
  int stop;
  vector[J_mis] Mj_mis;
  vector[J_pop] Mj_pop_est;
  real M_tot_est;
  vector[J_pop] logMj_pop_est;
  vector[J_pop] y_new;
  vector[J_mis] beta0_new;
  vector[J_mis] beta1_new;
  real ybar_new;

  // draw counts of how many times each distinct cluster size is repeated
  n_star = multinomial_rng(phi_star, J_mis);

  // populate the vector of unsampled sizes using the appropriate number of times from n_star
  count = 1;
  for (m in 1:M) {
    start = count;
    stop = start + n_star[m] - 1; // so that stop - start has n_star[m] elements inclusive
    for (p in start:stop) {
      Mj_mis[p] = Mj_sam[m];
    }
    count = start + n_star[m];
  }
  Mj_pop_est = append_row(Mj_sam, Mj_mis);
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
