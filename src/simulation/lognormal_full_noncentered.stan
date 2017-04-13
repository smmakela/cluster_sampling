functions {
  real ybar_new_ln_rng(int J, int K, vector beta0, vector beta1, vector Nj_sample,
                       vector xbar_pop, real mu_star, real sigma,
                       real alpha0, real gamma0, real sigma_beta0,
                       real alpha1, real gamma1, real sigma_beta1, real sigma_y) {
    // here we get an estimate of y_bar, the finite-pop mean
    // for sampled clusters, we can use the estimate
    // for unsampled clusters, need to draw beta0, beta1, Nj, and then calculate
    // y_bar
    vector[J] beta0_new;
    vector[J] beta1_new;
    vector[J] Nj_new;
    vector[J] log_Nj_new;
    vector[J] yj_new;
    real ybar_new;
    
    beta0_new[1:K] = beta0;
    beta1_new[1:K] = beta1;
    Nj_new[1:K] = Nj_sample;
    yj_new[1:K] = beta0 + beta1 .* xbar_pop[1:K];
    
    for (k in (K+1):J) {
      Nj_new[k] = lognormal_rng(mu_star, sigma);
    }
    log_Nj_new = log(Nj_new) - mean(log(Nj_new));
    for (j in (K+1):J) {
      beta0_new[j] = normal_rng(alpha0 + gamma0 * log_Nj_new[j], sigma_beta0);
      beta1_new[j] = normal_rng(alpha1 + gamma1 * log_Nj_new[j], sigma_beta1);
      yj_new[j]  = normal_rng(beta0_new[j] + beta1_new[j] * xbar_pop[j],
                              sigma_y/sqrt(Nj_new[j]));
    }
    
    ybar_new = sum(yj_new .* to_vector(Nj_new)) / sum(Nj_new);
    
    return ybar_new;
    
  } # end est_ybar_new function
} # end functions block
data {
  int J; // number of pop clusters
  int K; // number of sampled clusters
  int n; // total sample size
  vector[n] x;
  vector[n] y;
  int cluster_id[n]; // vector of cluster id's for each sampled unit
  vector[K] Nj_sample;
  vector[K] log_Nj_sample;
}
parameters {
  real<lower=0> sigma_beta0;
  real<lower=0> sigma_beta1;
  real<lower=0> sigma_y;
  real alpha0;
  real gamma0;
  real alpha1;
  real gamma1;
  vector[K] eta0;
  vector[K] eta1;
  real mu;
  real<lower=0> sigma;
}
transformed parameters {
  vector[K] beta0;
  vector[K] beta1;
  vector[n] ymean;
  real<lower=0> mu_star;
  
  beta0 = alpha0 + gamma0 * log_Nj_sample + eta0 * sigma_beta0;
  beta1 = alpha1 + gamma1 * log_Nj_sample + eta1 * sigma_beta1;
  
  mu_star = mu + sigma^2;
  
  for (i in 1:n) {
    ymean[i] = beta0[cluster_id[i]] + beta1[cluster_id[i]] * x[i];
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
  eta0 ~ normal(0, 1);
  eta1 ~ normal(0, 1);
  y ~ normal(ymean, sigma_y);
  Nj_sample ~ lognormal(mu_star, sigma);
}
