functions {
  // draw new cluster sizes using ln model
  vector Mj_new_ln_rng(int J, int K, vector Mj_sample,
                       real mu, real sigma) {

    vector[J] Mj_new;

    Mj_new[1:K] = Mj_sample;
    for (k in (K+1):J) {
      Mj_new[k] = lognormal_rng(mu, sigma);
    }

    return Mj_new;
  }
  // estimate ybar using drawn cluster sizes
  real ybar_new_ln_rng(int J, int K, vector xbar_pop,
                       vector beta0, vector beta1,
                       real alpha0, real gamma0,
                       real alpha1, real gamma1,
                       real sigma_beta0, real sigma_beta1,
                       real sigma_y, vector Mj_new, vector Nj_new) {

    vector[J] beta0_new;
    vector[J] beta1_new;
    vector[J] log_Mj_new;
    vector[J] yj_new;
    real ybar_new;
    
    log_Mj_new = log(Mj_new) - mean(log(Mj_new));

    beta0_new[1:K] = beta0;
    beta1_new[1:K] = beta1;
    yj_new[1:K] = beta0 + beta1 .* xbar_pop[1:K];
  
    for (j in (K+1):J) {
      beta0_new[j] = normal_rng(alpha0 + gamma0 * log_Mj_new[j], sigma_beta0);
      beta1_new[j] = normal_rng(alpha1 + gamma1 * log_Mj_new[j], sigma_beta1);
      yj_new[j]  = normal_rng(beta0_new[j] + beta1_new[j] * xbar_pop[j],
                              sigma_y/sqrt(Nj_new[j]));
    }
    
    ybar_new = sum(yj_new .* to_vector(Nj_new)) / sum(Nj_new);

    return ybar_new; 
  } 
} # end functions block
data {
  int J; // number of pop clusters
  int K; // number of sampled clusters
  int n; // total sample size
  vector[n] x;
  vector[n] y;
  int cluster_id[n]; // vector of cluster id's for each sampled unit
  vector[K] Mj_sample;
  vector[K] log_Mj_sample;
}
transformed data {
  vector[K] log_Mj_sample_scaled;
  log_Mj_sample_scaled = (log(Mj_sample) - mean(log(Mj_sample))) / sd(log(Mj_sample));
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
  real mu_star_scaled;
  real<lower=0> sigma_scaled;
}
transformed parameters {
  vector[K] beta0;
  vector[K] beta1;
  vector[n] ymean;
  real mu;
  real<lower=0> sigma;
  real<lower=0> mu_star;

  mu_star = mu_star_scaled + mean(log(Mj_sample));
  sigma = sigma_scaled * sd(log(Mj_sample));
  mu = mu_star - sigma^2;
  //mu_star = mu + sigma^2;
  //mu_star_scaled = mu_star - mean(log(Mj_sample));
  //sigma_scaled = sigma / sd(log(Mj_sample));

  beta0 = alpha0 + gamma0 * log_Mj_sample + eta0 * sigma_beta0;
  beta1 = alpha1 + gamma1 * log_Mj_sample + eta1 * sigma_beta1;
  for (i in 1:n) {
    ymean[i] = beta0[cluster_id[i]] + beta1[cluster_id[i]] * x[i];
  }
}
model {
  mu_star_scaled ~ normal(0, 1);
  sigma_scaled ~ normal(0, 1);
  sigma_beta0 ~ cauchy(0, 2.5);
  sigma_beta1 ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 10);
  gamma0 ~ normal(0, 10);
  alpha1 ~ normal(0, 10);
  gamma1 ~ normal(0, 10);
  eta0 ~ normal(0, 1);
  eta1 ~ normal(0, 1);
  y ~ normal(ymean, sigma_y);
  log_Mj_sample_scaled ~ normal(mu_star_scaled, sigma_scaled);
}

