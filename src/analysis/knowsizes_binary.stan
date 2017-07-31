functions {
  // estimate ybar using drawn cluster sizes
  real ybar_new_know_rng(int J, int K,
                         vector beta0,
                         real alpha0, real gamma0,
                         real sigma_beta0, 
                         vector Nj_pop) {

    vector[J] beta0_new;
    vector[J] log_Nj_pop;
    vector[J] theta_new;
    real ybar_new;
    
    log_Nj_pop = log(Nj_pop) - mean(log(Nj_pop));

    beta0_new[1:K] = beta0;
    for (j in 1:K) {
      theta_new[j] = inv_logit(beta0[j]);
    }
  
    for (j in (K+1):J) {
      beta0_new[j] = normal_rng(alpha0 + gamma0 * log_Nj_pop[j], sigma_beta0);
      theta_new[j] = inv_logit(beta0_new[j]);
    }
    
    ybar_new = sum(theta_new .* to_vector(Nj_pop)) / sum(Nj_pop);

    return ybar_new; 
  } 
} # end functions block
data {
  int<lower=0> J;         // number of clusters in population
  int<lower=0> K;         // number of clusters in sample
  int<lower=0> n;         // sample size
  int cluster_id[n]; // cluster ids for sample units
  int<lower=0,upper=1> y[n];             // outcomes
  vector[K] Nj_sample;       // vector of cluster sizes for sampled clusters
  vector[K] log_Nj_sample;    // log of cluster sizes for sampled clusters
}
parameters {
  real<lower=0> sigma_beta0;
  real alpha0;
  real gamma0;
  vector[K] beta0;
}
transformed parameters {
  vector[n] y_prob;

  for (i in 1:n) {
    y_prob[i] = beta0[cluster_id[i]];
  }
}
model {
  sigma_beta0 ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 10);
  gamma0 ~ normal(0, 10);
  beta0 ~ normal(alpha0 + gamma0 * log_Nj_sample, sigma_beta0);
  y ~ bernoulli_logit(y_prob);
}

