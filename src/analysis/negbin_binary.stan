functions {
  // draw new cluster sizes using nb model
  vector Nj_new_nb_rng(int J, int K, vector Nj_sample,
                       real mu, real phi) {

    // need to make Nj_new_tmp b/c negbin rng function can only work w/ ints
    int Nj_new_tmp[J-K];
    vector[J] Nj_new;

    for (k in 1:(J-K)) {
      Nj_new_tmp[k] = neg_binomial_2_rng(mu, phi);
      while (Nj_new_tmp[k] == 0) {
        Nj_new_tmp[k] = neg_binomial_2_rng(mu, phi);
      }
    }

    Nj_new[1:K] = Nj_sample;
    Nj_new[(K+1):J] = to_vector(Nj_new_tmp);

    return Nj_new;
  }
  // estimate ybar using drawn cluster sizes
  real ybar_new_nb_rng(int J, int K,
                       vector beta0,
                       real alpha0, real gamma0,
                       real sigma_beta0,
                       vector Nj_new) {

    vector[J] beta0_new;
    vector[J] log_Nj_new;
    vector[J] theta_new;
    real ybar_new;
    
    log_Nj_new = log(Nj_new) - mean(log(Nj_new));

    beta0_new[1:K] = beta0;
    for (j in 1:K) {
      theta_new[j] = inv_logit(beta0[j]);
    }
  
    for (j in (K+1):J) {
      beta0_new[j] = normal_rng(alpha0 + gamma0 * log_Nj_new[j], sigma_beta0);
      theta_new[j] = inv_logit(beta0_new[j]);
    }
    
    ybar_new = sum(theta_new .* to_vector(Nj_new)) / sum(Nj_new);

    return ybar_new; 
  } 
} # end functions block
data {
  int J; // number of pop clusters
  int K; // number of sampled clusters
  int n; // total sample size
  int<lower=0,upper=1> y[n];             // outcomes
  int cluster_id[n]; // vector of cluster id's for each sampled unit
  int Nj_sample[K];
  vector[K] log_Nj_sample;
}
transformed data {
  int Nj_sample_minus1[K];
  for (k in 1:K) {
    Nj_sample_minus1[k] = Nj_sample[k] - 1;
  }
}
parameters {
  real<lower=0> sigma_beta0;
  real alpha0;
  real gamma0;
  vector[K] beta0;
  real<lower=0> mu;
  // reparameterize phi in terms of sqrt of CV of gamma dist
  real<lower=0> recip_sqrt_alpha_nb;
}
transformed parameters {
  vector[n] y_prob;
  real<lower=0> phi;
  real<lower=0> mu_star;
  real<lower=0> phi_star;

  phi = 1/(recip_sqrt_alpha_nb)^2;
  phi_star = phi + 1;
  mu_star = mu + (mu / phi);

  for (i in 1:n) {
    y_prob[i] = beta0[cluster_id[i]];
  }
}
model {
  recip_sqrt_alpha_nb ~ exponential(1);
  sigma_beta0 ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 10);
  gamma0 ~ normal(0, 10);
  beta0 ~ normal(alpha0 + gamma0 * log_Nj_sample, sigma_beta0);
  y ~ bernoulli_logit(y_prob);
  Nj_sample_minus1 ~ neg_binomial_2(mu_star, phi_star);
}

