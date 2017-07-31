functions {
  // draw new cluster sizes using ln model
  vector Nj_new_ln_rng(int J, int K, vector Nj_sample,
                       real mu, real sigma) {

    vector[J] Nj_new;

    Nj_new[1:K] = Nj_sample;
    for (k in (K+1):J) {
      Nj_new[k] = lognormal_rng(mu, sigma);
    }

    return Nj_new;
  }
  // estimate ybar using drawn cluster sizes
  real ybar_new_ln_rng(int J, int K,
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
  vector[K] Nj_sample;
  vector[K] log_Nj_sample;
}
transformed data {
  vector[K] log_Nj_sample_scaled;
  log_Nj_sample_scaled = (log(Nj_sample) - mean(log(Nj_sample))) / sd(log(Nj_sample));
}
parameters {
  real<lower=0> sigma_beta0;
  real alpha0;
  real gamma0;
  vector[K] beta0;
  real mu_star_scaled;
  real<lower=0> sigma_scaled;
}
transformed parameters {
  vector[n] y_prob;
  real<lower=0> mu_star;
  real mu;
  real<lower=0> sigma;

  mu_star = mu_star_scaled + mean(log(Nj_sample));
  sigma = sigma_scaled * sd(log(Nj_sample));
  mu = mu_star - sigma^2;


  for (i in 1:n) {
    y_prob[i] = beta0[cluster_id[i]];
  }
}
model {
  sigma_scaled ~ normal(0, 1);
  mu_star_scaled ~ normal(0, 1);
  sigma_beta0 ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 10);
  gamma0 ~ normal(0, 10);
  beta0 ~ normal(alpha0 + gamma0 * log_Nj_sample, sigma_beta0);
  y ~ bernoulli_logit(y_prob);
  log_Nj_sample ~ normal(mu_star_scaled, sigma_scaled);
}

