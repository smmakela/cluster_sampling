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
parameters {
  real<lower=0> sigma_beta0;
  real alpha0;
  real gamma0;
  vector[K] eta0;
  real mu;
  real<lower=0> sigma;
}
transformed parameters {
  vector[K] beta0;
  vector[n] y_prob;
  real<lower=0> mu_star;

  mu_star = mu + sigma^2;

  beta0 = alpha0 + gamma0 * log_Nj_sample + eta0 * sigma_beta0;

  for (i in 1:n) {
    y_prob[i] = beta0[cluster_id[i]];
  }
}
model {
  sigma_beta0 ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 1);
  gamma0 ~ normal(0, 1);
  eta0 ~ normal(0, 1);
  y ~ bernoulli_logit(y_prob);
  Nj_sample ~ lognormal(mu_star, sigma);
}
// generated quantities {
//   // here we get an estimate of y_bar, the finite-pop mean
//   // for sampled clusters, we can use the estimate 
//   // for unsampled clusters, need to draw beta0, beta1, Nj, and then calculate
//   // y_bar
//   vector[J] beta0_new;
//   vector[J] beta1_new;
//   vector[J] Nj_new;
//   vector[J] yj_new;
//   real ybar_new;
//   
//   beta0_new[1:K] = beta0;
//   beta1_new[1:K] = beta1;
//   Nj_new[1:K]    = Nj_sample;
//   yj_new[1:K] = beta0 + beta1 .* xbar_pop[1:K];
//   
//   for (j in (K+1):J) {
//     Nj_new[j] = lognormal_rng(mu_star, sigma);
//     beta0_new[j] = normal_rng(alpha0 + gamma0 * log(Nj_new[j]), sigma_beta0);
//     beta1_new[j] = normal_rng(alpha1 + gamma1 * log(Nj_new[j]), sigma_beta1);
//     yj_new[j]  = normal_rng(beta0_new[j] + beta1_new[j] * xbar_pop[j], sigma_y/sqrt(Nj_new[j]));
//   }
//   
//   ybar_new = sum(yj_new .* Nj_new) / sum(Nj_new);
//   
// }
// 
