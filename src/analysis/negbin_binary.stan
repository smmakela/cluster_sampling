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
  vector[K] eta0;
  real<lower=0> mu;
  real<lower=0> phi;
}
transformed parameters {
  vector[K] beta0;
  vector[n] y_prob;
  real<lower=0> mu_star;
  real<lower=0> phi_star;

  phi_star = phi + 1;
  mu_star = mu + (mu / phi);

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
  Nj_sample_minus1 ~ neg_binomial_2(mu_star, phi_star);
}
// generated quantities {
//   // here we get an estimate of y_bar, the finite-pop mean
//   // for sampled clusters, we can use the estimate 
//   // for unsampled clusters, need to draw beta0, beta1, Nj, and then calculate
//   // y_bar
//   vector[J] beta0_new;
//   vector[J] beta1_new;
//   int Nj_new[J];
//   vector[J] yj_new;
//   real ybar_new;
//   
//   beta0_new[1:K] = beta0;
//   beta1_new[1:K] = beta1;
//   Nj_new[1:K]    = Nj_sample;
//   yj_new[1:K] = beta0 + beta1 .* xbar_pop[1:K];
//   
//   for (j in (K+1):J) {
//     Nj_new[j] = neg_binomial_2_rng(mu_star, phi_star);
//     beta0_new[j] = normal_rng(alpha0 + gamma0 * log(Nj_new[j]), sigma_beta0);
//     beta1_new[j] = normal_rng(alpha1 + gamma1 * log(Nj_new[j]), sigma_beta1);
//     yj_new[j]  = normal_rng(beta0_new[j] + beta1_new[j] * xbar_pop[j], sigma_y/sqrt(Nj_new[j]));
//   }
//   
//   ybar_new = sum(yj_new .* to_vector(Nj_new)) / sum(Nj_new);
//   
// }
