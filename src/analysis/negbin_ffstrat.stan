functions {
  // draw new cluster sizes using nb model
  matrix Mj_new_nb_rng(int J, int K, vector Mj_sample, vector Nj_sample,
                       int[] is_strat2, vector mu, vector phi, real a, real b) {

    // need to make Mj_new_tmp b/c negbin rng function can only work w/ ints
    int Mj_new_tmp[J-K];
    vector[J] Mj_new;
    vector[J] Nj_new;
    matrix[J, 2] Mj_Nj;

    for (k in 1:(J-K)) {
      Mj_new_tmp[k] = neg_binomial_2_rng(mu[is_strat2[k]+1], phi[is_strat2[k]+1]);
      while (Mj_new_tmp[k] == 0) {
        Mj_new_tmp[k] = neg_binomial_2_rng(mu[is_strat2[k]+1], phi[is_strat2[k]+1]);
      }
    }

    Mj_new[1:K] = Mj_sample;
    Mj_new[(K+1):J] = to_vector(Mj_new_tmp);

    Nj_new[1:K] = Nj_sample;
    Nj_new[(K+1):J] = a + b*Mj_new[(K+1):J];

    Mj_Nj[,1] = Mj_new;
    Mj_Nj[,2] = Nj_new;

    return Mj_Nj;
  }
  // estimate ybar using drawn cluster sizes
  real ybar_new_nb_rng(int J, int K, vector xbar_pop, vector is_strat2,
                       vector beta0, vector beta1,
                       real alpha0, real gamma0, real kappa0,
                       real alpha1, real gamma1, real kappa1,
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
      beta0_new[j] = normal_rng(alpha0 + gamma0 * log_Mj_new[j] + kappa0 * is_strat2[j], sigma_beta0);
      beta1_new[j] = normal_rng(alpha1 + gamma1 * log_Mj_new[j] + kappa1 * is_strat2[j], sigma_beta1);
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
  int Mj_sample[K];
  int Nj_sample[K];
  vector[K] log_Mj_sample;
  vector[K] is_strat2; // 0/1 indicators for being in strat 2
}
transformed data {
  int Mj_sample_minus1[K];
  for (k in 1:K) {
    Mj_sample_minus1[k] = Mj_sample[k] - 1;
  }
}
parameters {
  real<lower=0> sigma_beta0;
  real<lower=0> sigma_beta1;
  real<lower=0> sigma_y;
  real alpha0;
  real gamma0;
  real kappa0;
  real alpha1;
  real gamma1;
  real kappa1;
  vector[K] beta0;
  vector[K] beta1;
  vector<lower=0>[2] mu;
  // reparameterize phi in terms of sqrt of CV of gamma dist
  vector<lower=0>[2] recip_sqrt_alpha_nb;
  // parameters for regressing Nj on Mj
  real a;
  real b;
  real<lower=0> tau;
}
transformed parameters {
  vector[n] ymean;
  vector<lower=0>[2] phi;
  vector<lower=0>[2] mu_star;
  vector<lower=0>[2] phi_star;

  for (s in 1:2) {
    phi[s] = 1/(recip_sqrt_alpha_nb[s])^2;
    phi_star[s] = phi[s] + 1;
    mu_star[s] = mu[s] + (mu[s] / phi[s]);
  }

  for (i in 1:n) {
    ymean[i] = beta0[cluster_id[i]] + beta1[cluster_id[i]] * x[i];
  }
}
model {
  recip_sqrt_alpha_nb ~ exponential(1); // it can be that this is close to 1
  tau ~ cauchy(0, 2.5);
  sigma_beta0 ~ cauchy(0, 2.5);
  sigma_beta1 ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 10);
  gamma0 ~ normal(0, 10);
  kappa0 ~ normal(0, 10);
  alpha1 ~ normal(0, 10);
  gamma1 ~ normal(0, 10);
  kappa1 ~ normal(0, 10);
  beta0 ~ normal(alpha0 + gamma0 * log_Mj_sample + kappa0 * is_strat2, sigma_beta0);
  beta1 ~ normal(alpha1 + gamma1 * log_Mj_sample + kappa1 * is_strat2, sigma_beta1);
  y ~ normal(ymean, sigma_y);
  for (s in 1:2) {
    Mj_sample_minus1[s] ~ neg_binomial_2(mu_star[s], phi_star[s]);
  }
  to_vector(Nj_sample) ~ normal(a + b * to_vector(Mj_sample), tau);
}

