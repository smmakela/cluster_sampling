functions {
  // estimate ybar using drawn cluster sizes
  real ybar_new_inds_rng(int J, int K, vector xbar_pop,
                         vector beta0, vector beta1,
                         real alpha0, real alpha1,
                         real sigma_beta0, real sigma_beta1,
                         real sigma_y, vector Nj_pop) {
    vector[J] beta0_new;
    vector[J] beta1_new;
    vector[J] yj_new;
    real ybar_new;
    
    beta0_new[1:K] = beta0;
    beta1_new[1:K] = beta1;
    yj_new[1:K] = beta0 + beta1 .* xbar_pop[1:K];
  
    for (j in (K+1):J) {
      beta0_new[j] = normal_rng(alpha0, sigma_beta0);
      beta1_new[j] = normal_rng(alpha1, sigma_beta1);
      yj_new[j]  = normal_rng(beta0_new[j] + beta1_new[j] * xbar_pop[j],
                              sigma_y/sqrt(Nj_pop[j]));
    }
    
    ybar_new = sum(yj_new .* to_vector(Nj_pop)) / sum(Nj_pop);

    return ybar_new; 
  } 
} # end functions block
data {
  int<lower=0> J;         // number of clusters in population
  int<lower=0> K;         // number of clusters in sample
  int<lower=0> n;         // sample size
  int cluster_id[n]; // cluster ids for sample units
  vector[n] x;            // individual-level covar ("age")
  vector[n] y;            // outcomes
}
parameters {
  real<lower=0> sigma_beta0;
  real<lower=0> sigma_beta1;
  real<lower=0> sigma_y;
  real alpha0;
  real alpha1;
  vector[K] beta0;
  vector[K] beta1;
}
transformed parameters {
  vector[n] ymean;

  for (i in 1:n) {
    ymean[i] = beta0[cluster_id[i]] + beta1[cluster_id[i]]*x[i];
  }
}
model {
  sigma_beta0 ~ cauchy(0, 2.5);
  sigma_beta1 ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 10);
  alpha1 ~ normal(0, 10);
  beta0 ~ normal(alpha0, sigma_beta0);
  beta1 ~ normal(alpha1, sigma_beta1);
  y ~ normal(ymean, sigma_y);
}

