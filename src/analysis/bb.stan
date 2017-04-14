functions {
  vector Nj_new_bb_rng(int J, int K, int M,
                       vector Nj_sample, vector Nj_unique,
                       real mu_star, vector phi_star) {

    int J_mis;
    int n_star[M]; // counts of sizes for nonsampled clusters
    int count;
    int start;
    int stop;
    vector[J_mis] Nj_mis;
    vector[J] Nj_new;
    real N_tot_est;
  
    // number of unsampled clusters
    J_mis = J - K;

    // draw counts of how many times each distinct cluster size is repeated
    n_star = multinomial_rng(phi_star, J_mis);
  
    // populate the vector of unsampled sizes using the counts in n_star
    // because R-like subsetting on the lhs is tricky in stan, we have to
    // manually set the start and stop indices
    count = 1;
    for (m in 1:M) {
      start = count;
      // define stop such that stop - start has n_star[m] elements inclusive
      stop = start + n_star[m] - 1;
      for (p in start:stop) {
        Nj_mis[p] = Nj_unique[m];
      }
      count = start + n_star[m];
    }
    Nj_new = append_row(Nj_sample, Nj_mis);

    return(Nj_new);
  }
  real ybar_new_bb_rng(int J, int K, vector xbar_pop,
                       vector beta0, vector beta1,
                       real alpha0, real gamma0,
                       real alpha1, real gamma1,
                       real sigma_beta0, real sigma_beta1, real sigma_y,
                       vector Nj_new) {
    int J_mis;
    real N_tot_new;
    vector[J] log_Nj_new;
    vector[J] y_new;
    vector[J_mis] beta0_new;
    vector[J_mis] beta1_new;
    real ybar_new;

    J_mis = J - K;
  
    N_tot_new = sum(Nj_new);
    log_Nj_new = log(Nj_new) - mean(log(Nj_new));
  
    // for the sampled clusters, use poste=rior means of beta0, beta1
    y_new[1:K] = beta0 + beta1 .* xbar_pop[1:K];
  
    // for unsampled clusters, need to first draw new beta0, beta1 from their posteriors
    for (j in 1:J_mis) {
      beta0_new[j] = normal_rng(alpha0 + gamma0 * logNj_new[K + j], sigma_beta0);
      beta1_new[j] = normal_rng(alpha1 + gamma1 * logNj_new[K + j], sigma_beta1);
    }
    y_new[(K + 1):J] = beta0_new + beta1_new .* xbar_pop[(K + 1):J];
  
    ybar_new = sum(y_new .* Nj_new) / sum(Nj_new);

    return(ybar_new);
  }
} # end functions block
data {
  int<lower=0> J;          // number of clusters
  int<lower=0> K;          // number of clusters
  int<lower=0> n;          // sample size
  int<lower=0> N;          // population size
  int<lower=0> M;	   // number of unique cluster sizes
  int cluster_id[n];       // renumbered cluster ids, length = n
  vector[n] x;             // individual-level covar ("age")
  vector[n] y;             // outcomes
  vector[K] Nj_sample;     // vector of cluster sizes for sampled clusters
  vector[K] logNj_sample;  // log of cluster sizes
  vector[M] Nj_unique;     // the unique sampled cluster sizes
  vector[J] xbar_pop;      // cluster avgs of xij
  real Tx;                 // total pop size (of exp(log(Mj)-mean(log(Mj))))
  int M_counts[M];    	   // counts of unique cluster sizes among sampled clusters
}
transformed data {
  vector[M] alpha_phi; // for prior on phi

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
  vector[K] beta0;
  vector[K] beta1;
  simplex[M] phi;
}
transformed parameters {
  vector[n] yhat;
  vector[M] phi_star_unnorm;
  vector[M] phi_star;
  real cee; // normalizer for phi_star
  real pii;

  for (i in 1:n) {
    yhat[i] = beta0[cluster_id[i]] + x[i] * beta1[cluster_id[i]];
  }
  for (m in 1:M) {
    pii = K * Nj_unique[m] / Tx;
    phi_star_unnorm[m] = phi[m] * (1 - pii) / pii;
  }
  cee = sum(phi_star_unnorm); // calculate normalizing constant for phi_star
  phi_star = phi_star_unnorm / cee; // normalize phi_star
}
model {
  sigma_beta0 ~ cauchy(0, 2.5);
  sigma_beta1 ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 1);
  gamma0 ~ normal(0, 1);
  alpha1 ~ normal(0, 1);
  gamma1 ~ normal(0, 1);
  beta0 ~ normal(alpha0 + gamma0 * logNj_sample, sigma_beta0);
  beta1 ~ normal(alpha1 + gamma1 * logNj_sample, sigma_beta1);
  y ~ normal(yhat, sigma_y);
  phi ~ dirichlet(alpha_phi);
  M_counts ~ multinomial(phi);
}

