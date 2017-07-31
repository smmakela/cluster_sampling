functions {
  vector Nj_new_bb_rng(int J, int K, int M,
                       vector Nj_sample, vector Nj_unique,
                       vector phi_star) {

    int n_star[M]; // counts of sizes for nonsampled clusters
    int count;
    int start;
    int stop;
    vector[J-K] Nj_mis;
    vector[J] Nj_new;
    real N_tot_est;
  
    // draw counts of how many times each distinct cluster size is repeated
    n_star = multinomial_rng(phi_star, J - K);
  
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
  real ybar_new_bb_rng(int J, int K,
                       vector beta0,
                       real alpha0, real gamma0,
                       real sigma_beta0,
                       vector Nj_new) {
    real N_tot_new;
    vector[J] log_Nj_new;
    vector[J] theta_new;
    vector[J] beta0_new;
    real ybar_new;

    N_tot_new = sum(Nj_new);
    log_Nj_new = log(Nj_new) - mean(log(Nj_new));
  
    beta0_new[1:K] = beta0;
    for (j in 1:K) {
      theta_new[j] = inv_logit(beta0[j]);
    }
 
    // for unsampled clusters, need to first draw new beta0, beta1 from their posteriors
    for (j in (K+1):J) {
      beta0_new[j] = normal_rng(alpha0 + gamma0 * log_Nj_new[j], sigma_beta0);
      theta_new[j] = inv_logit(beta0_new[j]); 
    }

    ybar_new = sum(theta_new .* Nj_new) / sum(Nj_new);

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
  int<lower=0,upper=1> y[n];             // outcomes
  vector[K] Nj_sample;     // vector of cluster sizes for sampled clusters
  vector[K] log_Nj_sample;  // log of cluster sizes
  vector[M] Nj_unique;     // the unique sampled cluster sizes
  int M_counts[M];    	   // counts of unique cluster sizes among sampled clusters
}
transformed data {
  vector[M] alpha_phi; // for prior on phi

  alpha_phi = rep_vector(1, M); // create vector of M zeros
}
parameters {
  real<lower=0> sigma_beta0;
  real alpha0;
  real gamma0;
  vector[K] beta0;
  simplex[M] phi;
}
transformed parameters {
  vector[n] y_prob;
  vector[M] phi_star_unnorm;
  simplex[M] phi_star;
  real cee; // normalizer for phi_star
  real pii;

  for (i in 1:n) {
    y_prob[i] = beta0[cluster_id[i]];
  }
  for (m in 1:M) {
    pii = K * Nj_unique[m] / N;
    phi_star_unnorm[m] = phi[m] * (1 - pii) / pii;

  }
  cee = sum(phi_star_unnorm); // calculate normalizing constant for phi_star
  phi_star = phi_star_unnorm / cee; // normalize phi_star
}
model {
  sigma_beta0 ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 10);
  gamma0 ~ normal(0, 10);
  beta0 ~ normal(alpha0 + gamma0 * log_Nj_sample, sigma_beta0);
  y ~ bernoulli_logit(y_prob);
  phi ~ dirichlet(alpha_phi);
  M_counts ~ multinomial(phi);
}

