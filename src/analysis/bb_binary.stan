functions {
  vector Mj_new_bb_rng(int J, int K, int num_uniq_sz,
                       vector Mj_sample, vector vec_uniq_sz,
                       vector phi_star) {

    int n_star[num_uniq_sz]; // counts of sizes for nonsampled clusters
    int count;
    int start;
    int stop;
    vector[J-K] Mj_mis;
    vector[J] Mj_new;
    real N_tot_est;
  
    // draw counts of how many times each distinct cluster size is repeated
    n_star = multinomial_rng(phi_star, J - K);
  
    // populate the vector of unsampled sizes using the counts in n_star
    // because R-like subsetting on the lhs is tricky in stan, we have to
    // manually set the start and stop indices
    count = 1;
    for (m in 1:num_uniq_sz) {
      start = count;
      // define stop such that stop - start has n_star[m] elements inclusive
      stop = start + n_star[m] - 1;
      for (p in start:stop) {
        Mj_mis[p] = vec_uniq_sz[m];
      }
      count = start + n_star[m];
    }
    Mj_new = append_row(Mj_sample, Mj_mis);

    return(Mj_new);
  }
  real ybar_new_bb_rng(int J, int K,
                       vector beta0,
                       real alpha0, real gamma0,
                       real sigma_beta0,
                       vector Mj_new, vector Nj_new) {
    real N_tot_new;
    vector[J] log_Mj_new;
    vector[J] theta_new;
    vector[J] beta0_new;
    real ybar_new;

    N_tot_new = sum(Mj_new);
    log_Mj_new = log(Mj_new) - mean(log(Mj_new));
  
    beta0_new[1:K] = beta0;
    for (j in 1:K) {
      theta_new[j] = inv_logit(beta0[j]);
    }
 
    // for unsampled clusters, need to first draw new beta0, beta1 from their posteriors
    for (j in (K+1):J) {
      beta0_new[j] = normal_rng(alpha0 + gamma0 * log_Mj_new[j], sigma_beta0);
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
  int<lower=0> M_tot;          // population size
  int<lower=0> num_uniq_sz;	   // number of unique cluster sizes
  int cluster_id[n];       // renumbered cluster ids, length = n
  int<lower=0,upper=1> y[n];             // outcomes
  vector[K] Mj_sample;     // vector of cluster sizes for sampled clusters
  vector[K] log_Mj_sample;  // log of cluster sizes
  vector[num_uniq_sz] vec_uniq_sz;     // the unique sampled cluster sizes
  int cts_uniq_sz[num_uniq_sz];    	   // counts of unique cluster sizes among sampled clusters
}
transformed data {
  vector[num_uniq_sz] alpha_phi; // for prior on phi

  alpha_phi = rep_vector(1, num_uniq_sz); // create vector of num_uniq_sz zeros
}
parameters {
  real<lower=0> sigma_beta0;
  real alpha0;
  real gamma0;
  vector[K] beta0;
  simplex[num_uniq_sz] phi;
}
transformed parameters {
  vector[n] y_prob;
  vector[num_uniq_sz] phi_star_unnorm;
  simplex[num_uniq_sz] phi_star;
  real cee; // normalizer for phi_star
  real pii;

  for (i in 1:n) {
    y_prob[i] = beta0[cluster_id[i]];
  }
  for (m in 1:num_uniq_sz) {
    pii = K * vec_uniq_sz[m] / M_tot;
    phi_star_unnorm[m] = phi[m] * (1 - pii) / pii;

  }
  cee = sum(phi_star_unnorm); // calculate normalizing constant for phi_star
  phi_star = phi_star_unnorm / cee; // normalize phi_star
}
model {
  sigma_beta0 ~ cauchy(0, 2.5);
  alpha0 ~ normal(0, 10);
  gamma0 ~ normal(0, 10);
  beta0 ~ normal(alpha0 + gamma0 * log_Mj_sample, sigma_beta0);
  y ~ bernoulli_logit(y_prob);
  phi ~ dirichlet(alpha_phi);
  cts_uniq_sz ~ multinomial(phi);
}

