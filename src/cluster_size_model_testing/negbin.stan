data {
  int J; // number of pop clusters
  int K; // number of sampled clusters
  int Nj_sample_minus1[K];
}
parameters {
  real<lower=0> mu;
  real<lower=0> phi;
}
transformed parameters {
  real<lower=0> mu_star;
  real<lower=0> phi_star;
  
  phi_star = phi + 1;
  mu_star = mu + (mu / phi);
}
model {
  Nj_sample_minus1 ~ neg_binomial_2(mu_star, phi_star);
}
generated quantities {
  int Nj_pop_new[J];
  int Nj_sample_new[K];
  for (j in 1:J) {
    Nj_pop_new[j] = neg_binomial_2_rng(mu, phi);
  }
  for (k in 1:K) {
    Nj_sample_new[k] = neg_binomial_2_rng(mu_star, phi_star);
  }
}
// unbiasedness is really in terms of drawing from the prior and averaging that way
// negbin as mixture of poissons -- see paper w/ Tian
