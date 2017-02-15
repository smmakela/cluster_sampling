data {
  int J;
  int K;
  real Nj_sample[K];
}
parameters {
  real mu;
  real<lower=0> sigma;
}
transformed parameters {
  real mu_star;
  mu_star = mu + sigma^2;
}
model {
  Nj_sample ~ lognormal(mu_star, sigma);
}
generated quantities {
  real Nj_pop_new[J];
  real Nj_sample_new[K];
  for (j in 1:J) {
    Nj_pop_new[j] = lognormal_rng(mu, sigma);
  }
  for (k in 1:K) {
    Nj_sample_new[k] = lognormal_rng(mu_star, sigma);
  }
}
