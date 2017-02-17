data {
  int K;
  real Nj_sample[K];
}
parameters {
  real mu;
  real<lower=0> sigma;
}
transformed parameters {
  real mu_sb;
  mu_sb = mu + sigma^2;
}
model {
  mu ~ normal(100, 100);
  sigma ~ cauchy(0, 2.5);
  Nj_sample ~ lognormal(mu_sb, sigma);
}
