data {
  int K;
  vector[K] Nj_sample_minus1;
}
parameters {
  real beta; # hyperprior for mu_star
  real<lower=0> tau; # hyperprior for mu_star
  real eta;
  real<lower=0> sigma;
}
transformed parameters {
  real<lower=0> mu_star;
  mu_star = beta + eta * tau; # mu_star ~ N(beta, tau)
}
model {
  beta ~ normal(0, 1);
  tau ~ cauchy(0, 2.5);
  eta ~ normal(0, 1);
  sigma ~ cauchy(0, 5);
  Nj_sample_minus1 ~ normal(mu_star, sigma);
}
generated quantities {
  real phi_star; # no restrictions so stan doesn't crash in case they don't hold during warmup
  real phi;
  real mu;
  phi_star = mu_star^2 / (sigma^2 - mu_star);
  phi = phi_star - 1;
  mu = mu_star * (phi_star - 1) / phi_star;
}
