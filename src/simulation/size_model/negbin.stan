data {
  int K;
  int Nj_sample[K];
}
parameters {
  real<lower=0> phi;
  real<lower=0> mu;
}
transformed parameters {
  real<lower=0> mu_star;
  real<lower=0> phi_star;
  real<lower=0> k;
  real<lower=0> p;
  
  k = phi;
  p = phi / (phi + mu);
  
  phi_star = k + 1;
  mu_star = (k + 1) * (1 - p) / p;
}
model {
  Nj_sample ~ neg_binomial_2(mu_star, phi_star);
}
