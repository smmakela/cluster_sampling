data {
  int<lower=0> N;  # number of data points
  vector[N] y;     # outcome
  vector[N] x;     # covariate
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  y ~ normal(alpha + beta * x, sigma);
}
