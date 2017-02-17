data {
  int<lower=0> N;  # number of data points
  int<lower=0> J;  # number of groups
  int inds[N];     # vector of group indices
  vector[J] Mj;    # size of groups
  vector[J] logMj; # log size of groups, centered
  vector[N] y;     # outcome
  vector[N] x;     # covariate
  vector[J] xbar;  # mean of x by group
}
parameters {
  real alpha0;
  real alpha1;
  real gamma0;
  real gamma1;
  vector[J] beta0;
  vector[J] beta1;
  real<lower=0> sigma_beta0;
  real<lower=0> sigma_beta1;
  real<lower=0> sigma_y;
}
model {
  alpha0 ~ normal(0, 1);
  gamma0 ~ normal(0, 1);
  sigma_beta0 ~ cauchy(0, 2.5);
  alpha1 ~ normal(0, 1);
  gamma1 ~ normal(0, 1);
  sigma_beta1 ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);
  beta0 ~ normal(alpha0 + gamma0 * logMj, sigma_beta0);
  beta1 ~ normal(alpha1 + gamma1 * logMj, sigma_beta1);
  for (i in 1:N) {
    y[i] ~ normal(beta0[inds[i]] + beta1[inds[i]] * x[i], sigma_y);
  }
}
generated quantities {
  vector[J] y_new;
  real ybar;
  
  y_new = beta0 + beta1 .* xbar;
  ybar = sum(Mj .* y_new) / sum(Mj);
}
