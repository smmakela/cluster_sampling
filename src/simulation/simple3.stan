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
  vector[J] eta0;
  vector[J] eta1;
  real<lower=0> sigma_beta0;
  real<lower=0> sigma_beta1;
  real<lower=0> sigma_y;
}
  transformed parameters {
    vector[J] beta0;
    vector[J] beta1;
    real ymean[N];
    
    beta0 = alpha0 + gamma0 * logMj + eta0 * sigma_beta0;
    beta1 = alpha1 + gamma1 * logMj + eta1 * sigma_beta1;
    for (i in 1:N) {
      ymean[i] = beta0[inds[i]] + beta1[inds[i]] * x[i];
      //ymean[i] = alpha0 + gamma0 * logMj[inds[i]] + sigma_beta0 * eta0[inds[i]] + alpha1 * x[i] + gamma1 * logMj[inds[i]] * x[i] + sigma_beta1 * eta1[inds[i]]*x[i];
    }
}
  model {
    alpha0 ~ normal(0, 1);
    gamma0 ~ normal(0, 1);
    sigma_beta0 ~ cauchy(0, 2.5);
    alpha1 ~ normal(0, 1);
    gamma1 ~ normal(0, 1);
    sigma_beta1 ~ cauchy(0, 2.5);
    sigma_y ~ cauchy(0, 2.5);
    eta0 ~ normal(0, 1);
    eta1 ~ normal(0, 1);
    y ~ normal(ymean, sigma_y);
  }
generated quantities {
  vector[J] y_new;
  real ybar;
  
  y_new = alpha0 + gamma0 * logMj + alpha1 * xbar + gamma1 * xbar .* logMj;
  ybar = sum(Mj .* y_new) / sum(Mj);
}
