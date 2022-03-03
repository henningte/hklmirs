// generated with brms 2.13.0
functions {

  /* Efficient computation of the horseshoe prior
   * see Appendix C.1 in https://projecteuclid.org/euclid.ejs/1513306866
   * Args:
   *   z: standardized population-level coefficients
   *   lambda: local shrinkage parameters
   *   tau: global shrinkage parameter
   *   c2: slap regularization parameter
   * Returns:
   *   population-level coefficients following the horseshoe prior
   */
  vector horseshoe(vector z, vector lambda, real tau, real c2) {
    int K = rows(z);
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return z .* lambda_tilde * tau;
  }
}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for the horseshoe prior
  real<lower=0> hs_df;  // local degrees of freedom
  real<lower=0> hs_df_global;  // global degrees of freedom
  real<lower=0> hs_df_slab;  // slab degrees of freedom
  real<lower=0> hs_scale_global;  // global prior scale
  real<lower=0> hs_scale_slab;  // slab prior scale
  int prior_only;  // should the likelihood be ignored?
  // data for peat samples
  int<lower = 1> N1; // no of peat samples
  int<lower = 1> N1_cores; // no of peat cores
  matrix[N1, K] X1; // spectral data
  vector[N1_cores] latitudes; // latitudes of peat cores
  int index_core[N1]; // indicator assigning samples to cores

  // parameter values for priors
  real<lower = 0> a_latitude_sigma;
  real a_latitude_mu;
  real<lower = 0> b_latitude_sigma;
  real<lower = 0> phi_pop_a;
  real<lower = 0> phi_pop_b;
  real<lower = 0> phi_pop_scale;
  real<lower = 0> phi_core_a;
  real<lower = 0> phi_core_b;
  real<lower = 0> phi_core_scale;
  real<lower = 0> phi_a;
  real<lower = 0> phi_b;
  real<lower = 0> phi_scale;
  real<lower = 0> Intercept_sigma;
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering

  matrix[N1, Kc] Xc1;

  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
    Xc1[, i - 1] = X1[, i] - means_X[i - 1];
  }
}
parameters {
  // local parameters for horseshoe prior
  vector[Kc] zb;
  vector<lower=0>[Kc] hs_local;
  real Intercept;  // temporary intercept for centered predictors
  // horseshoe shrinkage parameters
  real<lower=0> hs_global;  // global shrinkage parameters
  real<lower=0> hs_slab;  // slab regularization parameter
  real<lower=0> phi;  // precision parameter
  // parameters for core averages
  real<lower = 0, upper = 1> Y1[N1]; // predictions for peat samples
  vector<lower = 0>[N1_cores] phi_core; // standard deviation for content for each peat core
  // parameters for model for latitudinal gradient
  real a_latitude;
  real b_latitude;
  vector<lower = 0, upper = 1>[N1_cores] mu_core;
  real<lower = 0> phi_pop;

}
transformed parameters {

  vector[N1] mu_pred;
  vector[N] mu;
  vector[N1_cores] mu_pop;
  //vector[N1_cores] mu_pop; // average content for each peat core

  vector[Kc] b;  // population-level effects
  // compute actual regression coefficients
  b = horseshoe(zb, hs_local, hs_global, hs_scale_slab^2 * hs_slab);

  mu_pred = Intercept * Intercept_sigma + Xc1 * b;
  for(n in 1:N1) {
    // apply the inverse link function
    mu_pred[n] = inv_logit(mu_pred[n]);
  }

  mu = Intercept * Intercept_sigma + Xc * b;
  mu_pop = a_latitude * a_latitude_sigma + a_latitude_mu + latitudes * b_latitude * b_latitude_sigma;
  for (n in 1:N) {
    // apply the inverse link function
    mu[n] = inv_logit(mu[n]);
  }
  for(n in 1:N1_cores) {
    // apply the inverse link function
    mu_pop[n] = inv_logit(mu_pop[n]);
  }

}
model {
  // priors including all constants
  target += std_normal_lpdf(zb);
  target += student_t_lpdf(hs_local | hs_df, 0, 1)
    - rows(hs_local) * log(0.5);
  target += normal_lpdf(Intercept | 0, 1);
  target += student_t_lpdf(hs_global | hs_df_global, 0, hs_scale_global)
    - 1 * log(0.5);
  target += inv_gamma_lpdf(hs_slab | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
  target += gamma_lpdf(phi | phi_a, phi_b);

  //target += gamma_lpdf(phi1 | 0.01, 0.01);
  // likelihood including all constants
  if (!prior_only) {
    target += beta_lpdf(Y | mu * phi * phi_scale, (1 - mu) * phi * phi_scale);
    Y1 ~ beta(mu_pred * phi * phi_scale, (1 - mu_pred) * phi * phi_scale);
  }

 // average contents for each core
  phi_core ~ gamma(phi_core_a, phi_core_b); // 50, 10
  phi_pop ~ gamma(phi_pop_a, phi_pop_b); // 90, 20
  for(n in 1:N1) {
    Y1[n] ~ beta(mu_core[index_core[n]] * phi_core[index_core[n]] * phi_core_scale, (1 - mu_core[index_core[n]]) * phi_core[index_core[n]] * phi_core_scale);
  }

  // model for latitudinal gradient
  a_latitude ~ normal(0, 1);
  b_latitude ~ normal(0, 1);
  mu_core ~ beta(mu_pop * phi_pop * phi_pop_scale, (1 - mu_pop) * phi_pop * phi_pop_scale);

}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept * Intercept_sigma - dot_product(means_X, b);
  real Y_rep[N] = beta_rng(mu * phi * phi_scale, (1 - mu) * phi * phi_scale);
}
