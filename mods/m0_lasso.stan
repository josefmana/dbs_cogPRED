// generated with brms 2.17.0
functions {
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
  /* compute the logm1 link
   * Args:
   *   p: a positive scalar
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real logm1(real y) {
     return log(y - 1);
   }
  /* compute the inverse of the logm1 link
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a positive scalar
   */
   real expp1(real y) {
     return exp(y) + 1;
   }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=-1,upper=2> cens[N];  // indicates censoring
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for the lasso prior
  real<lower=0> lasso_df;  // prior degrees of freedom
  real<lower=0> lasso_scale;  // prior scale
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  // lasso shrinkage parameter
  real<lower=0> lasso_inv_lambda;
  real<lower=0> sigma;  // dispersion parameter
  real<lower=1> nu;  // degrees of freedom or shape
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_1 = r_1[, 1];
  r_1_2 = r_1[, 2];
  lprior += normal_lpdf(Intercept | 0.3, 0.1);
  lprior += chi_square_lpdf(lasso_inv_lambda | lasso_df);
  lprior += exponential_lpdf(sigma | 1);
  lprior += gamma_lpdf(nu | 2, 0.1)
    - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(sd_1[1] | 0, 0.1)
    - 1 * normal_lccdf(0 | 0, 0.1);
  lprior += normal_lpdf(sd_1[2] | 0, 0.1)
    - 1 * normal_lccdf(0 | 0, 0.1);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 2);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];
    }
    for (n in 1:N) {
    // special treatment of censored data
      if (cens[n] == 0) {
        target += student_t_lpdf(Y[n] | nu, mu[n], sigma);
      } else if (cens[n] == 1) {
        target += student_t_lccdf(Y[n] | nu, mu[n], sigma);
      } else if (cens[n] == -1) {
        target += student_t_lcdf(Y[n] | nu, mu[n], sigma);
      }
    }
  }
  // priors including constants
  target += lprior;
  target += double_exponential_lpdf(b | 0, lasso_scale * lasso_inv_lambda);
  target += std_normal_lpdf(to_vector(z_1));
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // additionally sample draws from priors
  real prior_Intercept = normal_rng(0.3,0.1);
  real prior_lasso_inv_lambda = chi_square_rng(lasso_df);
  real prior_sigma = exponential_rng(1);
  real prior_nu = gamma_rng(2,0.1);
  real prior_sd_1__1 = normal_rng(0,0.1);
  real prior_sd_1__2 = normal_rng(0,0.1);
  real prior_cor_1 = lkj_corr_rng(M_1,2)[1, 2];
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
  // use rejection sampling for truncated priors
  while (prior_lasso_inv_lambda < 0) {
    prior_lasso_inv_lambda = chi_square_rng(lasso_df);
  }
  while (prior_sigma < 0) {
    prior_sigma = exponential_rng(1);
  }
  while (prior_nu < 1) {
    prior_nu = gamma_rng(2,0.1);
  }
  while (prior_sd_1__1 < 0) {
    prior_sd_1__1 = normal_rng(0,0.1);
  }
  while (prior_sd_1__2 < 0) {
    prior_sd_1__2 = normal_rng(0,0.1);
  }
}
