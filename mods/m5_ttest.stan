//
// Bayesian "t-test" for between-group mean differences with unequel variance
//

data {
  
  // observations
  int<lower=1> N0; // number of observations in group without MRI
  int<lower=1> N1; // number of observations in group with MRI
  vector[N0] y0; // response vector of group without MRI
  vector[N1] y1; // response vector of group with MRI

}

parameters {

  real mu_0;
  real mu_1;
  real<lower=0> sigma_0;
  real<lower=0> sigma_1;

}

model {

  // likelihoods
  target += normal_lpdf( y0 | mu_0, sigma_0 );
  target += normal_lpdf( y1 | mu_1, sigma_1 );
  
  // priors
  target += normal_lpdf( mu_0 | 0, 1 );
  target += normal_lpdf( mu_1 | 0, 1 );
  target += exponential_lpdf( sigma_0 | 1 );
  target += exponential_lpdf( sigma_1 | 1 );
  
}