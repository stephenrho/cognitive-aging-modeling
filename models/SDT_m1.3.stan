// parsimonious signal detection theory model for a rating experiment (based on Selker et al. 2020)
// all parameters (d, s, a, b) are allowed to vary by age group (coded in X)

// modification of model 1 to allow groups to differ in variability

data {
  int<lower=0> N;               // n observations
  int y[N];                     // ratings
  int<lower=0> J;               // n participants
  int<lower=1,upper=J> id[N];   // participant ids
  matrix[N, 2] X;               // design matrix for fixed (group-level) effects
  vector[N] sig_trial;          // indicator for signal (1) or noise (0) trial
  int<lower=2> K;               // n categories
  int group[J];                 // indicator of group membership (0 or 1)
}

transformed data {
  vector[K-1] unb_c;
  for (k in 1:(K-1)){
    unb_c[k] = -log( (1 - (exp(log(k) - log(K))))/(exp(log(k) - log(K))) );
  }
}

parameters {
  // d
  vector[2] B_d;
  real b_d[J];
  real<lower=0> tau_d[2];

  // c (shift [b] and scale [a])
  vector[2] B_a;
  real b_a[J];
  real<lower=0> tau_a;
  
  vector[2] B_b;
  real b_b[J];
  real<lower=0> tau_b;
  
  // s
  vector[2] B_s;
  real b_s[J];
  real<lower=0> tau_s;
}

transformed parameters {
  real<lower=0> d[N];
  real<lower=0> a[N]; 
  real b[N];
  vector[K-1] c[N];
  real<lower=0> s[N];
  simplex[K] theta[N];
  
  for (i in 1:N){
    // observation level parameters
    d[i] = exp( dot_product(X[i,], B_d) + b_d[id[i]] );
    a[i] = exp( dot_product(X[i,], B_a) + b_a[id[i]] );
    b[i] = dot_product(X[i,], B_b) + b_b[id[i]];
    s[i] = exp( dot_product(X[i,], B_s) + b_s[id[i]] );
    
    c[i] = a[i]*unb_c + b[i];
    
    // rating probabilities under SDT
    theta[i,1] = normal_cdf(c[i,1], d[i]*sig_trial[i], (1 + (-1 + s[i])*sig_trial[i]));
    for (k in 2:(K-1)){
      theta[i,k] = normal_cdf(c[i,k], d[i]*sig_trial[i], (1 + (-1 + s[i])*sig_trial[i])) - sum(theta[i,1:(k-1)]);
    }
    theta[i,K] = 1 - sum(theta[i,1:(K-1)]);
  }
}

model {
  // priors
  B_d[1] ~ normal(0, 1);
  B_d[2] ~ normal(0, 0.5);
  B_a[1] ~ normal(0, 1);
  B_a[2] ~ normal(0, 0.5);
  B_b[1] ~ normal(0, 2);
  B_b[2] ~ normal(0, 1);
  B_s[1] ~ normal(0, 0.5);
  B_s[2] ~ normal(0, 0.25);
  
  tau_d ~ cauchy(0, 1);
  tau_a ~ cauchy(0, 1);
  tau_b ~ cauchy(0, 2);
  tau_s ~ cauchy(0, 0.5);
  
  // individual level deviations
  for (j in 1:J){
    b_d[j] ~ normal(0, tau_d[group[j]+1]);
  }
  b_a ~ normal(0, tau_a);
  b_b ~ normal(0, tau_b);
  b_s ~ normal(0, tau_s);
  
  // likelihood
  for (i in 1:N){
    y[i] ~ categorical(theta[i]);
  }
}

generated quantities {
  vector[N] log_lik;      // log likelihood matrix
  vector[N] y_rep;        // posterior predictions
  for (i in 1:N){
    log_lik[i] = categorical_lpmf(y[i] | theta[i]);
    y_rep[i] = categorical_rng(theta[i]);
  }
}
