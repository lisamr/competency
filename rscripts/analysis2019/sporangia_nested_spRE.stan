data {
  int<lower=0> N; //observations
  int J; //leafIDs 
  int S; //species
  int count[N];
  int leafID[N];
  int spID[J];
  vector[N] leaf_area_cm2;
}
parameters {
  real a0;
  real<lower=0> sdID;
  real<lower=0> sdSp;
  vector[S] zSp;
  vector[J] zID;
}
transformed parameters{
  vector[S] aSp;
  vector[J] aID;
  aSp = sdSp*zSp;
  aID = aSp[spID] + sdID*zID;
}
model {
  vector[N] mu;
  a0 ~ normal(2.2, 1);
  zSp ~ normal(0, 1);
  zID ~ normal(0,1);
  sdSp ~ exponential(1);
  sdID ~ exponential(1);
  
  for(i in 1:N){
    mu[i] = exp(a0 + aID[leafID[i]] + log(leaf_area_cm2[i]));
  }
  count ~ poisson(mu);
}
generated quantities{
  vector[N] y_rep;
  vector[N] log_lik;
  
  for(i in 1:N){
    y_rep[i] = poisson_log_rng((a0 + aID[leafID[i]] + log(leaf_area_cm2[i])));
    log_lik[i] = poisson_log_lpmf(count[i]|a0 + aID[leafID[i]] + log(leaf_area_cm2[i]));
  }
}
