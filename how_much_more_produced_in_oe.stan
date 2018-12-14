data {
  int<lower=1> T1; //how many WT data pts
  int<lower=1> T2; //how many OE data pts
  real x[T1]; //wt observations
  real y[T2]; //oe observations
  vector[T1] ts1; //time series for observations
  vector[T2] ts2;
}
parameters {
  real<lower=0> alpha[2]; //growth rates
  real<lower=0> sigma[2];
}
transformed parameters { 
  real gamma = alpha[2]/alpha[1];
}
model {
  //priors
  for (i in 1:2){
    alpha[i] ~ normal(0,10) T[0,];
    sigma[i] ~ cauchy(0,3) T[0,];
  }
  x ~ normal(15*alpha[1]*ts1,sigma[1]);
  y ~ normal(15*alpha[2]*ts2,sigma[2]);
}

