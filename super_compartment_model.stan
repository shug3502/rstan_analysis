
data {
  int<lower=1> N; //number of observations
  real x1_obs[N];
  real x2_obs[N];
}

parameters {
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0,upper=1> nu;
  real<lower=0> sigma;
}

transformed parameters {
  real theta;
  theta = 2*nu - 1;
}

model {
  //set priors
  a ~ normal(0,100);
  b ~ normal(0,100);
  nu ~ beta(1,1);
  sigma ~ normal(0,100);
  
  //from aggregated two super-compartment model in long time limit
  for (i in 1:N){
    x2_obs[i] ~ normal(8*a+b*(1-nu)*x1_obs[i]/(b*nu), sigma);
  }
}
