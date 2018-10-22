data {
  int<lower=1> T;
  int transcription_obs[T,16];
}
parameters {
  real<lower=0> a[16];
  real<lower=0,upper=1> p_missing;
}
model {
  a ~ cauchy(0,5);
  p_missing ~ beta(1,1);
  for (t in 1:T) {
    for (j in 1:16){
      if (transcription_obs[t,j] == 0){
        target += log_sum_exp(bernoulli_lpmf(1 | p_missing),
                      bernoulli_lpmf(0 | p_missing)
                        + poisson_lpmf(transcription_obs[t,j] | a[j]));
      } else {
        target += bernoulli_lpmf(0 | p_missing)
                + poisson_lpmf(transcription_obs[t,j] | a[j]);  
      }
    }
  }
}

