functions {
  real[] mrnatransport(real t,
                       real[] y,
                       real[] theta,
                       real[] x_r,
                       int[] x_i
  ) {
    vector[16] dydt;
    matrix[16,16] B;
    vector[16] producers;
    int count;
    producers = rep_vector(1,16);
    producers[1] = 0;
    B = to_matrix(x_r,16,16);
    dydt = theta[1] * B * to_vector(y) + theta[2] * producers;
    return to_array_1d(dydt);
  }
}
data {
  int<lower=1> T;
  real y[T,16];
  real y0[16];
  real t0;
  real ts[T];
  real B1[16*16];
  real B2[16*16];
}
transformed data {
  real x_r[16*16*2];
  int x_i[0];
}
parameters {
  real<lower=0> sigma; //noise param
  real<lower=0,upper=1> phi; //difference between particles in NCs and in Oocyte
  real<lower=0> theta[3];
  real<lower=0,upper=1> q;
}
model {
  real z[T,16];
  int cell_indices[5]; //only model these cells
  sigma ~ normal(0,50) T[0,]; 
  phi ~ normal(0.289,0.0285) T[0,1];
  theta[1] ~ normal(0,10) T[0,];
  theta[2] ~ normal(0,100) T[0,];
  theta[3] ~ beta(0.1,0.1) T[0,1];
  q ~ uniform(0,1);
  if (q<theta[3]){
    z = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, B1, x_i);
  } else {
    z = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, B2, x_i);
  }
  cell_indices[1] = 1;
  cell_indices[2] = 2;
  cell_indices[3] = 3;
  cell_indices[4] = 9;
  cell_indices[5] = 5;
  for (t in 1:T){
    for (j in 1:5) {
      if (j>1){
        //y[t,j] ~ neg_binomial_2(z[t,cell_indices[j]], sigma);
        //y[t,j] ~ poisson(z[t,cell_indices[j]]);
        y[t,j] ~ normal(z[t,cell_indices[j]], sigma);
      } else {
        //y[t,j] ~ neg_binomial_2(phi*z[t,cell_indices[j]], sigma);    
        //y[t,j] ~ poisson(phi*z[t,cell_indices[j]]);    
        y[t,j] ~ normal(phi*z[t,cell_indices[j]], sigma);
      }
    }
  }
}
generated quantities {
  real y_pred[T,16];
  real y_ode[T,16];
  int q_pred;
  q_pred = bernoulli_rng(theta[3]);
  if (q_pred>0){
  y_ode = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, B1, x_i );
  } else {
  y_ode = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, B2, x_i );
  }
  for (t in 1:T){
    for (j in 1:16){
      if (j>1){
        //y_pred[t,j] = neg_binomial_2_rng(y_ode[t,j], sigma);
        //y_pred[t,j] = poisson_rng(y_ode[t,j]);
    	  y_pred[t,j] = normal_rng(y_ode[t,j], sigma);
      } else {
        //y_pred[t,j] = neg_binomial_2_rng(phi*y_ode[t,j], sigma);    
        //y_pred[t,j] = poisson_rng(phi*y_ode[t,j]);    
    	  y_pred[t,j] = normal_rng(phi*y_ode[t,j], sigma);
      }
    }
  }
}

 
