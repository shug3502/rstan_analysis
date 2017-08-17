functions {
  real[] mrnatransport(real t,
                       real[] y,
                       real[] theta,
                       real[] x_r,
                       int[] x_i
  ) {
    vector[16] dydt;
    real B[16,16];
    int count;
    if (theta[3]<0){
      count = 0; //real either first or last part of the B matrix to give each model
    } else {
      count = 256;
    }
    for (j in 1:16){
      for (i in 1:16){
        count = count + 1;
        B[i,j] = x_r[count]; //must be a better way to do this!!
      }
    }
    dydt = theta[1] * to_matrix(B) * to_vector(y) + theta[2];
    return to_array_1d(dydt);
  }
}
data {
  int<lower=1> T;
  real y[T,16];
  real y0[16];
  real t0;
  real ts[T];
  real bothB[16*16*2]; // unidirectional RC connectivity matrix
}
transformed data {
  real x_r[16*16*2];
  int x_i[0];
}
parameters {
  real<lower=0> sigma;
  real theta[3];
}
model {
  real z[T,16];
  sigma ~ cauchy(0,3.5);
  theta[1] ~ normal(0,2);
  theta[2] ~ normal(0,2);
  theta[3] ~ normal(0,0.2);
  z <- integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, bothB, x_i);
  for (t in 1:T){
    y[t] ~ normal(z[t], sigma);
  }
}
generated quantities {
  real y_pred[T,16];
  y_pred = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, bothB, x_i);
  for (t in 1:T){
    for (j in 1:16){
      y_pred[t,j] = y_pred[t,j] + normal_rng(0,sigma);
    }
  }
}
//only relies on the oocyte observation, not all of the nurse cells
