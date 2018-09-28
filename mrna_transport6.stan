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
    producers = rep_vector(1,16);
    producers[1] = 0;
    B = to_matrix(x_r,16,16);
    dydt = theta[1] * B * to_vector(y) + theta[2] * producers;
    return to_array_1d(dydt);
  }
}
data {
  int<lower=1> T;
  real y0[16];
  real t0;
  real ts[T];
  real theta[2];
  real sigma;
  real phi;
  real B[16*16];
}
transformed data {
  real x_r[16*16];
  int x_i[0];
}
parameters {
}
model {
}
generated quantities {
  int y_hat[T,16];
  real y_ode[T,16];
  y_ode = integrate_ode_rk45(mrnatransport, y0, 0, ts, theta, B, x_i );
  for (t in 1:T){
    for (j in 1:16){
      if (j>1){
        y_hat[t,j] = neg_binomial_2_rng(y_ode[t,j], sigma);
      } else {
        y_hat[t,j] = neg_binomial_2_rng(phi*y_ode[t,j], sigma);    
      }
    }
  }
}
