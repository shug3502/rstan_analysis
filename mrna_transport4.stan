functions {
  real[,] getmatrix(int N, real a){
    real d[N,N];
    int v[16];
    v[1] = 0;
    v[2] = 1;
    v[3] = 1;
    v[4] = 2;
    v[5] = 1;
    v[6] = 2;
    v[7] = 3;
    v[8] = 4;
    v[9] = 1;
    v[10] = 2;
    v[11] = 3;
    v[12] = 4;
    v[13] = 5;
    v[14] = 6;
    v[15] = 7;
    v[16] = 8;
    for (j in 1:N){
      for (i in 1:N){
        d[i,j] = 0;
      }
      if (j>1){
        d[v[j],j] = 1;
        d[j,j] = -1;
      }
    }
    return d;
  }
  real[] mrnatransport(real t,
                       real[] y,
                       real[] theta,
                       real[] x_r,
                       int[] x_i
  ) {
    vector[16] dydt;
    real B[16,16];
    B = getmatrix(16,0);
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
  /*
    real y0[16];
  //real z[T];
  real t0;
  real ts[T];
  real theta[2];
  real sigma[16];
  */
}
transformed data {
  real x_r[0];
  int x_i[0];
}
parameters {
  //  real y0[16];
  real<lower=0> sigma;
  real<lower=0> theta[2];
}
model {
  real z[T,16];
  sigma ~ cauchy(0,2.5);
  theta ~ normal(0,2);
  //  y0 ~ normal(0,1);
  z <- integrate_ode(mrnatransport, y0, t0, ts, theta, x_r, x_i);
  for (t in 1:T)
    y[t] ~ normal(z[t], sigma);
}
generated quantities {
  real y_pred[T,16];
  y_pred = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, x_r, x_i );
  for (t in 1:T){
    for (j in 1:16){
      y_pred[t,j] = y_pred[t,j] + normal_rng(0,sigma);
    }
  }
}

