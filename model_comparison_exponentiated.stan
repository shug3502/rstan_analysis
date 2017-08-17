functions {
  real[,] mrna_transport_model(real t0, real[] ts, int T, real[] y0, matrix B, real[] theta){
    //use matrix exponentiation rather than solving an ode
    matrix[17, 17] K; //append constant to B matrix
    vector[16] unit_col;
    real x[T,17];
    unit_col = rep_vector(1,16);
    unit_col[1] = 0;
    K = append_row(append_col(theta[1]*B,theta[2]*unit_col),rep_row_vector(0,17));
    for (t in 1:T){
      x[t] = to_array_1d(matrix_exp((ts[t] - t0) * K) * to_vector(y0));
    }
    return x;
    }
}
data {
  int<lower=1> T;
  real y[T,17];
  real y0[17];
  real t0;
  real ts[T];
  matrix[16,16] biB;
  matrix[16,16] uniB;
}
parameters {
  real<lower=0> sigma;
  real theta[3];
}
transformed parameters {
  matrix[16,16] B; 
  if (theta[3]<0){
    B = biB;
  } else{
    B = uniB;
  }    
}
model {
  real z[T,17];
  sigma ~ cauchy(0,2.5);
  theta[1] ~ normal(0,2) T[0,];
  theta[2] ~ normal(0,2) T[0,];
  theta[3] ~ normal(0,0.2);
  z = mrna_transport_model(t0, ts, T, y0, B, theta);
  for (t in 1:T){
    y[t] ~ normal(z[t], sigma);
  }
}
generated quantities {
  real y_pred[T,17];
  y_pred =  mrna_transport_model(t0, ts, T, y0, B, theta);
  for (t in 1:T){
    for (j in 1:17){
      y_pred[t,j] = y_pred[t,j] + normal_rng(0,sigma);
    }
  }
}

