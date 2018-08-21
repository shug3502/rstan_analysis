functions {
  matrix construct_matrix(real th, real gamma_scaled) {
  //add decay at rate gamma
  real nu;
  matrix[16,16] B;
  nu = (1+th)/2;
  B = rep_matrix(0,16,16);
  B[1,1] = -4*(1-nu);
  B[2,1] = nu;
  B[3,1] = nu;
  B[5,1] = nu;
  B[9,1] = nu;
  B[1,2] = 1 - nu;
  B[2,2] = -3*(1-nu) - nu;
  B[4,2] = nu;
  B[6,2] = nu;
  B[10,2] = nu;
  B[1,3] = 1-nu;
  B[3,3] = -2*(1-nu) - nu;
  B[7,3] = nu;
  B[11,3] = nu;
  B[2,4] = 1-nu;
  B[4,4] = -2*(1-nu) - nu;
  B[8,4] = nu;
  B[12,4] = nu;
  B[1,5] = 1-nu;
  B[5,5] = -(1-nu) - nu;
  B[13,5] = nu;
  B[2,6] = 1-nu;
  B[6,6] = -(1-nu) - nu;
  B[14,6] = nu;
  B[3,7] = 1-nu;
  B[7,7] = -(1-nu) - nu;
  B[15,7] = nu;
  B[4,8] = 1-nu;
  B[8,8] = -(1-nu) - nu;
  B[16,8] = nu;
  B[1,9] = 1-nu;
  B[9,9] = -nu;
  B[2,10] = 1-nu;
  B[10,10] = -nu;
  B[3,11] = 1-nu;
  B[11,11] = -nu;
  B[4,12] = 1-nu;
  B[12,12] = -nu;
  B[5,13] = 1-nu;
  B[13,13] = -nu;
  B[6,14] = 1-nu;
  B[14,14] = -nu;
  B[7,15] = 1-nu;
  B[15,15] = -nu;
  B[8,16] = 1-nu;
  B[16,16] = -nu;
//need to return transpose
  return   B' - diag_matrix(rep_vector(gamma_scaled, 16));
  }

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
    B = construct_matrix(x_r[1],theta[3]);
    dydt = theta[1] * B * to_vector(y) + theta[2] * producers;
    return to_array_1d(dydt);
  }
}
data {
  int<lower=1> T1;
  int<lower=1> T2;
  int y[T1,16];
  real y0[16];
  real t0;
  real ts1[T1]; //times for 'training data'
  real ts2[T2];  //times for 'test' data
  real nu;
//  real gamma; //decay rate
}
transformed data {
  real x_r[1];
  int x_i[0];
  x_r[1]=nu;
//  x_r[2]=gamma;
}
parameters {
  real<lower=0> sigma; //noise param
  real<lower=0,upper=1> phi; //difference between particles in NCs and in Oocyte
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> gamma;
}
transformed parameters {
real theta[3];
theta[1] = b;
theta[2] = a;
theta[3] = gamma;
}
model {
  real z[T1,16];
  int cell_indices[16]; //only model these cells
  sigma ~ normal(0,10) T[0,]; //cauchy(0,2.5) T[0,]; //normal(1.0,0.25) T[0,]; 
  phi ~ normal(0.289,0.0285) T[0,1];
  a ~ normal(0,100) T[0,];
  b ~ normal(0,100) T[0,];
  gamma ~normal(0,1) T[0,];
  z = integrate_ode_rk45(mrnatransport, y0, t0, ts1, theta, x_r, x_i);
  cell_indices[1] = 1;
  cell_indices[2] = 2;
  cell_indices[3] = 3;
  cell_indices[4] = 9;
  cell_indices[5] = 5;
  for (t in 1:T1){
    for (j in 1:5) {
      //cell_indices[j]=j;
      if (j>1){
        y[t,cell_indices[j]] ~ neg_binomial_2(z[t,cell_indices[j]], sigma);
        //y[t,j] ~ poisson(z[t,cell_indices[j]]);
	      //y[t,cell_indices[j]] ~ normal(z[t,cell_indices[j]], sigma);
      } else {
        y[t,cell_indices[j]] ~ neg_binomial_2(phi*z[t,cell_indices[j]], sigma);    
        //y[t,j] ~ poisson(phi*z[t,cell_indices[j]]);    
	      //y[t,cell_indices[j]] ~ normal(phi*z[t,cell_indices[j]], sigma);
      }
    }
  }
}
generated quantities {
  int y_pred[T2,16];
  real y_ode[T2,16];
  real y_lik_ode[T1,16];
  vector[T1] log_lik;
  real theta_overexpression[3];
  int cell_indices[16]; //only model these cells
  for (j in 1:16){
    cell_indices[j] = j;
  }
  theta_overexpression[1] = theta[1];
  theta_overexpression[2] = 2*theta[2]; //double the rate of production
  theta_overexpression[3] = theta[3];
  y_ode = integrate_ode_rk45(mrnatransport, y0, t0, ts2, theta_overexpression, x_r, x_i );
  for (t in 1:T2){
    for (j in 1:16){
      if (j>1){
    	  y_pred[t,j] = neg_binomial_2_rng(y_ode[t,j],sigma); //normal_rng(y_ode[t,j], sigma);
      } else {
    	  y_pred[t,j] = neg_binomial_2_rng(phi*4.81/3.46*y_ode[t,j],sigma); //normal_rng(phi*4.81/3.46*y_ode[t,j], sigma);
      }
    }
  }
    //compute log likelihood for model comparison via loo
  log_lik = rep_vector(0,T1);
  y_lik_ode = integrate_ode_rk45(mrnatransport, y0, t0, ts1, theta, x_r, x_i );
  for (t in 1:T1){
    for (j in 1:5){
      if (j>1) {
        log_lik[t] = log_lik[t] + normal_lpdf(y[t,cell_indices[j]] | y_lik_ode[t,cell_indices[j]]/phi,sigma);
      } else {
        log_lik[t] = log_lik[t] + normal_lpdf(y[t,cell_indices[j]] | y_lik_ode[t,cell_indices[j]],sigma);
      }
    }
  }
}
 

 