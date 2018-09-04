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
                       real[] u,
                       real[] theta,
                       real[] x_r,
                       int[] x_i
  ) {
    real dydt[16];
    real dudt[16];
    matrix[16,16] B;
    B = construct_matrix(theta[4],theta[3]);
    dydt = to_array_1d(theta[1] * B * exp(to_vector(u)) + theta[2] * to_vector(x_r));
    for (i in 1:16){
      dudt[i] = dydt[i]*exp(-u[i]);
    }
    return dudt;
  }
}
data {
  int<lower=1> T1;
  int<lower=1> T2;
  int<lower=1> T3;
  int y[T1,16];
  real y0[16];
  real t0;
  real ts1[T1]; //times for 'training data'
  real ts2[T2];  //times for 'test' data
  real ts3[T3]; //times for 'OE' test predictions
  matrix[T3,16] OE_producers; //matrix of how much RNA each cell produces in OE mutant, for predictions only
}
transformed data {
  real x_r[16];
  int x_i[0];
  real y0_trans[16];
  x_r[1] = 0;
  for (j in 2:16){
    x_r[j] = 1; //in WT,assume all nurse cells produce RNA equally, but none from oocyte
  }
  for (i in 1:16){
    y0_trans[i] = -4;
  }
}
parameters {
  real q;
}
model {
  2 ~ normal(q,1);
}
generated quantities {
  int y_pred[T1,16]; //predictions for WT
  real y_ode[T1,16]; 
  real theta[4]; //parameters for each WT egg chamber
  real b;
  real a;
  real gamma;
  real nu;
  real sigma;
  real phi;
  
  // for wild type
  
  b = fabs(normal_rng(0,1));
  a = fabs(normal_rng(0,100));
  gamma = fabs(normal_rng(0,0.01));
  nu = beta_rng(1,1);
  sigma = fabs(normal_rng(0,10)); 
  phi = fabs(normal_rng(0.57,0.118));
  theta[1]=b;
  theta[2]=a;
  theta[3]=gamma;
  theta[4]=nu;
  y_ode = integrate_ode_rk45(mrnatransport, y0_trans, t0, ts1, theta, x_r, x_i);
  for (t in 1:T1){
    for (j in 1:16) {
      if (j>1){
        y_pred[t,j] = neg_binomial_2_rng(exp(y_ode[t,j]), sigma);
      } else {
        y_pred[t,j] = neg_binomial_2_rng(exp(y_ode[t,j])*phi, sigma);  //treat observations in oocyte differently due to aggregation of rna complexes  
      }
    }
  }
}
 

 
