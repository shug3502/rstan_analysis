
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
    B = construct_matrix(theta[4],theta[3]);
    dydt = theta[1] * B * to_vector(y) + theta[2] * to_vector(x_r);
    return to_array_1d(dydt);
  }
  real[] mrnatransport_altered(real t,
                       real[] y,
                       real[] theta,
                       real[] x_r,
                       int[] x_i
  ) {
    vector[16] dydt;
    matrix[16,16] B;
    B = construct_matrix(theta[4],theta[3]);
    // alter the matrix to remove a certain entry
    B[x_i[2],x_i[2]] = B[x_i[2],x_i[2]] + B[x_i[1],x_i[2]]; 
    B[x_i[1],x_i[1]] = B[x_i[1],x_i[1]] + B[x_i[2],x_i[1]]; 
    B[x_i[1],x_i[2]] = 0;
    B[x_i[2],x_i[1]] = 0;    
    dydt = theta[1] * B * to_vector(y) + theta[2] * to_vector(x_r);
    return to_array_1d(dydt);
  }
}
data {
  int<lower=1> T1;
  real y0[16];
  real t0;
  real ts1[T1]; //times for 'training data'
  real a;
  real b;
  real gamma;
  real nu;
  real sigma;
  real phi;
  int alter_matrix; //logical whether to alter matrix
  int b_i;
  int b_j;
}
transformed data {
  real x_r[16];
  int x_i[2];
  x_r[1] = 0;
  for (j in 2:16){
    x_r[j] = 1; //in WT,assume all nurse cells produce RNA equally, but none from oocyte
  }
  x_i[1] = b_i;
  x_i[2] = b_j;
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
  
  // for wild type
  theta[1]=b;
  theta[2]=a;
  theta[3]=gamma;
  theta[4]=nu;
  y_ode = integrate_ode_rk45(mrnatransport, y0, 0, ts1, theta, x_r, x_i);
  for (t in 1:T1){
    for (j in 1:16) {
      if (j>1){
        y_pred[t,j] = neg_binomial_2_rng(y_ode[t,j], sigma);
      } else {
        y_pred[t,j] = neg_binomial_2_rng(y_ode[t,j]*phi, sigma);  //treat observations in oocyte differently due to aggregation of rna complexes  
      }
    }
  }
}
 

 
