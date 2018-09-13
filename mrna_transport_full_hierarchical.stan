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
  x_r[1] = 0;
  for (j in 2:16){
    x_r[j] = 1; //in WT,assume all nurse cells produce RNA equally, but none from oocyte
  }
}
parameters {
  real<lower=0> sigma; //noise param
  real<lower=0,upper=1> phi; //difference between particles in NCs and in Oocyte
  real mu[3]; //use this as popn mean
  real<lower=0> tau[3];
  real transformed_theta[T1,3];
}

transformed parameters{
  real theta[T1,4];
  real a; //want these for ease of plotting
  real b;
  real nu;
  a = mu[2];
  b = mu[1];
  nu = (1+inv_logit(mu[3]))/2;
  for (i in 1:T1){
    for (j in 1:2){
      theta[i,j] = fabs(mu[j] + tau[j]*transformed_theta[i,j]); //truncated normal
    }
    theta[i,3] = 0;
    theta[i,4] = inv_logit(mu[3] + tau[3]*transformed_theta[i,3]);
  }
}

model {
  real z[T1,16];
  real aux[T1,16];
  sigma ~ normal(0,10) T[0,]; 
  phi ~ normal(0.57,0.0118) T[0,1];
  for (i in 1:2){
    mu[i] ~ normal(0,10) T[0,];    
  }
  mu[3] ~ normal(2,2);
  for (i in 1:3){
    tau[i] ~ normal(0,1) T[0,];
  }
  for (t in 1:T1){
    transformed_theta[t,1] ~ normal(0,1);
    transformed_theta[t,2] ~ normal(0,1);    
//    transformed_theta[t,3] ~ normal(0,1);
    transformed_theta[t,3] ~ normal(0,1); //logit for nu   
    aux = integrate_ode_rk45(mrnatransport, y0, 0, ts1, theta[t], x_r, x_i);
    z[t] = aux[t];
    for (j in 1:16) {
      if (j>1){
        y[t,j] ~ neg_binomial_2(z[t,j], sigma);
      } else {
        y[t,j] ~ neg_binomial_2(phi*z[t,j], sigma);  //treat observations in oocyte differently due to aggregation of rna complexes  
      }
    }
  }
}
generated quantities {
  int y_pred[T2,16]; //predictions for WT
  real y_ode[T2,16]; 
  real theta_wt[T2,4]; //parameters for each WT egg chamber
  int y_pred_OE[T3,16]; //predictions for OE
  real y_ode_OE[T3,16];  
  real theta_OE[T3,4]; //parameters for each OE egg chamber
  real y_lik_ode[T1,16];
  vector[T1] log_lik;
  real aux_gq[max(max(T1,T2),T3),16];
  // for wild type
  for (t in 1:T2){
    for (j in 1:2){
      theta_wt[t,j] = fabs(normal_rng(mu[j],tau[j]));
    }
    theta_wt[t,3] = 0;
    theta_wt[t,4] = inv_logit(normal_rng(mu[3],tau[3]));    
    aux_gq[1:T2] = integrate_ode_rk45(mrnatransport, y0, 0, ts2, theta_wt[t], x_r, x_i);
    y_ode[t] = aux_gq[t];
    for (j in 1:16){
      if (j>1){
        y_pred[t,j] = neg_binomial_2_rng(y_ode[t,j], sigma);
      } else {
        y_pred[t,j] = neg_binomial_2_rng(phi*y_ode[t,j], sigma);    
      }
    }
  }
  // for the overexpression mutant
  for (t in 1:T3){
    theta_OE[t,1] = fabs(normal_rng(mu[1],tau[1]));
    theta_OE[t,2] = 2*fabs(normal_rng(mu[2],tau[2]));    
//    theta_OE[t,3] = fabs(normal_rng(mu[3],tau[3]));
    theta_OE[t,3] = 0;
    theta_OE[t,4] = inv_logit(normal_rng(mu[3],tau[3]));
    aux_gq[1:T3] = integrate_ode_rk45(mrnatransport, y0, 0, ts3, theta_OE[t], to_array_1d(OE_producers[t,]), x_i);
    y_ode_OE[t] = aux_gq[t];
    for (j in 1:16){
      if (j>1){
        y_pred_OE[t,j] = neg_binomial_2_rng(y_ode_OE[t,j], sigma);
      } else {
        y_pred_OE[t,j] = neg_binomial_2_rng(phi*y_ode_OE[t,j], sigma);    
      }
    }
  }
    //compute log likelihood for model comparison via loo
  log_lik = rep_vector(0,T1);
  for (t in 1:T1){
    aux_gq[1:T1] = integrate_ode_rk45(mrnatransport, y0, 0, ts1, theta[t], x_r, x_i );
    y_lik_ode[t] = aux_gq[t];
    for (j in 1:16){
      if (j>1) {
        log_lik[t] = log_lik[t] + neg_binomial_2_lpmf(y[t,j] | y_lik_ode[t,j],sigma);
      } else {
        log_lik[t] = log_lik[t] + neg_binomial_2_lpmf(y[t,j] | y_lik_ode[t,j]*phi,sigma);
      }
    }
  }
}
 

 
