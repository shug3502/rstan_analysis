//runs inference with the mrna transport model without decay
//fits also the model for timescales together with this to forward propogate uncertainty about time points 
//forward simulates from a model with blocking of ring canals for the overexpressor 
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
  int[] get_RC_from_dict(int w){
    int RCs[32];
    int blocked_cells[2];
    RCs[1] = 0;
    RCs[2] = 0;
    RCs[3] = 1;
    RCs[4] = 2;
    RCs[5] = 1;
    RCs[6] = 3;
    RCs[7] = 1;
    RCs[8] = 5;
    RCs[9] = 1;
    RCs[10] = 9;
    RCs[11] = 2;
    RCs[12] = 4;
    RCs[13] = 2;
    RCs[14] = 6;
    RCs[15] = 2;
    RCs[16] = 10;
    RCs[17] = 3;
    RCs[18] = 7;
    RCs[19] = 3;
    RCs[20] = 11;
    RCs[21] = 4;
    RCs[22] = 8;
    RCs[23] = 4;
    RCs[24] = 12;
    RCs[25] = 5;
    RCs[26] = 13;
    RCs[27] = 6;
    RCs[28] = 14;
    RCs[29] = 7;
    RCs[30] = 15;
    RCs[31] = 8;
    RCs[32] = 16;
    blocked_cells[1] = RCs[2*w - 1];
    blocked_cells[2] = RCs[2*w];
    return(blocked_cells);
  }
  matrix alter_matrix(matrix A,
                      int[] x_i) {
     // alter the matrix to remove a certain entry
     //argument x_i are the ring canal indices
    matrix[16,16] B = A;
    B[x_i[2],x_i[2]] = B[x_i[2],x_i[2]] + B[x_i[1],x_i[2]];
    B[x_i[1],x_i[1]] = B[x_i[1],x_i[1]] + B[x_i[2],x_i[1]];
    B[x_i[1],x_i[2]] = 0;
    B[x_i[2],x_i[1]] = 0;
    return(B);
  }

  real[] mrnatransport(real t,
                       real[] y,
                       real[] theta,
                       real[] x_r,
                       int[] x_i
  ) {
    vector[16] dydt;
    matrix[16,16] B;
    int blocked_cells[2];
    B = construct_matrix(theta[3],0);
    for (j in 1:3){
      blocked_cells = get_RC_from_dict(x_i[j]);
      if (sum(blocked_cells)>0){ //otheriwse no alterations to matrix
        B = alter_matrix(B,blocked_cells);
      }
    }
    dydt = theta[1] * B * to_vector(y) + theta[2] * to_vector(x_r);
    return to_array_1d(dydt);
  }
  
    
  real[,] forward_simulate_OE_rng(int T,
                          real[] ts, 
                          real t0,
                          real[] y0,
                          real[] theta,
                          matrix producers,
                          int[,] OE_blocked){
//theta here is theta as above, but also augemented with phi and any other parameters
  real y_ode_OE[T,16];
  real y_pred_OE[T,16];
  int OE_x_i[3];
  real phi = theta[5];
  real sigma = theta[6];
  for (t in 1:T){
    for (i in 1:3){
      OE_x_i[i] = OE_blocked[t,i];
    }
    y_ode_OE = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, to_array_1d(producers[t,]), OE_x_i ); //use overexpression producers and specify any blocked RCs
    y_ode_OE[t,1] = y_ode_OE[t,1]*phi;
    for (j in 1:16){
      y_pred_OE[t,j] = neg_binomial_2_rng(y_ode_OE[t,j],sigma);
    }
  }
  return y_pred_OE;
}
  matrix my_log_lik(int[,] y_OE,
                    int T,
                    real[] ts, 
                    real t0,
                    real[] y0,
                    real[] theta,
                    matrix producers,
                    int[,] OE_blocked){
  matrix[T,16] ll = rep_matrix(0,T,16);
  real y_ode_OE[T,16];
  real phi = theta[5];
  real sigma = theta[6];
  int OE_x_i[3];  
  for (t in 1:T){
    for (i in 1:3){
      OE_x_i[i] = OE_blocked[t,i];
    }
    y_ode_OE = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, to_array_1d(producers[t,]), OE_x_i ); //use overexpression producers and specify which RCs are blocked
    y_ode_OE[t,1] = y_ode_OE[t,1]*phi;
    for (j in 1:16){
      ll[t,j] = ll[t,j] + neg_binomial_2_lpmf(y_OE[t,j] | y_ode_OE[t,j],sigma);
    }
  }
return ll;    
}
  
}
data {
  int<lower=1> T1;
  int<lower=1> T2;
  int<lower=1> T3;
  int y[T1,16];
  real y0[16];
  real t0; //intercept , log(a0) should be less than youngest egg chamber or will have to omit data
  real areas1[T1]; //times for 'training data'
  real areas2[T2];  //times for 'test' data
  real areas3[T3]; //times for 'OE' test predictions
  real obs_times1[T1]; //based on manual staging
  real obs_times2[T2];
  real obs_times3[T3];
  int OE_blocked[T3,3]; //matrix of indices for which ring canals are blocked, maximum of 3 blocked ring canals per egg chamber
  matrix[T3,16] OE_producers; //matrix of how much RNA each cell produces in OE mutant, for predictions only
  int y_OE[T3,16];
}
transformed data {
  real la1[T1] = sort_asc(log(areas1));
  real la2[T2] = sort_asc(log(areas2));
  real la3[T3] = sort_asc(log(areas3));  
  int ind1[T1] = sort_indices_asc(log(areas1));
  int ind2[T2] = sort_indices_asc(log(areas2));
  int ind3[T3] = sort_indices_asc(log(areas3));
  real x_r[16];
  int x_i[3];
  x_r[1] = 0;
  for (j in 2:16){
    x_r[j] = 1; //in WT,assume all nurse cells produce RNA equally, but none from oocyte
  }
  for (i in 1:3){
    x_i[i] = 1;
  }    
}
parameters {
  real<lower=0> sigma; //noise param
  real<lower=0> phi; //difference between particles in NCs and in Oocyte
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0,upper=1> nu;
  real<lower=0> tau;
  real<lower=0> epsilon;
}
transformed parameters {
  real theta[3];
  real ts1[T1]; //could write as matrix or array (note: not equal number of columns)
  real ts2[T2];
  real ts3[T3];
  theta[1] = b/tau; //due to scaled timescales
  theta[2] = a/tau; //due to scaled timescales
  theta[3] = nu;
  //assuming the same scaling for all times
  for (t in 1:T1){
    ts1[t] = t0 + log(areas1[t])*tau;
  }
  for (t in 1:T2){
    ts2[t] = t0 + log(areas2[t])*tau;
  }
  for (t in 1:T3){
    ts3[t] = t0 + log(areas3[t])*tau;
  }
}
model {
  real z[T1,16];
  tau ~ normal(0.06,0.02);
  epsilon ~ normal(0,1) T[0,];
  sigma ~ normal(0,10) T[0,]; 
  phi ~ normal(0.345,0.047) T[0,];
  a ~ normal(0,10) T[0,];
  b ~ normal(0,10) T[0,];
  nu ~ beta(1,1) T[0,1];
  //fit timescales model
  for (t in 1:T1){
    obs_times1 ~ normal(ts1[t],epsilon);
  }
  for (t in 1:T2){
    obs_times2 ~ normal(ts2[t],epsilon);
  }
  for (t in 1:T3){
    obs_times3 ~ normal(ts3[t],epsilon);
  }  
  z = integrate_ode_rk45(mrnatransport, y0, t0, la1, theta, x_r, x_i);
  for (t in 1:T1){
    for (j in 1:16) {
      if (j>1){
        y[ind1[t],j] ~ neg_binomial_2(z[t,j], sigma);
      } else {
        y[ind1[t],j] ~ neg_binomial_2(phi*z[t,j], sigma);  //treat observations in oocyte differently due to aggregation of rna complexes  
      }
    }
  }
}
generated quantities {
  int y_pred[T2,16]; //predictions for WT
  real y_ode[T2,16]; 
  int y_pred_OE[T3,16]; //predictions for OE
  real y_ode_OE[T3,16];  
  real theta_OE[3];
  real y_lik_ode[T3,16];
  int OE_x_i[3];
  vector[T3] log_lik;
  // for wild type
  y_ode = integrate_ode_rk45(mrnatransport, y0, t0, la2, theta, x_r, x_i );
  for (t in 1:T2){
    for (j in 1:16){
      if (j>1){
        y_pred[ind2[t],j] = neg_binomial_2_rng(y_ode[t,j], sigma);
      } else {
        y_pred[ind2[t],j] = neg_binomial_2_rng(phi*y_ode[t,j], sigma);    
      }
    }
  }
  // for the overexpression mutant
  theta_OE = theta;
  theta_OE[2] = theta[2]; //double the rate of production in the overexpressor, but do so via producers
  for (t in 1:T3){
    for (i in 1:3){
      OE_x_i[i] = OE_blocked[ind3[t],i]; //need to input as integers
    }
    y_ode_OE = integrate_ode_rk45(mrnatransport, y0, t0, la3, theta_OE, to_array_1d(OE_producers[ind3[t],]), OE_x_i); //use overexpression producers
    for (j in 1:16){
      if (j>1){
        y_pred_OE[ind3[t],j] = neg_binomial_2_rng(y_ode_OE[t,j], sigma);
      } else {
        y_pred_OE[ind3[t],j] = neg_binomial_2_rng(phi*y_ode_OE[t,j], sigma);    
      }
    }
  }

    //compute log likelihood for model comparison via loo
  log_lik = rep_vector(0,T3);
  for (t in 1:T3){
    for (i in 1:3){
      OE_x_i[i] = OE_blocked[ind3[t],i]; //need to input as integers
    }
    y_lik_ode = integrate_ode_rk45(mrnatransport, y0, t0, la3, theta_OE, to_array_1d(OE_producers[t,]), OE_x_i );
    for (j in 1:16){
      if (j>1) {
        log_lik[ind3[t]] = log_lik[ind3[t]] + neg_binomial_2_lpmf(y_OE[ind3[t],j] | y_lik_ode[t,j],sigma);
      } else {
        log_lik[ind3[t]] = log_lik[ind3[t]] + neg_binomial_2_lpmf(y_OE[ind3[t],j] | y_lik_ode[t,j]*phi,sigma);
      }
    }
  }
}
 

 
