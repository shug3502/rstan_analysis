functions {
  matrix construct_matrix(real th) {
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
  return B';
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
  real get_first_nonzero_entry(vector v){
    int z;
    z = 1;
    while (fabs(v[z])<10^-15){
      z = z+1;
      if (z>16){
        return(0.0);
      }
    }
    return(v[z]);
  }

 vector get_k2(real th, int[] x_i){
   //alter matrix for blocking and compute QR decomposition from this
    matrix[16,16] B;
    matrix[16,16] Q;
    vector[16] N;
    vector[16] N_tilde;
    vector[16] N_bar;
    real s;
    real nz;
    real nz_tilde;
    int y_i[2];
    
    if (sum(x_i) > 0){
      // then we will block a ring canal
      B = alter_matrix(construct_matrix(th), x_i); 
      Q = qr_Q(B'); //compute qr decomposition
    //take last n-r cols of Q as basis of null space
    //here null spcae is of dimension 2 when we have altered the matrix
      N = Q[1:16,16];
      s = get_first_nonzero_entry(N);
      for (i in 1:16){ //normalise
          N[i] = N[i]/s;
      }
      N_tilde = Q[1:16,15];
      s = get_first_nonzero_entry(N_tilde);
      for (i in 1:16){ //normalise
          N_tilde[i] = N_tilde[i]/s;
      }
      //need to combine and weight by number of nurse cells
      nz=0;
      nz_tilde=0;
      for (j in 1:16){
        if (N[j] > 10^-15){
          nz = nz+1;
        }
        if (N_tilde[j] > 10^-15){
          nz_tilde = nz_tilde + 1;
        }
      }
      if (nz>nz_tilde){
        nz = nz-1; //account for ooxyte which does not produce
      } else {
        nz_tilde = nz_tilde-1;
      }
      //combine and weight
      N_bar = nz/(nz+nz_tilde)*N + nz_tilde/(nz+nz_tilde)*N_tilde; 
    } else {
    
      // no blocking
      y_i[1]=15;
      y_i[2]=16;
      B = alter_matrix(construct_matrix(th), y_i); 
      Q = qr_Q(B'); //compute qr decomposition
      N = Q[1:16,16];
      s = N[1];
      for (i in 1:16){ //normalise
          N[i] = N[i]/s;
      }
  N_bar = Q[1:16,16];
  s = N_bar[1];
  for (i in 1:16){ //normalise
      N_bar[i] = N_bar[i]/s;
  }
 }
  return(N_bar);
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
}
data {
  int<lower=1> T1;
  int<lower=0> T2;
  real y_obs[T1,16];
  //int x_i[2];
}
parameters {
  simplex[16] theta; //probability of blocking for each ring canal
  real<lower=0,upper=1> nu;
  real<lower=0> xi;
  real<lower=0> phi;
}
model {
  real y_stst[T1,16];
  vector[16] alpha;
  int z[T1];
  int blocked_cells[2];
  alpha = rep_vector(1,16);
  print(alpha);
  theta ~ dirichlet(alpha);
  print(theta);
  xi ~ normal(0,0.1) T[0,]; 
  nu ~ beta(1,1) T[0,1];
  phi ~ normal(0.30,0.036) T[0,1]; 
  for (t in 1:T1){
    z[t] ~ categorical(theta);
    blocked_cells = get_RC_from_dict(z[t]);
    y_stst[t] = to_array_1d(get_k2(nu,blocked_cells));
    for (j in 2:16){
      y_obs[t,j] ~ normal(y_stst[t,j]/phi,xi) T[0,];
    }
  }
}
generated quantities {
  real y_pred[(T1+T2),16];
  real y_sim[(T1+T2),16];
  vector[T1] log_lik;
  real y_stst[T1,16];
  int z_sim[T1];
  int blocked_cells_sim[2];
  for (t in 1:(T1+T2)) {
    z_sim[t] = categorical_rng(theta);    
    blocked_cells_sim = get_RC_from_dict(z_sim[t]);
    y_pred[t] = to_array_1d(get_k2(nu,blocked_cells_sim));
    y_sim[t,1] = 1; 
    for (j in 2:16){
      y_sim[t,j] = fabs(normal_rng(y_pred[t,j]/phi, xi));
    }
  }
    //compute log likelihood for model comparison via loo
  log_lik = rep_vector(0,T1);
  for (t in 1:T1){
    y_stst[t] = to_array_1d(get_k2(nu,blocked_cells_sim));
  }
}
