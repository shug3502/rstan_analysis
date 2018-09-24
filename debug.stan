functions {
  matrix qr_Q_stan(matrix A) {
    return qr_Q(A); //reveal QR decomposition from Stan
  }
  matrix qr_R_stan(matrix A) {
    return qr_R(A); //reveal QR decomposition from Stan
  }
  int find_zero(vector v, int l){
    // integer l is length of v
    // will assume the final column will be a zero and want to find the other
    int k = 0;
    int j = 0;
    real tol = 10^-14;
    while (j < l){
      j += 1;
      if (fabs(v[j])<tol){
        k=j;
        break;
      }   
    }
    return(k);
  }
  int[] get_index_for_Q_cols(matrix A){
    int zero_indices[2];
    zero_indices[1] = find_zero(diagonal(qr_R(A)), 16);
    zero_indices[2] = 16;
    return(zero_indices); 
  } 
  //the following is to reproduce the matrix that I have been having problems with
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
}
