functions {
  matrix qr_Q_stan(matrix A) {
    return qr_Q(A); // qr.Q(qr, complete = TRUE) in R
  }
  matrix qr_R_stan(matrix A) {
    return qr_R(A); // qr.R(qr, complete = TRUE) in R
  }
  matrix lkj_cor_rng(int K, real eta) {
    return lkj_corr_rng(K, eta);
  }
}
