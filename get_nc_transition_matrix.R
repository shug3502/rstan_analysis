#easier to write functions to get transition matrices in R and feed in as data
get_nc_transition_matrix <- function(isBidirectional){
  if (isBidirectional){
    v <- list(c(2,3,5,9),c(1,4,6,10),c(1,7,11),c(2,8,12),c(1,13),c(2,14),c(3,15),c(4,16),1,2,3,4,5,6,7,8)
  } else {
    v <- c(0,1,1,2,1,2,3,4,1,2,3,4,5,6,7,8)
  }
  B = matrix(rep(0,16^2),nrow=16)
  for (j in 1:16){
    B[v[[j]],j] = 1;
    B[j,j] = ifelse(isBidirectional,-length(v[[j]]),ifelse(j>1,-1,0));
  }
  if (isBidirectional){
    B[,1] = rep(0,16)
  }
  return(B)
}