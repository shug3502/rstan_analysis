get_blocked_indices <- function(nTestOE){
  blocked = matrix(data=NA,nrow=nTestOE,ncol=3)
  for (i in seq_len(nTestOE)){
    path = paste('data/Overexpression',i,'/blocked.csv',sep='')
    if (file.exists(path)){
      blocked[i,] = as.integer(read.csv(path,head=FALSE)[1,]) #read manually determined data
    } else {
      warning('path does not exist for blocking info')
      blocked[i,] = as.integer(rep(1,3)) #use default
    }
  }  
  return(blocked)
}

