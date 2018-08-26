get_producers <- function(nTestOE){
  producers = matrix(data=NA,nrow=nTestOE,ncol=16)
  for (i in seq_len(nTestOE)){
    path = paste('data/Overexpression',i,'/binary.csv',sep='')
    if (file.exists(path)){
      producers[i,] = read.csv(path,head=FALSE)[,2] #read manually determined data
    } else {
      producers[i,] = rep(2,16) #use default
      producers[i,1] = 0
    }
  }  
  return(producers)
}
