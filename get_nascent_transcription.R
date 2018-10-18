get_nascent_transcription <- function(nSamples,phenotype='OE'){

  root=dplyr::case_when(phenotype=='OE' ~ 'Overexpression',
			phenotype=='WT' ~ 'Example',
			TRUE ~ '') #TODO: better error handling
  nascent = matrix(data=NA,nrow=nSamples,ncol=16)
  for (i in seq_len(nSamples)){
    path = paste('data/',root,i,'/nascent_transcription.csv',sep='')
    if (file.exists(path)){
      nascent[i,] = read.csv(path,head=TRUE)[,2] #read manually determined data
    } else {
      nascent[i,] = rep(1,16) #use default
      nascent[i,1] = 0
    }
  }
    return(nascent)
}
