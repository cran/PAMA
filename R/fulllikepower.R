fulllikepower= function(dat,I,phi,smlgamma){
  # this returns log-likelihood of data
  #source('BARDMallowslikepower.R')
  output=lapply(seq_len(ncol(dat)), function(i) BARDMallowslikepower(dat[,i],I,phi,smlgamma[i]))
  a=unlist(output)
  return(sum(a))# log likelihood
}
