fulllikepower= function(dat,I,phi,smlgamma){
  # this returns log-likelihood of data
  #source('PAMAlike.R')
  output=lapply(seq_len(ncol(dat)), function(i) PAMAlike(dat[,i],I,phi,smlgamma[i]))
  a=unlist(output)
  return(sum(a))# log likelihood
}
