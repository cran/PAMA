gprimegammak=function(bsrkr,I,phi,smlgamma){
  # this function returns log-likelihood of observation
  # bsrkr is the observed base ranker
  # I is the true classification of entities
  # I_r is the true ranking of the relative entities
  # phi is the disperse parameter in Mallows model
  # smlgamma is the parameter to distinguish relative and background entities. smlgammak.
  gprime=0
  rank_RE=bsrkr[I>0] #find out the relative entities
  nRe=length(rank_RE) # d
  n=length(I)
  gprime=gprime+(-phi)*distance(rank(rank_RE),I[I>0])
  
  
  tau01=conditionalranking(I,bsrkr)+1 # return tau01 (power law 1:nRe+1)
  n01=unlist(lapply(c(1:(nRe+1)), function(i) sum(tau01==i)))
  gprime=gprime-sum(n01*log(c(1:(nRe+1)))) # part2
  
  #E
  for(i in 2:nRe){
    gprime = gprime- i*phi*exp(-i*phi*smlgamma)/(1-exp(-i*phi*smlgamma))
  }
  
  for(i in 2:nRe){
    gprime = gprime+ 1*phi*exp(-phi*smlgamma)/(1-exp(-phi*smlgamma))
  }
  
  
  
  
  Cgamma=sum(c(1:(nRe+1))^(-smlgamma)) # C(.) normalizing constant
  normi=0
  normi=sum(c(1:(nRe+1))^(-smlgamma)*log(c(1:(nRe+1))))
  gprime=gprime+(n-nRe)*normi/Cgamma
  return(gprime)
}


