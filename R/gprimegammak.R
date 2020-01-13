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
  gprime=gprime+(-phi)*PerMallows::distance(rank(rank_RE),I[I>0])


  tau01=conditionalranking(I,bsrkr)+1 # return tau01 (power law 1:nRe+1)
  n01=unlist(lapply(c(1:(nRe+1)), function(i) sum(tau01==i)))
  gprime=gprime-sum(n01*log(c(1:(nRe+1)))) # part2

  for (i in 1:(nRe-1)){
    denorm=0
    normi=0
    for (j in 0:i){
      denorm=denorm+exp(-j*phi*smlgamma)
      normi=normi+exp(-j*phi*smlgamma)*j*phi
    }
    gprime=gprime+normi/denorm
  }
  Cgamma=sum(c(1:(nRe+1))^(-smlgamma)) # C(.) normalizing constant
  normi=0
  normi=sum(c(1:(nRe+1))^(-smlgamma)*log(c(1:(nRe+1))))
  gprime=gprime+(n-nRe)*normi/Cgamma
  return(gprime)
}


