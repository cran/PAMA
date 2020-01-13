gprimephi=function(bsrkr,I,phi,smlgamma){
  # this function returns log-likelihood of observation
  # bsrkr is the observed base ranker. This should be a matrix
  # I is the true classification of entities
  # I_r is the true ranking of the relative entities
  # phi is the disperse parameter in Mallows model
  # smlgamma is the parameter to distinguish relative and background entities

  n=length(I)
  m=length(smlgamma)
  gprime=0
  for(k in 1:m){
    rank_RE= rank(bsrkr[I>0,k]) #find out the relative entities and calculate the relative rank

    gprime=gprime - sum(smlgamma*PerMallows::distance(rank_RE,I[I>0]))
  }
  nRe=length(rank_RE)
  #tau01=conditionalranking(I,bsrkr) # return tau01 (power law 1:nRe+1)

  for (k in 1:m){
    for (i in 1:(nRe-1)){
      denorm=0
      normi=0
      for (j in 0:i){
        denorm=denorm+exp(-j*phi*smlgamma[k])
        normi=normi+exp(-j*phi*smlgamma[k])*j*smlgamma[k]
      }
      gprime=gprime+normi/denorm
    }
  }

  return(gprime)
}


