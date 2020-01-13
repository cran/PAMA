gprimeprimephi=function(bsrkr,I,phi,smlgamma){
  # this function returns log-likelihood of observation
  # bsrkr is the observed base ranker
  # I is the true classification of entities
  # I_r is the true ranking of the relative entities
  # phi is the disperse parameter in Mallows model
  # smlgamma is the parameter to distinguish relative and background entities
  m=length(smlgamma)
  nRe=sum(I>0)
  gprimeprime=0
  for(k in 1:m){
    for (i in 1:(nRe-1)){
      denorm=0
      normi1=0
      normi2=0
      for (j in 0:i){
        denorm=denorm+exp(-j*phi*smlgamma[k])
        normi1=normi1+exp(-j*phi*smlgamma[k])*j*smlgamma[k]
        normi2=normi2+exp(-j*phi*smlgamma[k])*j^2*smlgamma[k]^2
      }
      gprimeprime=gprimeprime+(normi1^2-denorm*normi2)/denorm^2
    }
  }
  return(gprimeprime)
}


