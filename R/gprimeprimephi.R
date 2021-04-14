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
    for (i in 2:(nRe)){
      gprimeprime=gprimeprime+(i^2)*(smlgamma[k]^2)*exp(-i*phi*smlgamma[k])/((1-exp(-i*phi*smlgamma[k]))^2)
     
      
    }
  }
  
  for(k in 1:m){
    for (i in 2:(nRe)){
      gprimeprime=gprimeprime-1*(smlgamma[k]^2)*exp(-phi*smlgamma[k])/((1-exp(-phi*smlgamma[k]))^2)
      
      
    }
  }
  return(gprimeprime)
}


