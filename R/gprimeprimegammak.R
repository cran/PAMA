gprimeprimegammak=function(bsrkr,I,phi,smlgamma){
  # this function returns log-likelihood of observation
  # bsrkr is the observed base ranker
  # I is the true classification of entities
  # I_r is the true ranking of the relative entities
  # phi is the disperse parameter in Mallows model
  # smlgamma is the parameter to distinguish relative and background entities
  gprimeprime=0
  nRe=sum(I>0)
 n=length(I)
  
  for (i in 1:(nRe-1)){
    denorm=0
    normi1=0
    normi2=0
    for (j in 0:i){
      denorm=denorm+exp(-j*phi*smlgamma)
      normi1=normi1+exp(-j*phi*smlgamma)*j*phi
      normi2=normi2+exp(-j*phi*smlgamma)*j^2*phi^2
    }
    gprimeprime=gprimeprime+(normi1^2-denorm*normi2)/denorm^2
  }
  Cgamma=sum(c(1:(nRe+1))^(-smlgamma)) # C(.) normalizing constant
  normi1=0
  normi=sum(c(1:(nRe+1))^(-smlgamma)*(log(c(1:(nRe+1)))^2))
  normi2=0
  normi2=(sum(c(1:(nRe+1))^(-smlgamma)*log(c(1:(nRe+1)))))^2
  gprimeprime=gprimeprime-(n-nRe)/((Cgamma)^2)*(Cgamma*normi1+normi2)
  return(gprimeprime)
}


