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
 
 for (i in 2:(nRe)){
   gprimeprime=gprimeprime+(i^2)*(phi^2)*exp(-i*phi*smlgamma)/((1-exp(-i*phi*smlgamma))^2)
 }
 
 for (i in 2:(nRe)){
   gprimeprime=gprimeprime-1*(phi^2)*exp(-phi*smlgamma)/((1-exp(-phi*smlgamma))^2)
 }
 
  Cgamma=sum(c(1:(nRe+1))^(-smlgamma)) # C(.) normalizing constant
  normi1=0
  normi=sum(c(1:(nRe+1))^(-smlgamma)*(log(c(1:(nRe+1)))^2))
  normi2=0
  normi2=(sum(c(1:(nRe+1))^(-smlgamma)*log(c(1:(nRe+1)))))^2
  gprimeprime=gprimeprime-(n-nRe)/((Cgamma)^2)*(Cgamma*normi1+normi2)
  return(gprimeprime)
}


