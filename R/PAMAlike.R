PAMAlike=function(bsrkr,I,phi,smlgamma){
  #' this function returns the log-likelihood of an observed ranking list given the parameters.
  #' @export
  #' @export
  #' @import PerMallows
  #' @import stats
  #' @import mc2d
  #' @import ExtMallows
  #' @param bsrkr It is a observed ranking list.
  #' @param I It is the true classification of entities. 0 denotes the corresponding entity is a background entity. The positive integer denotes the relative rankings of a relevant entity.
  #' @param phi It is a positive number. It is the common disperse parameter in Mallows model for all the rankers
  #' @param  smlgamma A positive number. It is the quality parameter of the ranker. It is used to distinguish relative and background entities
  #' @return The lon-likelihood of barkr given I, phi and smlgamma
  #' @examples
  #' dat=t(PerMallows::rmm(10,1:20,0.5))
  #' I=c(1:10,rep(0,10))
  #' like=PAMAlike(bsrkr=dat[,1],I=I,phi=0.2,smlgamma=1)
  # source('conditionalranking.R')

  rank_RE= bsrkr[I>0] #find out the relative entities
  mallowlike=PerMallows::dmm(rank(rank_RE),I[I>0], phi*smlgamma) # likelihood of reletive entities using Mallows

  log.mallowlike = log(mallowlike)
  tau01=conditionalranking(I,bsrkr)+1 # return tau01 (power law 1:nRe+1)
  possiblepos=c(1:(length(rank_RE)+1))
  C_gamma=sum((possiblepos^(-smlgamma))) # normalizing constant
  log.tau01like=sum( log(tau01^(-smlgamma))-log(C_gamma) ) #log.likelihood of tau01
  logfactorial= function(n){
    if(n>0&&n<=100){
      s=log(prod(c(1:n)))
    }else if(100<n && n<=200){
      s=log(prod(c(1:100)))+log(prod(c(101:n)))
    }else if(200<n){s=log(prod(c(1:100)))+log(prod(c(101:200)))+log(prod(c(201:n)))}
    else if (n==0){s=0}
    return(s)
  }

   possiblecomb=lapply(c(1:(length(rank_RE)+1)), function(i) logfactorial(sum(tau01==i))) # calculate the combinations of background entities

    log.backgroundlike=-sum((unlist(possiblecomb)))
  fulllike= log.mallowlike +log.tau01like  +log.backgroundlike
  return(fulllike)
}


