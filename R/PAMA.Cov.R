
PAMA.Cov=function(datfile,Covdatfile,nRe,iter){
  #' This function implements Bayesian inference of PAMA model with covariates.
  #'
  #' @export
  #' @import PerMallows
  #' @import stats
  #' @import mc2d
  #' @param datfile A matrix or dataframe. This is the data where our algorithm will work on. Each colomn denotes a ranker's ranking. The data should be in entity-based format.
  #' @param Covdatfile  A matrix or dataframe. Each column denotes a covariate.
  #' @param nRe A number. Number of relevant entities
  #' @param iter A number. Numner of iterations of MCMC. Defaulted as 1000.
  #' @return List. It contains Bayesian posterior samples of all the parameters and log-likelihood.
  #' \enumerate{
  #'   \item I.mat: posterior samples of I
  #'   \item phi.mat: posterior samples of phi
  #'   \item smlgamma.mat: posterior samples of gamma
  #'   \item l.mat: posterior samples of log-likelihood.
  #'   \item theta.mat: posterior samples of coefficients of covariates.
  #' }
  #' @details The covariates are incoporated in the PAMA framework as indicators of groupmember. That is covariates are associated to group members via a logistic regression.
  #' @examples
  #' a=NBANFL()
  #' PAMA.Cov(t(a$NFLdata),a$NFLcov,nRe=10,iter=10)
  #' \donttest{PAMA.Cov(t(a$NFLdata),a$NFLcov,nRe=10,iter=1000)}

  # this function implements Bayesian inference of PAMA model with covariates.
  # The input
  # parameter 'datfile': this is the data where our algorithm will work on. Each row denotes a ranker's ranking. The data should be in entity-based format.
  # parameter 'Covdatfile': this is data of covariates.
  # parameter 'nRe' : number of relevant entities
  # parameter 'iter' : numner of iterations of MCMC
  #parameter 'threshold': the stopping threshold in determining convergence of MLE. if the two consecutive iterations of log-likelihood is smaller than threshold, then the convergence achives.

  # The output
  # Output: I.mat: posterior samples of I
  # Output: phi.mat: posterior samples of \phi
  # Output: smlgamma.mat: posterior samples of \gamma
  # Output: l.mat: posterior samples of log-likelihood
  #' @author Wanchuang Zhu, Yingkai Jiang, Jun S. Liu, Ke Deng
  #' @references Wanchuang Zhu, Yingkai Jiang, Jun S. Liu, Ke Deng (2021) Partition-Mallows Model and Its Inference for Rank Aggregation *Journal of the American Statistical Association*

  #source('conditionalranking.R')
  #source('BARDMallowslikepower.R')
  #source('fulllikepower.R')
  adaptation=iter*0.25
  dat=datfile
  n=dim(dat)[1]
  m=dim(dat)[2]
  covinfo=Covdatfile
  covinfofull=cbind(1,covinfo)
  nc=dim(covinfo)[2]

  gamma.hyper=rep(0.1,m)
  smlgamma.upper=10
  smlgamma.lower=0.1
  ## starting point
  phi.start=0.3
  phi.hyper=0.05
  I.start=PerMallows::lmm(t(dat))$mode
  I.start[I.start>nRe]=0
  theta.start=runif(nc+1,-1,1)
  theta.hyper=rep(0.01,nc+1)

  smlgamma.start=rep(1,m)
  I.mat=matrix(NA,n,iter)
  theta.mat=matrix(NA,nc+1,iter)

  I.r.list=list()
  smlgamma.mat=matrix(NA,m,iter)
  phi.mat=matrix(NA,iter,1)
  l.mat=matrix(NA,iter,1)
  for(i in 1:iter){
    ## update 0s in I
    zerolist=which(I.start==0)
    for(j in zerolist){
      nonzerolist=which(I.start==nRe)
      zeropvec=c() # for zeros, it is a (nRe+1)-dimensional multinomial distribution.

      I.new=replace(I.start,c(j,nonzerolist),I.start[c(nonzerolist,j)])
      zeropvec=c(zeropvec,fulllikepower(dat = dat,I = I.new,phi = phi.start,smlgamma = smlgamma.start) + theta.start%*%(covinfofull[j,]))

      zeropvec=c(zeropvec,fulllikepower(dat = dat,I = I.start,phi = phi.start,smlgamma = smlgamma.start)+ theta.start%*%(covinfofull[nonzerolist,] )) # keep the zero unchanged
      zeropvec=exp(zeropvec-(max(zeropvec)))
      zeropvec=zeropvec/sum(zeropvec)
      gibbsrlz=mc2d::rmultinomial(1,1,zeropvec) # multinomial sampling
      pos=which(gibbsrlz==1)

      if(pos==1){
        I.start=I.new
      }
    }
    ## update 1s I
    nonzerolist=which(I.start>0)
    for(j in (nRe):2){
      pos1=which(I.start==j)
      pos2=which(I.start==(j-1))
      I.new=replace(I.start,c(pos1,pos2),I.start[c(pos2,pos1)])

      nonzeropvec=rep(NA,2) # for nonzeros, it is a n-dimensional multinomial distribution

      nonzeropvec[1]=fulllikepower(dat = dat,I = I.new,phi = phi.start,smlgamma = smlgamma.start)
      nonzeropvec[2]=fulllikepower(dat = dat,I = I.start,phi = phi.start,smlgamma = smlgamma.start)

      nonzeropvec=exp(nonzeropvec-(max(nonzeropvec)))
      nonzeropvec=nonzeropvec/sum(nonzeropvec)
      gibbsrlz=mc2d::rmultinomial(1, 1,nonzeropvec) # multinomial sampling
      pos=which(gibbsrlz==1)

      if(pos==1){ #
        I.start=I.new
      }
    }

    ## update I_r and phi
    ## the posterior likelihood can be written out,
    Mallowsdat=dat[I.start>0,]
    Mallowsdat=apply(Mallowsdat,2,rank)

    ##################################
    phi.new=phi.start+ phi.hyper* rnorm(1)
    if (phi.new>0){
      log.prob.start <- lapply(seq_len(ncol(Mallowsdat)), function(i) log(dmm(Mallowsdat[,i],I.start[I.start>0], phi.start*smlgamma.start[i])) )
      # log.prob.start <- apply(Mallowsdat,2,FUN=function(x){log(dmm(x,I.start[I.start>0], -log(phi.start)))})
      # log.prob.new <- apply(Mallowsdat,2,FUN=function(x){log(dmm(x,I.start[I.start>0], -log(phi.new)))})
      log.prob.new <- lapply(seq_len(ncol(Mallowsdat)), function(i) log(dmm(Mallowsdat[,i],I.start[I.start>0], phi.new*smlgamma.start[i])) )
      if (sum(unlist(log.prob.new))-sum(unlist(log.prob.start)) >log(runif(1))){
        phi.start=phi.new
      }
    }

    # phi.start=exp(-lmm(t(Mallowsdat))$theta)
    ## update smallgamma. MH is in use #############
    for(j in 1:m){
      smlgamma.tem=smlgamma.start[j]
      smlgamma.tem=smlgamma.tem+gamma.hyper[j]*rnorm(1)
      if(smlgamma.tem>smlgamma.lower && smlgamma.tem< smlgamma.upper){
        like.start=BARDMallowslikepower(dat[,j],I.start,phi.start,smlgamma.start[j])
        like.tem=BARDMallowslikepower(dat[,j],I.start,phi.start,smlgamma.tem)
        if((like.tem-like.start) > log(runif(1))){
          smlgamma.start[j]=smlgamma.tem
        }
      }
    }
    ## update theta's
    for(j in 1:(nc+1)){
      theta.tem=theta.start[j]+theta.hyper[j]*rnorm(1)
      theta.new=theta.start
      theta.new[j]=theta.tem
      mmat.start=t(apply(covinfofull,1,function(x) theta.start*x))
      mmat.start=apply(mmat.start,2,sum)
      mmat.new=t(apply(covinfofull,1,function(x) theta.new*x))
      mmat.new=apply(mmat.new,2,sum)
      like.start=theta.start%*%(apply(covinfofull[I.start>0,], 2, sum))-sum(log(1+exp(mmat.start)))
      like.new=theta.new%*%((apply(covinfofull[I.start>0,], 2, sum)))-sum(log(1+exp(mmat.new)))
      if(like.new==-Inf){ ## sometimes exp will lead to Inf
        like.new=theta.new%*%((apply(covinfofull[I.start>0,], 2, sum)))-sum(mmat.new[mmat.new>0])
      }
      if(like.start==-Inf){
        like.start=theta.start%*%((apply(covinfofull[I.start>0,], 2, sum)))-sum(mmat.start[mmat.start>0])
      }
      if((like.new-like.start)> log(runif(1))){
        theta.start=theta.new
      }

    }

    I.mat[,i]=I.start
    phi.mat[i]=phi.start
    smlgamma.mat[,i]=smlgamma.start
    l.mat[i]=fulllikepower(dat = dat,I = I.start,phi = phi.start,smlgamma = smlgamma.start)
    theta.mat[,i]=theta.start
    if(i==adaptation){
      gamma.hyper=sqrt(diag(cov(t(smlgamma.mat[,(adaptation-adaptation*0.6):adaptation]))))
      phi.hyper=sd((phi.mat[(adaptation-adaptation*0.6):adaptation]))
      theta.hyper=sqrt(diag(cov(t(theta.mat[,(adaptation-adaptation*0.6):adaptation]))))
    }
  }
  return(list(I.mat=I.mat,phi.mat=phi.mat,smlgamma.mat=smlgamma.mat,l.mat=l.mat,theta.mat=theta.mat))
}
