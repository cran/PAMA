PAMA.PL=function(datfile,PLdatfile,nRe,iter){
  #' This function implements Bayesian inference of PAMA model with partial lists.
  #'
  #' @export
  #' @import PerMallows
  #' @import stats
  #' @import mc2d
  #' @param datfile A matrix or dataframe. This is the data where our algorithm will work on. Each colomn denotes a ranker's ranking. The data should be in entity-based format.
  #' @param PLdatfile  A matrix or dataframe.  It contains all the partial lists. Each colomn denotes a partial list.
  #' @param nRe A number. Number of relevant entities.
  #' @param iter A number. Numner of iterations of MCMC. Defaulted as 1000.
  #' @return List. It contains Bayesian posterior samples of all the parameters and log-likelihood.
  #' \enumerate{
  #'   \item I.mat: posterior samples of I
  #'   \item phi.mat: posterior samples of phi
  #'   \item smlgamma.mat: posterior samples of gamma
  #'   \item l.mat: posterior samples of log-likelihood.
  #' }
  #' @details The partial lists are handle by Data Augmentation strategy.
  #' @examples
  #' a=NBANFL()
  #' PAMA.PL(a$NBA,a$NBAPL,nRe=10,iter=1)
  #' \donttest{PAMA.PL(a$NBA,a$NBAPL,nRe=10,iter=100)}
  #'
  # this function implements Bayesian inference of PAMA model with Partial lists.
  # The input
  # parameter 'datfile': this is the data where our algorithm will work on. Each colomn denotes a ranker's ranking. The data should be in entity-based format.
  # parameter 'PLdatfile': this is data of partial lists.
  # parameter 'nRe' : number of relevant entities
  # parameter 'iter' : numner of iterations of MCMC
  #parameter 'threshold': the stopping threshold in determining convergence of MLE. if the two consecutive iterations of log-likelihood is smaller than threshold, then the convergence achives.
  #' @author Wanchuang Zhu, Yingkai Jiang, Jun S. Liu, Ke Deng
  #' @references Wanchuang Zhu, Yingkai Jiang, Jun S. Liu, Ke Deng (2021) Partition-Mallows Model and Its Inference for Rank Aggregation. Journal of the American Statistical Association

  # The output
  #  Output: I.mat: posterior samples of I
  # Output: phi.mat: posterior samples of \phi
  # Output: smlgamma.mat: posterior samples of \gamma
  # Output: l.mat: posterior samples of log-likelihood

  #source('conditionalranking.R')
  #source('BARDMallowslikepower.R')
  #source('fulllikepower.R')
  #source('pl2fullV2.R')


  adaptation=0.25*iter
  dat=datfile
  datPL=PLdatfile
  n=dim(dat)[1]
  m=dim(dat)[2]
  mPL=dim(PLdatfile)[2]
  ## the info about the data


  gamma.hyper=rep(0.15,m+mPL)
  smlgamma.upper=5
  smlgamma.lower=0.01
  ## starting point
  phi.start=0.3
  phi.hyper=0.1
  # startiing points of I
  mallowresult=PerMallows::lmm(t(dat),dist.name="kendall", estimation="approx")
  mallowresult=mallowresult$mode
  mallowresult[mallowresult>nRe]=0
  I.start=mallowresult

  smlgamma.start=c(rep(3,m),rep(0.3,mPL))
  I.mat=matrix(NA,n,iter)
  l.mat=matrix(NA,1,iter)
  I.r.list=list()
  smlgamma.mat=matrix(NA,m+mPL,iter)
  phi.mat=matrix(NA,iter,1)
  for(i in 1:iter){
    FdatPL=pl2fullV2(I.start,datPL,nRe,phi.start,smlgamma.start[-c(1:m)])
    zerolist=which(I.start==0)
    for(j in zerolist){
      nonzerolist=which(I.start==nRe) # positions of non-zeros
      zeropvec=c() # for zeros, it is a (nRe+1)-dimensional multinomial distribution.
      I.new=replace(I.start,c(j,nonzerolist),I.start[c(nonzerolist,j)])
      FdatPLnew=pl2fullV2(I.new,datPL,nRe,phi.start,smlgamma.start[-c(1:m)])
      zeropvec=c(zeropvec,fulllikepower(dat = cbind(dat,FdatPLnew),I = I.new,phi = phi.start,smlgamma = smlgamma.start))

      zeropvec=c(zeropvec,fulllikepower(dat = cbind(dat,FdatPL),I = I.start,phi = phi.start,smlgamma = smlgamma.start)) # keep the zero unchanged
      zeropvec=exp(zeropvec-(max(zeropvec)))# fulllike returns log-likelihood
      zeropvec=zeropvec/sum(zeropvec)

      gibbsrlz=mc2d::rmultinomial(1,1,zeropvec) # multinomial sampling
      pos=which(gibbsrlz==1)
      # }
      if(pos==1){
        I.start=I.new
      }
      FdatPL=pl2fullV2(I.start,datPL,nRe,phi.start,smlgamma.start[-c(1:m)])
    }
    ## update 1s I
    FdatPL=pl2fullV2(I.start,datPL,nRe,phi.start,smlgamma.start[-c(1:m)])
    for(j in (nRe-1):2){ # down or up for 1 step
      pos1=which(I.start==j)
      pos2=which(I.start==(j-1))
      pos3=which(I.start==(j+1))
      poses=c(pos2,pos3)

      nonzeropvec=rep(NA,3) # for nonzeros, it is a 2-dimensional multinomial distribution. rank j vs (j-1)
      I.new=replace(I.start,c(pos1,pos2),I.start[c(pos2,pos1)])

      nonzeropvec[1]=fulllikepower(cbind(dat,FdatPL),I = I.new,phi = phi.start,smlgamma = smlgamma.start)
      I.new=replace(I.start,c(pos1,pos3),I.start[c(pos3,pos1)])

      nonzeropvec[2]=fulllikepower(cbind(dat,FdatPL),I = I.new,phi = phi.start,smlgamma = smlgamma.start)
      nonzeropvec[3]=fulllikepower(cbind(dat,FdatPL),I = I.start,phi = phi.start,smlgamma = smlgamma.start)
      nonzeropvec=exp(nonzeropvec-(max(nonzeropvec)))
      nonzeropvec=nonzeropvec/sum(nonzeropvec)
      gibbsrlz=mc2d::rmultinomial(1, 1,nonzeropvec) # multinomial sampling
      pos=which(gibbsrlz==1)
      if(pos<3){ #
        I.start=replace(I.start,c(pos1,poses[pos]),I.start[c(poses[pos],pos1)])
      }
    }
    ## update I_r and phi
    ## the posterior likelihood can be written out,
    FdatPL=pl2fullV2(I.start,datPL,nRe,phi.start,smlgamma.start[-c(1:m)])
    tem=cbind(dat,FdatPL)
    Mallowsdat=tem[I.start>0,]
    Mallowsdat=apply(Mallowsdat,2,rank)

    ##################################
    phi.new=phi.start + phi.hyper* rnorm(1)
    if (phi.new>0 & phi.new<1){
      log.prob.start <- lapply(seq_len(ncol(Mallowsdat)), function(i) log(PerMallows::dmm(Mallowsdat[,i],I.start[I.start>0], -log(phi.start)*smlgamma.start[i])) )

      log.prob.new <- lapply(seq_len(ncol(Mallowsdat)), function(i) log(PerMallows::dmm(Mallowsdat[,i],I.start[I.start>0], -log(phi.new)*smlgamma.start[i])) )
      if (sum(unlist(log.prob.new))-sum(unlist(log.prob.start)) >log(runif(1))){
        phi.start=phi.new
      }
    }
    ## update smallgamma. MH is in use #############
    dattem=cbind(dat,FdatPL)
    for(j in 1:(m+mPL)){
      smlgamma.tem=smlgamma.start[j]
      smlgamma.tem=smlgamma.tem+gamma.hyper[j]*rnorm(1)
      if(smlgamma.tem>smlgamma.lower && smlgamma.tem< smlgamma.upper){
        like.start=BARDMallowslikepower(dattem[,j],I.start,phi.start,smlgamma.start[j])
        like.tem=BARDMallowslikepower(dattem[,j],I.start,phi.start,smlgamma.tem)
        if((like.tem-like.start) > log(runif(1))){
          smlgamma.start[j]=smlgamma.tem
        }
      }
    }
    FdatPL=pl2fullV2(I.start,datPL,nRe,phi.start,smlgamma.start[-c(1:m)])
    l.mat[i]=fulllikepower(dat = cbind(dat,FdatPL),I = I.start,phi = phi.start,smlgamma = smlgamma.start)

    I.mat[,i]=I.start
    phi.mat[i]=phi.start
    smlgamma.mat[,i]=smlgamma.start
    if(i==adaptation){
      gamma.hyper=sqrt(diag(cov(t(smlgamma.mat[,(adaptation-adaptation*0.6):adaptation]))))
      phi.hyper=sd((phi.mat[(adaptation-adaptation*0.6):adaptation]))
    }
  }
  return(list(I.mat=I.mat,phi.mat=phi.mat,smlgamma.mat=smlgamma.mat,l.mat=l.mat))
}
