PAMA.B=function(datfile,nRe,iter=1000,init="EMM"){
  #' This function implements Bayesian inference of PAMA model.
  #'
  #' @export
  #' @import PerMallows
  #' @import stats
  #' @import mc2d
  #' @import ExtMallows
  #' @import rankdist

  #' @param datfile A matrix or dataframe. This is the data where our algorithm will work on. Each colomn denotes a ranker's ranking. The data should be in entity-based format.
  #' @param nRe A number. Number of relevant entities
  #' @param iter A number. Numner of iterations of MCMC
  #' @param init A string. This indicates which method is used to initiate the starting point of the aggregated ranking list. "mean" uses the sample mean. "EMM" uses the method from R package 'ExtMallows'.
  #' @return List. It contains Bayesian posterior samples of all the parameters and log-likelihood.
  #' \enumerate{
  #'   \item I.mat: posterior samples of I
  #'   \item phi.mat: posterior samples of phi
  #'   \item smlgamma.mat: posterior samples of gamma
  #'   \item l.mat: posterior samples of log-likelihood
  #' }
  #' @examples
  #' dat=t(PerMallows::rmm(10,1:20,0.5))
  #' PAMA.B(dat,10,iter=10)
  #' \donttest{PAMA.B(dat,10,iter=1000)}
  # this function implements Bayesian inference of PAMA model ().
  # The input
# parameter 'datfile': this is the data where our algorithm will work on. Each row denotes a ranker's ranking. The data should be in entity-based format.
  # parameter 'nRe' : number of relevant entities
  # parameter 'iter' : numner of iterations of MCMC

  # The output
  #  Output: I.mat: posterior samples of I
  # Output: phi.mat: posterior samples of \phi
  # Output: smlgamma.mat: posterior samples of \gamma
  # Output: l.mat: posterior samples of log-likelihood
  #' @author Wanchuang Zhu, Yingkai Jiang, Jun S. Liu, Ke Deng
  #' @references Wanchuang Zhu, Yingkai Jiang, Jun S. Liu, Ke Deng (2021) Partition-Mallows Model and Its Inference for Rank Aggregation. Journal of the American Statistical Association


  #source('conditionalranking.R')
  #source('PAMAlike.R')
  #source('fulllikepower.R')
  adaptation=0.25*iter
  dat=datfile
  m=dim(dat)[2] # the number of rankers
  n=dim(dat)[1] # number of entities

  ## hyperparameters
  sigma.square=m^(-1/2)
  gamma.hyper=rep(0.1,m)
  smlgamma.upper=10
  smlgamma.lower=0
  ## starting point
  phi.start=0.5
  phi.hyper=0.05
  # starting point of I
  if(init=="EMM"){
    mallowsinfer=ExtMallows::EMM((dat))
    mallowsinfer=as.numeric(mallowsinfer$op.pi0)
    mallowsinfer[mallowsinfer>nRe]=0
    I.start=mallowsinfer
  } else if(init=="mean"){
    mallowsinfer=PerMallows::lmm(t(dat),dist.name="kendall",estimation="approx")
    mallowsinfer=mallowsinfer$mode
    mallowsinfer[mallowsinfer>nRe]=0
    I.start=mallowsinfer
  }


  # starting point of gamma
  smlgamma.start=rep(runif(1,smlgamma.lower,smlgamma.upper),m)

  ## create matrix to store all the MCMC results
  I.mat=matrix(NA,n,iter)
  I.r.list=list()
  smlgamma.mat=matrix(NA,m,iter)
  phi.mat=matrix(NA,iter,1)
  l.mat=matrix(NA,iter,1)

  for(i in 1:iter){
    zerolist=which(I.start==0) # position of zeros
    for(j in zerolist){
      nonzerolist=which(I.start==nRe) # positions of non-zeros
      zeropvec=c() # for zeros, it is a 2-dimensional multinomial distribution. Either unchanged or change to nRe
      I.new=replace(I.start,c(j,nonzerolist),I.start[c(nonzerolist,j)])
      zeropvec=c(zeropvec,fulllikepower(dat = dat,I = I.new,phi = phi.start,smlgamma = smlgamma.start))
      zeropvec=c(zeropvec,fulllikepower(dat = dat,I = I.start,phi = phi.start,smlgamma = smlgamma.start)) # keep the zero unchanged
      zeropvec=exp(zeropvec-(max(zeropvec)))# fulllike returns log-likelihood
      if(sum(is.nan(exp(zeropvec-(min(zeropvec)))))>1){
        zeropvec=c(0.5,0.5)
      }else{
        zeropvec=zeropvec/sum(zeropvec)
      }
      gibbsrlz=mc2d::rmultinomial(1,1,zeropvec) # multinomial sampling
      pos=which(gibbsrlz==1)
      if(pos==1){ # zero should change to nRe
        I.start=I.new
      }
    }
    ## update 1s I
    pos=c()
    # nonzerolist=which(I.start>1)
    for(j in nRe:2){
      # newnonzerolist=which(I.start>0)
      pos1=which(I.start==j)
      pos2=which(I.start==(j-1))
      I.new=replace(I.start,c(pos1,pos2),I.start[c(pos2,pos1)])
      nonzeropvec=PAMA.likediff2R(dat,I.new,I.start,phi.start,smlgamma.start,nRe)

      if(nonzeropvec>log(runif(1))){ # zero should change to nRe
        I.start=I.new
      }
    }

    ## update I_r and phi
    ## the posterior likelihood can be written out,
    Mallowsdat=dat[I.start>0,]
    Mallowsdat=apply(Mallowsdat,2,rank)

    ##################################
    phi.new=phi.start+(phi.hyper* rnorm(1))
    if (phi.new>0 ){
      log.prob.start <- sapply(seq_len(ncol(Mallowsdat)), function(i){-DistancePair(Mallowsdat[,i],I.start[I.start>0])*phi.start*smlgamma.start[i] - logZ.MM(phi.start*smlgamma.start[i],nRe)} )
      log.prob.new <- sapply(seq_len(ncol(Mallowsdat)), function(i){-DistancePair(Mallowsdat[,i],I.start[I.start>0])*phi.new*smlgamma.start[i] - logZ.MM(phi.new*smlgamma.start[i],nRe)}  )
      if ((sum(log.prob.new)) >sum(log.prob.start)+log(runif(1))){
        phi.start=phi.new
      }
    }

    ## update smallgamma. MH is in use #############
    for(j in 1:m){
      smlgamma.tem=smlgamma.start[j]
      smlgamma.tem=smlgamma.tem+gamma.hyper[j]*rnorm(1)
      if(smlgamma.tem>0 ){
        like.start=PAMAlike(dat[,j],I.start,phi.start,smlgamma.start[j])
        like.tem=PAMAlike(dat[,j],I.start,phi.start,smlgamma.tem)
        if((like.tem) > like.start+log(runif(1))){
          smlgamma.start[j]=smlgamma.tem
        }
      }
    }

    I.mat[,i]=I.start
    phi.mat[i]=phi.start
    smlgamma.mat[,i]=smlgamma.start
    l.mat[i]=fulllikepower(dat = dat,I = I.start,phi = phi.start,smlgamma = smlgamma.start)
    if(i==adaptation){
      gamma.hyper=sqrt(diag(cov(t(smlgamma.mat[,(adaptation-adaptation*0.6):adaptation]))))
      phi.hyper=sd((phi.mat[(adaptation-adaptation*0.6):adaptation]))
    }
  }
  return(list(I.mat=I.mat,phi.mat=phi.mat,smlgamma.mat=smlgamma.mat,l.mat=l.mat))
}
