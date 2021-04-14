PAMA.F=function(datfile,nRe,threshold,iter=1000){
  #' This function implements Maximum Likelihood estimation of PAMA model.
  #'
  #' @export
  #' @import PerMallows
  #' @import stats
  #' @import mc2d
  #' @param datfile A matrix or dataframe. This is the data where our algorithm will work on. Each row denotes a ranker's ranking. The data should be in entity-based format.
  #' @param nRe A number. Number of relevant entities.
  #' @param iter A number. Numner of iterations of MCMC.
  #' @param threshold A number(positive). The stopping threshold in determining convergence of MLE. if the two consecutive iterations of log-likelihood is smaller than threshold, then the convergence achives.

  #' @return List. It contains MLE of all the parameters and log-likelihood.
  #' \enumerate{
  #' \item I.mat:  samples of I
  #' \item phi.mat:  samples of phi.
  #' \item smlgamma.mat:  samples of gamma
  #' \item l.mat:  samples of log-likelihood
  #' }
  #' @examples
  #' a=NBANFL()
  #' PAMA.F(a$NBA,nRe=10,threshold=0.1,iter=100)
  #' @author Wanchuang Zhu, Yingkai Jiang, Jun S. Liu, Ke Deng
  #' @references Wanchuang Zhu, Yingkai Jiang, Jun S. Liu, Ke Deng (2021) Partition-Mallows Model and Its Inference for Rank Aggregation *Journal of the American Statistical Association*


  dat=datfile
  #source('gprimegammak.R')
  #source('gprimeprimegammak.R')
  #source('gprimephi.R')
  #source('gprimeprimephi.R')
  #source("BARDMallowslikepower.R")
  #source('conditionalranking.R')

  n=dim(dat)[1]
  m=dim(dat)[2]
  smlgamma.upper=10
  smlgamma.lower=0.01
  ## starting point of I
  phi.start=0.3
  mallowsinfer=PerMallows::lmm(t(dat),dist.name="kendall",estimation="approx")
  mallowsinfer=mallowsinfer$mode
  mallowsinfer[mallowsinfer>nRe]=0
  I.start=mallowsinfer

  smlgamma.start=rep(runif(1,smlgamma.lower,smlgamma.upper),m)
  ## create matrix to store all the MCMC results
  I.mat=matrix(NA,n,iter)
  I.r.list=list()
  smlgamma.mat=matrix(NA,m,iter)
  phi.mat=matrix(NA,iter,1)
  l.mat=matrix(NA,iter,1)
  alpha.gamma=0.01
  for(i in 1:iter){
    ## update smallgamma.
    for(j in 1:m){
      likeold=fulllikepower(dat = dat,I = I.start,phi = phi.start,smlgamma = smlgamma.start)
      tem=smlgamma.start[j]-alpha.gamma*gprimegammak(dat[,j],I.start,phi.start,smlgamma.start[j])/gprimeprimegammak(dat[,j],I.start,phi.start,smlgamma.start[j])
      smlgamma.tem=smlgamma.start
      smlgamma.tem[j]=tem

      while(tem<smlgamma.lower || tem > smlgamma.upper || (likeold>fulllikepower(dat = dat,I = I.start,phi = phi.start,smlgamma = smlgamma.tem))){
        alpha.gamma=alpha.gamma/2
        tem=smlgamma.start[j]-alpha.gamma*gprimegammak(dat[,j],I.start,phi.start,smlgamma.start[j])/gprimeprimegammak(dat[,j],I.start,phi.start,smlgamma.start[j])
        smlgamma.tem[j]=tem
      }
      smlgamma.start[j]=tem
      alpha.gamma=0.01
    }


    ################################## update phi
  alpha.phi=0.0001
  likeold=fulllikepower(dat = dat,I = I.start,phi = phi.start,smlgamma = smlgamma.start)
    tem=phi.start-alpha.phi* gprimephi(dat,I.start,phi.start,smlgamma.start)/gprimeprimephi(dat,I.start,phi.start,smlgamma.start)
    likenew=fulllikepower(dat = dat,I = I.start,phi = tem,smlgamma = smlgamma.start)
    if (is.na(likenew)){
      tem=phi.start
    }else{
      while(tem<0 || tem>1 || likeold>likenew){
        alpha.phi=alpha.phi/2
        tem=phi.start-alpha.phi* gprimephi(dat,I.start,phi.start,smlgamma.start)/gprimeprimephi(dat,I.start,phi.start,smlgamma.start)
        likenew=fulllikepower(dat = dat,I = I.start,phi = tem,smlgamma = smlgamma.start)
      }
    }

    phi.start=tem
    ############################  update I
  likeold=fulllikepower(dat = dat,I = I.start,phi = phi.start,smlgamma = smlgamma.start)
    zerolist=which(I.start==0) # position of zeros
    for(j in zerolist){
      nonzeros=which(I.start>0)

      nonzerolist=which(I.start==nRe)#sample(nonzeros)  # positions of non-zeros

        I.new=replace(I.start,c(j,nonzerolist),I.start[c(nonzerolist,j)])
        likenew=fulllikepower(dat = dat,I = I.new,phi = phi.start,smlgamma = smlgamma.start)
        if(likenew>likeold){ # zero should change to nRe
          I.start=I.new
          likeold=likenew
          break
        }
    }
    ## update 1s I
    pos=c()
    for(j in nRe:2){
      # newnonzerolist=which(I.start>0)
      pos1=which(I.start==j)
      pos2=which(I.start==(j-1))
      I.new=replace(I.start,c(pos1,pos2),I.start[c(pos2,pos1)])

      likenew=fulllikepower(dat = dat,I = I.new,phi = phi.start,smlgamma = smlgamma.start)
      if(likenew>likeold){ # zero should change to nRe
        I.start=I.new
        likeold=likenew
      }
    }

    I.mat[,i]=I.start
    phi.mat[i]=phi.start
    smlgamma.mat[,i]=smlgamma.start
    l.mat[i]=fulllikepower(dat = dat,I = I.start,phi = phi.start,smlgamma = smlgamma.start)

    if(i>1 &&l.mat[i]!=-Inf && l.mat[i-1]!=-Inf){
      if ((abs(l.mat[i]-l.mat[i-1])<threshold)){
        break
      }
    }
    if(l.mat[i]==-Inf && l.mat[i-1]==-Inf && i>iter){
      break
    }
  }

  return(list(I.mat=I.mat,phi.mat=phi.mat,smlgamma.mat=smlgamma.mat,l.mat=l.mat))
}
