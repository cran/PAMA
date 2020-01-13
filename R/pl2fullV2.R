pl2fullV2=function(I.start,PLdat,nRe,phi.start,smlgamma){
  # I.start is current true I
  # PLdat is partial list data  n*m1
  # nRe is number of relevant entities
  # phi.start is current phi of mallows model
  # smlgamma is the current gamma(quality parameter)
  # We assume people can only give the first several rankings, all the other are uniform distribution
  PLdat=as.matrix(PLdat)
  MHiter=1
  dat=PLdat
  n=dim(PLdat)[1]
  trank=I.start
  for(i in 1:dim(PLdat)[2]){
    ##########################   generate a compatible list

    P1=PLdat[,i]
    k=length(na.omit(P1))
    reP1=P1
    leftset=setdiff(1:n,as.vector(na.omit(reP1)))
    reP1[is.na(reP1)]=sample(leftset)

    #cat(is.permutation(reP1),'\t')
    ##########################   update the compatible list using MH
    knownP=P1[!is.na(P1)]
    for(j in 1:MHiter){
      for(jj in sample((k+2):n,1)){
        newlist=PerMallows::swap(reP1,which(reP1==jj),which(reP1==jj-1))
        like.start=fulllikepower(as.matrix(reP1),I = I.start,phi = phi.start,smlgamma = smlgamma[i])
        like.new=fulllikepower(as.matrix(newlist),I = I.start,phi = phi.start,smlgamma = smlgamma[i])
        if(like.new-like.start > log(runif(1)) ){
          reP1=newlist
        }
      }
    }

    dat[,i]=reP1
  }
  return(dat)
}

