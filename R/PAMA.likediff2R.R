PAMA.likediff2R=function(dat,I.new,I.start,phi,smlgamma,nRe){
  # This function returns the difference between the log-likelihood when two relevant entities are swapped.

  # #' @param dat Data matrix. Each column is a ranking list
  # #' @param I.new A ranking list. This is the ranking list after swapping two entities.
  # #' @param I.start A ranking list. This is the ranking list before swapping two entities.
  # #' @param phi A positive value. It is the common disperse parameter in Mallows model.
  # #' @param smlgamma A vector whose length is same as the numnber of rankers. This vector denotes the quality parameters for all the rankers.
  # #' @param alpha A value between [0,1]. This denotes the decay of the disperse parameter in the background entities compared to the relevant entities.
  # #' @param nRe An integer. This is the number of relevant entities.
  # #' @return
  # #'

  dat.re = dat[I.new>0,]
  dat.re=apply(dat.re,2,rank)
  mallows.start = sapply(seq_len(ncol(dat.re)), function(i) {-DistancePair(dat.re[,i],I.start[I.start>0])*phi*smlgamma[i]})
  mallows.start = sum(mallows.start)

  mallows.new = sapply(seq_len(ncol(dat.re)), function(i) {-DistancePair(dat.re[,i],I.new[I.new>0])*phi*smlgamma[i]})
  mallows.new = sum(mallows.new)
  return(mallows.new - mallows.start)
}
