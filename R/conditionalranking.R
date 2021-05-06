conditionalranking = function(I,bsrkr){
  # bsrkr is the observed base ranker
  # I is the true classification of entities
  rank_RE=bsrkr[I!=0]
  rank_BE=bsrkr[I==0]
  tau01=sapply(seq_len(length(rank_BE)), function(i) (sum(rank_BE[i]>rank_RE)))
  tau01=length(rank_RE)-tau01
  return(tau01)
}


