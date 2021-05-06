logZ.MM=function(phi,n){sum(log(1-exp(-c(2:n)*phi)))-((n-1)*log(1-exp(-phi)))} # log normalizing constant in the Mallows model.
# phi is the disperse parameter
# n is the number of entities
