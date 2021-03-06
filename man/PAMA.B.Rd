% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAMA.B.R
\name{PAMA.B}
\alias{PAMA.B}
\title{This function implements Bayesian inference of PAMA model.}
\usage{
PAMA.B(datfile, nRe, iter = 1000, init = "EMM")
}
\arguments{
\item{datfile}{A matrix or dataframe. This is the data where our algorithm will work on. Each colomn denotes a ranker's ranking. The data should be in entity-based format.}

\item{nRe}{A number. Number of relevant entities}

\item{iter}{A number. Numner of iterations of MCMC}

\item{init}{A string. This indicates which method is used to initiate the starting point of the aggregated ranking list. "mean" uses the sample mean. "EMM" uses the method from R package 'ExtMallows'.}
}
\value{
List. It contains Bayesian posterior samples of all the parameters and log-likelihood.
\enumerate{
  \item I.mat: posterior samples of I
  \item phi.mat: posterior samples of phi
  \item smlgamma.mat: posterior samples of gamma
  \item l.mat: posterior samples of log-likelihood
}
}
\description{
This function implements Bayesian inference of PAMA model.
}
\examples{
dat=t(PerMallows::rmm(10,1:20,0.5))
PAMA.B(dat,10,iter=10)
\donttest{PAMA.B(dat,10,iter=1000)}
}
\references{
Wanchuang Zhu, Yingkai Jiang, Jun S. Liu, Ke Deng (2021) Partition-Mallows Model and Its Inference for Rank Aggregation. Journal of the American Statistical Association
}
\author{
Wanchuang Zhu, Yingkai Jiang, Jun S. Liu, Ke Deng
}
