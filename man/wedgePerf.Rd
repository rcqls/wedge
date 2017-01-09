\name{wedge.perf}
\alias{wedge.perf}
\encoding{UTF-8}

\title{Performance of wedge probability computation
}
\description{
Returns running times for vectors of wedge probabilities under different implementations.
}
\usage{
wedge.perf(nbs.exp=c(10^6,10^7,10^8),tests=c("R","Rcpp","RcppParallel"),seed=10,N=3,nbs.threads=c(4,6,8),summary.print=FALSE)
}
\arguments{
  \item{nbs.exp}{vector of integer: numbers of wedge probabilities to be computed
}
  \item{tests}{vector of character: different implementations to be compared
}
  \item{seed}{numeric: random seed to fix
}
  \item{N}{integer: number of probabilities to be computed
}
  \item{nbs.threads}{vector of integer: numbers of threads
}
  \item{summary.print}{logical: if \code{TRUE} the summary of outputs is printed to check the results.
}
}
\details{
Compares the performance of different implementations for wedge probability computation.
}
\value{a list of system times.
}
\references{
B. Ycart and R. Drouilhet (2016) Computing wedge probabilities. arXiv:1612.05764
}
\author{
R\'emy Drouilhet and Bernard Ycart
}
\examples{
perfs <- wedge.perf()
perfs
}
