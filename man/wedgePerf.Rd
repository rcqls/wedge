\name{wedge.perf}
\alias{wedge.perf}
\encoding{UTF-8}

\title{Computes wedge probabilities
}
\description{
Computes the probability for a standard Brownian motion to stay in a wedge
of the plane limited by two straight lines, given by their slopes and intercepts.
}
\usage{
wedge.perf(nbs.exp=c(10^6,10^7,10^8),tests=c("R","Rcpp","RcppParallel"),seed=10,N=3,nbs.threads=c(4,6,8),summary.print=FALSE)
}
\arguments{
  \item{nbs.exp}{vector of integer: numbers of wedge probabilities to compute
}
  \item{tests}{vector of character: different functions to apply in the comparison
}
  \item{seed}{numeric: random seed to fix
}
  \item{N}{integer: number of terms to be computed in each series
}
  \item{nbs.threads}{vector of integer: numbers of threads
}
  \item{summary.print}{logical: if \code{TRUE} the summary is printed to check the results.
}
}
\details{
This function permits to compare the performance of the different functions computing the wedge probabilities.
}
\value{a list of system times.
}
\references{
B. Ycart (2016) Computing wedge probabilities. Arxiv #16??.????
}
\author{
Bernard Ycart
}
\examples{
perfs <- wedge.perf()
perfs
}
