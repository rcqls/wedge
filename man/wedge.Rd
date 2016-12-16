\name{wedge}
\alias{wedge}
\encoding{UTF-8}

\title{Computes wedge probabilities
}
\description{
Computes the probability for a standard Brownian motion to stay in a wedge
of the plane limited by two straight lines, given by their slopes and intercepts.
}
\usage{
wedge(a1, b1, a2, b2, N = 3, lower.tail = TRUE, type = c("R","Rcpp","RcppParallel"), nb.threads)
}
\arguments{
  \item{a1}{vector of numeric: slopes of lower boundaries
}
  \item{b1}{vector of numeric: intercepts of lower boundaries
}
  \item{a2}{vector of numeric: slopes of upper boundaries
}
  \item{b2}{vector of numeric: intercepts of upper boundaries
}
  \item{N}{integer: number of terms to be computed in each series
}
  \item{lower.tail}{logical: if \code{TRUE} the probabilities to stay inside the wedges
are returned; if \code{FALSE} the exit probabilities are returned.
}
  \item{type}{character: type of function to wrap
}
  \item{nb.threads}{integer: number of threads
}
}
\details{
Wedge probabilities are invariant through symmetry:
\code{wedge(a1,b1,a2,b2)} returns the same as
\code{wedge(b1,a1,b2,a2)} and
\code{wedge(a2,b2,a1,b1)}.
Depending on \code{N} a threshold is computed. If
(\code{a1}+\code{a2})*(\code{b1}+\code{b2})/4 is above the threshold, \code{N}
terms of Doob's formula are computed. If it is below the threshold, \code{N} terms
of the theta-function transform are computed.
}
\value{a vector of wedge probabilities.
}
\references{
B. Ycart and R. Drouilhet (2016) Computing wedge probabilities. Arxiv #16??.????
}
\author{
R\'emy Drouilhet and Bernard Ycart
}
\examples{
# random wedges
a1 <- 10*runif(10^6)^2
b1 <- 10*runif(10^6)^2
a2 <- 10*runif(10^6)^2
b2 <- 10*runif(10^6)^2
# wedge probabilities, 3 term approximation
pw3 <- wedge(a1,b1,a2,b2)
summary(pw3)
# wedge probabilities, 4 term approximation
pw4 <- wedge(a1,b1,a2,b2,N=4)
summary(abs(pw3-pw4))
#  exit probabilities, 3 term approximation
length(which(pw3==1))
pw3m <- wedge(a1,b1,a2,b2,lower.tail=FALSE)
summary(pw3+pw3m)
}
\keyword{distribution}
