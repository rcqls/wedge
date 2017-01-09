\name{wedge-package}
\alias{wedge-package}
\encoding{UTF-8}

\docType{package}
\title{
\packageTitle{wedge}
}
\description{
\packageDescription{wedge}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{wedge}
\packageIndices{wedge}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
\packageAuthor{wedge}

Maintainer: \packageMaintainer{wedge}
}
\references{
B. Ycart and R. Drouilhet (2016) Computing wedge probabilities. arXiv:1612.05764
}
\keyword{distribution}
\examples{
# random wedges
a1 <- 10*runif(10^6)^2
b1 <- 10*runif(10^6)^2
a2 <- 10*runif(10^6)^2
b2 <- 10*runif(10^6)^2
# wedge probabilities, 3 term approximation
pw3 <- wedge(a1,b1,a2,b2)
}
