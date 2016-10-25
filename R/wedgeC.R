wedgeC <-
function(a1,b1,a2,b2,N=3,lower.tail=TRUE){
  #  Computes the wedge probability:
  #         P[-a1*t -b1 < Wt < a2*t +b2, for all t>0]
  #  where Wt is the standard Brownian motion.
  #  variables:
  #  a1, a2, b1, b2: vectors of numeric
  #  N: integer, number of terms to be computed.
  #  lower.tail: logical. If FALSE, the complement to 1
  #  (exit probability) is returned.
  #  Value: vector of numeric
  #
  #----------------------------------------------- compute threshold
  bb <- function(x){
    res <- exp(-8*x*(N-1)^2)/(4*x*(N-1))
    res <- res-(2/pi)^(3/2)*sqrt(x)/N*exp(2*x)*exp(-pi^2*N^2/(2*x))
    return(res)
  }
  tau <- uniroot(bb,lower=0.5,upper=1.5)$root
  #cat("tau",tau,"\n")
  #----------------------------------------------- check input sizes
  size <- min(length(a1),length(a2),length(b1),length(b2))

  wedgeC0vect(a1,b1,a2,b2,N,tau,lower.tail,size)

}

wedgeCPar <- function(a1,b1,a2,b2,N=3,lower.tail=TRUE){
  #  Computes the wedge probability:
  #         P[-a1*t -b1 < Wt < a2*t +b2, for all t>0]
  #  where Wt is the standard Brownian motion.
  #  variables:
  #  a1, a2, b1, b2: vectors of numeric
  #  N: integer, number of terms to be computed.
  #  lower.tail: logical. If FALSE, the complement to 1
  #  (exit probability) is returned.
  #  Value: vector of numeric
  #
  #----------------------------------------------- compute threshold
  bb <- function(x){
    res <- exp(-8*x*(N-1)^2)/(4*x*(N-1))
    res <- res-(2/pi)^(3/2)*sqrt(x)/N*exp(2*x)*exp(-pi^2*N^2/(2*x))
    return(res)
  }
  tau <- uniroot(bb,lower=0.5,upper=1.5)$root
  #cat("tau",tau,"\n")
  #----------------------------------------------- check input sizes
  size <- min(length(a1),length(a2),length(b1),length(b2))

  wedgeC0Par(a1,b1,a2,b2,N,tau,lower.tail,size)

}
