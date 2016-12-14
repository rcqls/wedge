wedge.perf <- function(nbs.exp=c(10^6,10^7,10^8),tests=c("R","Rcpp","RcppParallel"),seed=10,N=3,nbs.threads=c(4,6,8),nbs.cores=c(4,6,8),summary.print=FALSE) {
  set.seed(seed)
  perfs<-list()
  for(nb.exp in nbs.exp) {
    a1 <- 10*runif(nb.exp)^2
    b1 <- 10*runif(nb.exp)^2
    a2 <- 10*runif(nb.exp)^2
    b2 <- 10*runif(nb.exp)^2

    # to be fair in the comparison do not take into an account the initialization
    size <- min(length(a1),length(a2),length(b1),length(b2))
    tau <- wedge:::wedge.threshold(N)

    if("Rcpp" %in% tests) {
      # wedge probabilities, 3 term approximation
      cat("Rcpp with nb.exp=",nb.exp,"\n")
      print(perfs[paste0("Rcpp-n=",nb.exp)] <- system.time(wpRcpp <- wedge(a1,b1,a2,b2,size,tau,N,type="Rcpp")))
      if(summary.print) print(summary(wpRcpp))
    }

    if("RcppParallel" %in% tests) {
      if("package:wedgeParallel" %in% search()) {
        for(nb.threads in nbs.threads) {
          RcppParallel::setThreadOptions(numThreads=nb.threads)
          cat("RcppParallel",nb.threads,"threads with nb.exp=",nb.exp,"\n")
          print(perfs[paste0("RcppParallel-n",nb.exp,"-t",nb.threads)] <- system.time(wpRcppParallel<- wedge(a1,b1,a2,b2,size,tau,N,type="RcppParallel",nb.threads=nb.threads)))
          if(summary.print) print(summary(wpRcppParallel))
        }
      } else {
        warning("Install wedgeParallel to parallelize the computation through RcppParallel.",immediate. = TRUE)
      }
    }

    if("R" %in% tests) {
      cat("Pure vectorized version with nb.exp=",nb.exp,"\n")
      if(nb.exp <= 10^7) {
        print(perfs[paste0("R-n",nb.exp)] <- system.time(wpR <- wedgeR(a1,b1,a2,b2,size,tau,N)))
        if(summary.print) print(summary(wpR))
      } else cat("Boom! We think that is a bit dangerous for your computer!\n")
    }
  }
}
