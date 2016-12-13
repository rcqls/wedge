wedge.perf <- function(nbs.exp=c(10^6,10^7,10^8),tests=c("R","Rcpp","RcppParallel"),seed=10,N=3,nbs.threads=c(4,6,8),nbs.cores=c(4,6,8),summary.print=FALSE) {
  set.seed(seed)
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
      print(system.time(wpRcpp <- wedge(a1,b1,a2,b2,size,tau,N,type="Rcpp")))
      if(summary.print) print(summary(wpRcpp))
    }

    if("RcppCores" %in% tests) {# This is in standby since we have to compute more carefully the number of block "b"
        require(parallel)
        for(nb.cores in nbs.cores) {
          cat("Rcpp",nb.cores,"cores with nb.exp=",nb.exp,"\n")
          print(system.time({
            size <- rep(nb.exp %/% nb.cores,nb.cores) + c(rep(1,r <- nb.exp %% nb.cores),rep(0,nb.cores - r))
            chunksSup <- cumsum(size)
            chunksInf <- 1 + c(0,chunksSup[-length(chunksSup)])
            wpRcppCore<-unlist(mclapply(1:nb.cores,function(i) wedge(a1[(chunksInf[i]:chunksSup[i])->ind],b1[ind],a2[ind],b2[ind],size[i],tau,type="Rcpp"),mc.cores=nb.cores))
          }))
          if(summary.print) print(summary(wpRcppCore))
        }
    }

    if("RcppParallel" %in% tests) {
      if("package:RcppParallel" %in% search()) {
        for(nb.threads in nbs.threads) {
          RcppParallel::setThreadOptions(numThreads=nb.threads)
          cat("RcppParallel",nb.threads,"threads with nb.exp=",nb.exp,"\n")
          print(system.time(wpRcppParallel<- wedge(a1,b1,a2,b2,size,tau,N,type="RcppParallel",nb.threads=nb.threads)))
          if(summary.print) print(summary(wpRcppParallel))
        }
      } else {
        warning("Install wedgeParallel to parallelize the computation through RcppParallel.",immediate. = TRUE)
      }
    }

    if("R" %in% tests) {
      cat("Pure vectorized version with nb.exp=",nb.exp,"\n")
      if(nb.exp <= 10^7) {
        print(system.time(wpR <- wedgeR(a1,b1,a2,b2,size,tau,N)))
        if(summary.print) print(summary(wpR))
      } else cat("Boom! We think that is a bit dangerous for your computer!\n")
    }
  }
}
