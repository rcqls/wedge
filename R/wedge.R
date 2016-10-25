wedge <-
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
a1r <- a1[1:size]
b1r <- b1[1:size]
a2r <- a2[1:size]
b2r <- b2[1:size]
#----------------------------------------------- initialize result
if(lower.tail){w <- rep(0,size)
   }else{w <- rep(1,size)}
#----------------------------------------------- positive values
indnn <- (a1r>0)&(a2r>0)&(b1r>0)&(b2r>0)
#----------------------------------------------- vectors of parameters
a1r <- a1r[indnn]
b1r <- b1r[indnn]
a2r <- a2r[indnn]
b2r <- b2r[indnn]
a1b1 <- a1r*b1r
a2b2 <- a2r*b2r
a1b2 <- a1r*b2r
a2b1 <- a2r*b1r
ap <- (a1r+a2r)/2
am <- (a1r-a2r)/2
bp <- (b1r+b2r)/2
bm <- (b1r-b2r)/2
abp <- ap*bp
abm <- am*bm
ga <- (a1b1-a2b2)/2
de <- (a1b2-a2b1)/2
#----------------------------------------------- which sum to use
ind1 <- (abp>tau)                              # sum 1
ind2 <- !ind1                                  # sum 2
#----------------------------------------------- terms for each sum
a1b1 <- a1b1[ind1]                             # for sum 1
a2b2 <- a2b2[ind1]
a1b2 <- a1b2[ind1]
a2b1 <- a2b1[ind1]
abp2 <- abp[ind2]                              # for sum2
abm2 <- abm[ind2]
de2 <- de[ind2]
ga2 <- ga[ind2]
#----------------------------------------------- initialize
w1 <- rep(0,sum(ind1))
w2 <- rep(0,sum(ind2))
#----------------------------------------------- main loop
for(n in 1:N){
#----------------------------------------------- sum 1
   n2<-n^2; nm12<-(n-1)^2; nnm1<-n*(n-1); nnp1<-n*(n+1)
   An <- n2*a2b2+nm12*a1b1+nnm1*(a2b1+a1b2)
   Bn <- nm12*a2b2+n2*a1b1+nnm1*(a2b1+a1b2)
   Cn <- n2*(a1b1+a2b2)+nnm1*a2b1+nnp1*a1b2
   Dn <- n2*(a1b1+a2b2)+nnp1*a2b1+nnm1*a1b2
   #cat("n2=",n2,"nm12=",nm12,"nnm1=",nnm1,"nnp1=",nnp1,"An=",An,"Bn=",Bn,"Cn=",Cn,"Dn=",Dn,"\n")

   w1 <- w1-exp(-2*An)-exp(-2*Bn)+exp(-2*Cn)+exp(-2*Dn)
   #cat("w1=",w1,"w1=",-exp(-2*An)-exp(-2*Bn)+exp(-2*Cn)+exp(-2*Dn),"\n")
#----------------------------------------------- sum 2
   exp1 <- exp(-(2*pi*n)^2/(8*abp2))
   exp1 <- exp1*(cos(pi*n*de2/abp2)-cos(pi*n*ga2/abp2))
   exp2 <- exp(-(pi*(2*n-1))^2/(8*abp2))
   exp2 <- exp2*(cos(pi*(n-0.5)*de2/abp2)+cos(pi*(n-0.5)*ga2/abp2))
   w2 <- w2+exp1+exp2
  #cat("w2=",w2,"w2=",exp1+exp2,"de=",de2,"ga=",ga2,"abp=",abp2,"exp1=",exp1,"exp2=",exp2,"\n")
}                                              # end for
w2 <- w2*sqrt(pi/(2*abp2))*exp(de2^2/(2*abp2))
#----------------------------------------------- result
if(lower.tail){
   w1 <- 1+w1
   }else{
   w1 <- (-w1); w2 <- 1-w2
   }
wr <- w[indnn]
wr[ind1] <- w1
wr[ind2] <- w2
w[indnn] <- wr
return(w)
}
