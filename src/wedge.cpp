#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


// #----------------------------------------------- compute threshold
// tau <- uniroot(bb,lower=0.5,upper=1.5)$root
// bb <- function(x){
// res <- exp(-8*x*(N-1)^2)/(4*x*(N-1))
// res <- res-(2/pi)^(3/2)*sqrt(x)/N*exp(2*x)*exp(-pi^2*N^2/(2*x))
// return(res)
// }

// [[Rcpp::export]]
double wedgeC0(double a1,double b1,double a2,double b2,double tau,int N=3,bool lower_tail=true) {
// #  Computes the wedge probability:
// #         P[-a1*t -b1 < Wt < a2*t +b2, for all t>0]
// #  where Wt is the standard Brownian motion.
// #  variables:
// #  a1, a2, b1, b2: vectors of numeric
// #  N: integer, number of terms to be computed.
// #  lower_tail: logical. If FALSE, the complement to 1
// #  (exit probability) is returned.
// #  Value: vector of numeric
// #
// #----------------------------------------------- check input sizes

  if (a1>0 && a2>0 && b1>0 && b2>0) {
    //----------------------------------------------- vectors of parameters

    double a1b1 = a1*b1, a2b2=a2*b2, a1b2 = a1*b2, a2b1 = a2*b1;
    double ap = (a1+a2)/2.0, am = (a1-a2)/2.0, bp = (b1+b2)/2.0, bm = (b1-b2)/2.0;
    double abp = ap*bp, abm = am*bm;
    double ga = (a1b1-a2b2)/2.0, de = (a1b2-a2b1)/2.0;


    //----------------------------------------------- initialize
    double w = 0.0;
    double n2,nm12,nnm1,nnp1;
    double An,Bn,Cn,Dn;

    double exp1,exp2;

    double pi=M_PI;//printf("pi=%lf\n",pi);
    //const double pi = std::atan(1.0)*4;

    if(abp>tau) {
      //------------------------------- sum 1
      for(int n=1;n<=N;n++) { //------------------------ main loop
        n2=n*n; nm12=n2-2*n+1; nnm1=n2-n; nnp1=n2+n;
        An = n2*a2b2+nm12*a1b1+nnm1*(a2b1+a1b2);
        Bn = nm12*a2b2+n2*a1b1+nnm1*(a2b1+a1b2);
        Cn = n2*(a1b1+a2b2)+nnm1*a2b1+nnp1*a1b2;
        Dn = n2*(a1b1+a2b2)+nnp1*a2b1+nnm1*a1b2;
        //printf("n2=%lf,nm12=%lf,nnm1=%lf,nnp1=%lf,An=%lf,Bn=%lf,Cn=%lf,Dn=%lf\n",n2,nm12,nnm1,nnp1,An,Bn,Cn,Dn);
        w += -exp(-2*An)-exp(-2*Bn)+exp(-2*Cn)+exp(-2*Dn);
        //printf("w1=%lf,w1=%lf\n",w,-exp(-2*An)-exp(-2*Bn)+exp(-2*Cn)+exp(-2*Dn));
      } //--------------------------------------------- end for loop
      return (lower_tail ? 1 + w : -w);
    } else {
      //------------------------------------ sum 2
      for(int n=1;n<=N;n++) { //------------------------ main loop
        exp1 = exp(-pow(2*pi*n,2.0)/(8.0*abp));
        exp1 = exp1*(cos(pi*n*de/abp)-cos(pi*n*ga/abp));
        exp2 = exp(-pow(pi*(2*n-1),2)/(8.0*abp));
        exp2 = exp2*(cos(pi*(n-0.5)*de/abp)+cos(pi*(n-0.5)*ga/abp));
        w += exp1+exp2;
        //printf("w2=%lf,w2=%lf,de=%lf,ga=%lf,abp=%lf,exp1=%lf,exp2=%lf\n",w,exp1+exp2,de,ga,abp,exp1,exp2);
      } //--------------------------------------------- end for loop
      w = w*sqrt(pi/(2*abp))*exp(pow(de,2)/(2*abp));
    }
    //----------------------------------------------- result
    return (lower_tail ? w : 1-w );
  } else {
    return (lower_tail ? 0.0 : 1.0 );
  }
}

// [[Rcpp::export]]
NumericVector wedgeC0vect(NumericVector a1,NumericVector b1,NumericVector a2,NumericVector b2,int N,double tau,bool lower_tail,int size_vect) {
  NumericVector w(size_vect);
  for(int i=0;i<size_vect;i++) w[i]=wedgeC0(a1[i],b1[i],a2[i],b2[i],tau,N,lower_tail);
  return w;
}
