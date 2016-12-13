#ifndef wedgeC_H_
#define wedgeC_H_

//--------------------------------------------------- Non vectorized version of wedge probability
double wedgeC(double a1,double b1,double a2,double b2,double tau,int N=3,bool lower_tail=true) {

  if (a1>0 && a2>0 && b1>0 && b2>0) {
    //----------------------------------------------- vectors of parameters

    double a1b1 = a1*b1, a2b2=a2*b2, a1b2 = a1*b2, a2b1 = a2*b1;
    double ap = (a1+a2)/2.0, bp = (b1+b2)/2.0;
    double abp = ap*bp;
    double ga = (a1b1-a2b2)/2.0, de = (a1b2-a2b1)/2.0;

    //----------------------------------------------- initialize
    double w = 0.0;
    double n2,nm12,nnm1,nnp1;
    double An,Bn,Cn,Dn;

    double exp1,exp2;

    double pi=M_PI;

    if(abp>tau) {
      //------------------------------- sum 1
      for(int n=1;n<=N;n++) { //------------------------ main loop
        n2=n*n; nm12=n2-2*n+1; nnm1=n2-n; nnp1=n2+n;
        An = n2*a2b2+nm12*a1b1+nnm1*(a2b1+a1b2);
        Bn = nm12*a2b2+n2*a1b1+nnm1*(a2b1+a1b2);
        Cn = n2*(a1b1+a2b2)+nnm1*a2b1+nnp1*a1b2;
        Dn = n2*(a1b1+a2b2)+nnp1*a2b1+nnm1*a1b2;
        w += -exp(-2*An)-exp(-2*Bn)+exp(-2*Cn)+exp(-2*Dn);
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
      } //--------------------------------------------- end for loop
      w = w*sqrt(pi/(2*abp))*exp(pow(de,2)/(2*abp));
    }
    //----------------------------------------------- result
    return (lower_tail ? w : 1-w );
  } else {
    return (lower_tail ? 0.0 : 1.0 );
  }
}

#endif
