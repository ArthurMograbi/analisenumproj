#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define VERBOSE 0
#define PRINTERVAL 1000

double simpson(double a, double b, double (*f)(double x)){
  double h = b-a;
  double ans = (h/6.0)*(f(a)+4.0*f((a+b)/2.0)+f(b));
  return ans;
}

int count = 0;

double adaptsimpson (double a, double b, double (*f)(double x), double tol){
  double c= (b+a)/2.0;
  double orig = simpson(a,b,f);
  double split =simpson(a,c,f)+simpson(c,b,f);
  double delta =  orig-split;
  if(VERBOSE && (count%PRINTERVAL==0)) printf("%d - A(%.8f) C(%.8f) B(%.8f)\t Delta = %.6f\n",count,a,c,b,delta);
  if(fabs(delta)<=tol*15.0){
    if(VERBOSE) printf("%.8f<%.10f*15 returning %.10f\n",delta,tol,split-delta/15.0);
    return split-delta/15.0;
    count = 0;
  }
  count++;
  return adaptsimpson(a,c,f,tol/2.0)+adaptsimpson(c,b,f,tol/2.0);
}


