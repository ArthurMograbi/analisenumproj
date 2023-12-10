#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "raiz.h"

#define VERBOSE 0

int bissecao (double a, double b, double (*f) (double x), double* r) {
  double fa = f(a);
  double fb = f(b);
  if(fa*fb>0){
    if(VERBOSE) printf("f(a)*f(b)>0!!\n");
    return 0;
  }
  
  *r = (a+b)/2;
  
  double emax = fabs((b-a)/2);
  if(emax<0.5e-8) {
    if(VERBOSE) printf("Erro máximo menor que 0.5e-8!\n");
    return 1;
  };
  
  double fres = f(*r);
  if(fabs(fres)<1e-12){
    if(VERBOSE) printf("|f(c)| menor que 1e-12!\nf(c)=%.10f\n",fres);
    return 1;
  };
  
  if(fa*fres<0)
     return bissecao(a,*r,f,r);    //raiz está entre a e ponto médio 
  return bissecao(*r,b,f,r);  //raiz está entre b e ponto médio 
}


void testwithvalues(double a,double b,double (*f)(double x)){
  double *r = (double *)malloc(sizeof(double));
  printf("Testing with:\na=%.2f\nb=%.2f\n",a,b);
  int rootfound = bissecao(a,b,f,r);
  if(rootfound){
    printf("Root found!\nAnswer = %.10f\nf(answer)=%.10f\n\n",*r,f(*r));
  } else {
    printf("Root not found!\n\n");
  }
  free(r);
}
