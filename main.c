#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "raiz.h"
#include "simpson.h"

//Arthur Saddy Mograbi - 1910430

#define VERBOSE 0

double kx[] = {0,0,0,0};
double ky[] = {0,0,0,0}; 
double s = 0;


double * bezier(double *kx,double *ky,double t,double *out){
  double _t = 1;
  for(int _n=0;_n<4;_n++){
    out[0]+=kx[_n]*_t;
    out[1]+=ky[_n]*_t;
    _t*=t;
  } 
}

void setBezierConst(double *_kx,double *_ky){
  for(int i=0;i<4;i++){
    kx[i]=_kx[i]; ky[i]=_ky[i];
  }
}

double * bezierConst(double t,double *out){
  double _t = 1;
  for(int _n=0;_n<4;_n++){
    out[0]+=kx[_n]*_t;
    out[1]+=ky[_n]*_t;
    _t*=t;
  } 
}

double bezierConstX(double t){ 
  double _t = 1;
  double x= 0;
  for(int _n=0;_n<4;_n++){
    x+=kx[_n]*_t;
    _t*=t;
  } 
  return x;
}

double bezierConstY(double t){ 
  double _t = 1;
  double y= 0;
  for(int _n=0;_n<4;_n++){
    y+=ky[_n]*_t;
    _t*=t;
  }
  return y; 
}

double bezierConstYL(double t){ 
  double _t = 1;
  double y= 0;
  for(int _n=1;_n<4;_n++){
    y+=ky[_n]*_t*_n;
    _t*=t;
  }
  return y; 
}

double bezierConstXL(double t){ 
  double _t = 1;
  double x= 0;
  for(int _n=1;_n<4;_n++){
    x+=kx[_n]*_t*_n;
    _t*=t;
  }
  return x; 
}

void printCoord(double * coords){
  printf("X:%.5f,\tY:%.5f\n",coords[0],coords[1]);
}

int compCoords(double * c1, double * c2,double tol){
  for(int i=0;i<2;i++)
    if(fabs(c1[i]-c2[i])>tol) return 0;
  return 1;
}

double lenFunc(double t){
  double x_ = bezierConstXL(t);
  double y_ = bezierConstYL(t);
  return sqrt(x_*x_+y_*y_);
}

double arcLen(double t1, double t2){
  double len = adaptsimpson(t1,t2,&lenFunc,1e-8);
  return len;
}


double bissectFunc(double t){
  return s - arcLen(0.0,t);
}


int sToT(double s_,double *r){
  s = s_;
  return bissecao(0,50.0,&bissectFunc,r);
}

void sToPoint(double s_,double *out){
  double *t = (double*) malloc(sizeof(double));
  if(sToT(s_,t)){
    if(VERBOSE) printf("Found %.5f!\n",*t);
    bezierConst(*t,out);
  } else {
    if(VERBOSE) printf("No root found!\n");
  }
  free(t);
}


int testrandoms(int n){
  srand(n);
  double * r = (double *) malloc(sizeof(double));
  double *coords = (double *) malloc(2*sizeof(double));
  double *coords2 = (double *) malloc(2*sizeof(double));
  int found;
  double l,end;
  for(int i=0;i<n;i++){
    end = ((double)rand())/((double)RAND_MAX);
    *r=0;coords[0]=0;coords[1]=0;coords2[0]=0;coords2[1]=0;
    l = arcLen(0.0,end); 
    found = sToT(l,r);
    sToPoint(l,coords);
    bezierConst(end,coords2);
    if(!compCoords(coords,coords2,1e-6)) return 0;
  }
  free(r); free(coords); free(coords2);
  return 1;
}




int main (void)
{
  double kx[]={0.5,0.3,3.9,-4.7};//{0,1,0,0};
  double ky[]={1.5,0.3,0.9,-2.7};//{0,1,2,3};
  double end = 1.0;

  double * r = (double *) malloc(sizeof(double));
  double *coords = (double *) malloc(2*sizeof(double));
  double *coords2 = (double *) malloc(2*sizeof(double));

  setBezierConst(kx,ky);
  
  double l = arcLen(0.0,end); 
  printf("Length:%.10f\n",l);
  
  int found = sToT(l,r);
  printf("Root %s!\n",found?"found":"not found");
  printf("S:%.5f\tT:%.5f\n",l,*r);
  
  sToPoint(l,coords);
  printCoord(coords);
  bezierConst(end,coords2);
  printCoord(coords2);

  printf("Results are %s\n",compCoords(coords,coords2,1e-6)?"equal":"different");

  free(r); free(coords); free(coords2);

  printf("Tests %s\n",testrandoms(100)?"passed":"failed");

  return 0; 
} 
