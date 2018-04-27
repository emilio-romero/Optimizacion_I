#include "benchmark.h"

/*
* Rosenbrock
* Funcion
* Gradiente
* Hessiano 
*/ 

double Rosenbrock(double *x, int n){
  double aux=0; 
  double sum1=0, sum2=0; 
  for(int i=0;i<(n-1);++i){
    sum1+=(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]);
    sum2+=(1.0-x[i])*(1.0-x[i]); 
  }
  aux=100.0*sum1+sum2; 
return(aux);}

int gRosenbrock(double *x, int n, double *out){
out[0]=-400.0*x[0]*(x[1]-x[0]*x[0])-2.0*(1-x[0]);
for(int i=1;i<(n-1);++i){
    out[i]=200.0*(x[i]-x[i-1]*x[i-1])-400.0*x[i]*(x[i+1]-x[i]*x[i])-2.0*(1-x[i]);
  }
out[n-1]=200.0*(x[n-1]-x[n-2]*x[n-2]); 
return(1);}

int hRosenbrock(double *x, int n, double **out){

return(1);}
