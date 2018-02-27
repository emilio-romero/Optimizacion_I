#include "funciones.h"

double Rosenbrock(double *x, int n){
double aux=0;
double sum1=0, sum2=0;
  for(int i=0;i<(n-1);i++){
    sum1+=(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]);
    sum2+=(1-x[i]*x[i])*(1-x[i]*x[i]);
  }
  aux=100*sum1 + sum2;
return(aux);}

double Wood(double *x, int n){
  double aux=0; 
  double t1,t2,t3,t4,t5; 
  t1=100*(x[0]*x[0]-x[1])*(x[0]*x[0]-x[1]); 
  t2=(x[0]-1)*(x[0]-1) + (x[2]-1)*(x[2]-1); 
  t3=90*(x[2]*x[2]-x[3])*(x[2]*x[2]-x[3]);
  t4=(x[1]-1)*(x[1]-1)+(x[3]-1)*(x[3]-1); t4=10.1*t4; 
  t5=19.8*((x[1]-1)*(x[3]-1));
  aux=t1+t2+t3+t4+t5;
return(aux);}

double SmoothingModel(double *x, int n, double lambda){
  double aux; 
  double suma1=0, suma2=0;
  for(int i=0;i<(n-1);i++){
    suma1+=(x[i]-x[i+n])*(x[i]-x[i+n]);
    suma2+=(x[i+1]-x[i])*(x[i+1]-x[i]);
  } 
  suma1+=(x[n-1]-x[n-1+n])*(x[n-1]-x[2*n-1]);
  aux=suma1+lambda*suma2;
return(aux);}
