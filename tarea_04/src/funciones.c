#include "funciones.h"
/*
 *
 *Rosenbrock, funcion, gradiente y Hessiano
 * 
 */
double Rosenbrock(datos md){
double aux=0;
double sum1=0, sum2=0;
  for(int i=0;i<(md.n-1);i++){
    sum1+=(md.x[i+1]-md.x[i]*md.x[i])*(md.x[i+1]-md.x[i]*md.x[i]);
    sum2+=(1-md.x[i]*md.x[i])*(1-md.x[i]*md.x[i]);
  }
  aux=100*sum1 + sum2;
return(aux);}

int gRosenbrock(double *x, int n, double *out){
  for(int i=1;i<(n-1);i++){
    out[i]=-400*x[i]*(x[i+1]-x[i]*x[i])-2*(1-x[i])+200*(x[i]-x[i-1]*x[i-1]);
  }
  out[0]=-400*(x[1]-x[0]*x[0])*x[0]-2*(1-x[0]);
  out[n-1]=200*(x[n-1]-x[n-2]*x[n-2]);
return(1);}


int hRosenbrock(double *x, int n, double **out){
  for(int i=1;i<(n-1);i++){
    out[i][i]=202+1200*x[i]*x[i]-400*x[i+1];
    out[i][i-1]=-400*x[i-1];
    out[i-1][i]=-400*x[i-1];
  }
  out[0][0]=1200*x[0]*x[0]-400*x[1]+2;
  out[n-1][n-1]=200; 
  out[n-2][n-1]=-400*x[n-2];
  out[n-1][n-2]=-400*x[n-2];
return(1);}



/*
 *
 *Wood, funcion, gradiente y Hessiano
 * 
 */

double Wood(datos md){
  double aux=0; 
  double t1,t2,t3,t4,t5; 
  t1=100*(md.x[0]*md.x[0]-md.x[1])*(md.x[0]*md.x[0]-md.x[1]); 
  t2=(md.x[0]-1)*(md.x[0]-1) + (md.x[2]-1)*(md.x[2]-1); 
  t3=90*(md.x[2]*md.x[2]-md.x[3])*(md.x[2]*md.x[2]-md.x[3]);
  t4=(md.x[1]-1)*(md.x[1]-1)+(md.x[3]-1)*(md.x[3]-1); t4=10.1*t4; 
  t5=19.8*((md.x[1]-1)*(md.x[3]-1));
  aux=t1+t2+t3+t4+t5;
return(aux);}

int gWood(double *x, int n, double *out){
out[0]=200*(x[0]*x[0]-x[1])*2*x[0]+2*(x[0]-1);
out[1]=-200*(x[0]*x[0]-x[1])+20.2*(x[1]-1)+19.8*(x[3]-1);
out[2]=2*(x[2]-1)+360*(x[2]*x[2]-x[3])*x[2];
out[3]=-180*(x[2]*x[2]-x[3])+20.2*(x[3]-1)+19.8*(x[1]-1);

return(1);}

int hWood(double *x, int n, double **out){
out[0][0]=1200*x[0]*x[0]+2-400*x[1];
out[0][1]=out[1][0]=-400*x[0];
out[1][1]=220.2;
out[2][2]=2+1080*x[2]*x[2]-360*x[3]; 
out[3][3]=200.2; 
out[3][2]=out[2][3]=-360*x[2];
out[3][1]=out[1][3]=-19.8;
return(1);}



/*
 * Smoothing model
 * gradiente y hessiano
 *
 */



double SmoothingModel(datos md){
  double aux; 
  double suma1=0, suma2=0;
  for(int i=0;i<(md.n-1);i++){
    suma1+=(md.x[i]-md.x[i+md.n])*(md.x[i]-md.x[i+md.n]);
    suma2+=(md.x[i+1]-md.x[i])*(md.x[i+1]-md.x[i]);
  } 
  suma1+=(md.x[md.n-1]-md.x[md.n-1+md.n])*(md.x[md.n-1]-md.x[2*md.n-1]);
  aux=suma1+md.param1*suma2;
return(aux);}

int gSModel(double *x, int n, double *out){
  
return(1);}

int hSModel(double *x, int n, double **out){
  
return(1);}
