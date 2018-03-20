#include <stdio.h>
#include "optimizadores06.h"
int main(int argc, char *argv[]){
Condiciones myc; 
myc.maxiter=10000; 
myc.tolg=1e-8; 
datos x0; 
x0.n=100; 
x0.x=(double*)malloc(x0.n*sizeof(double));

for(int i=0;i<x0.n;i++) x0.x[i]=1.0; 
x0.x[0]=-1.2; x0.x[x0.n-2]=-1.2; 
//x0.x[0]=-1.2; x0.x[1]=1.0; 

double *opti=MRC(Rosenbrock,gRosenbrock, hRosenbrock,x0,myc,1.0,0.085,0.22);
for(int i=0;i<10;i++) printf("%lf ",opti[i]);
printf("\n");
  printf("Su programa ha terminado\n");
return 0;}
