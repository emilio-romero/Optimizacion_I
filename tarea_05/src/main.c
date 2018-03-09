#include <stdio.h> 
#include "optimizadores.h"
#include <time.h>
int main(int argc, char *argv[]){
srand(time(NULL));
Condiciones myc; 
myc.tolg=1e-5; myc.maxiter=5000;
datos mx; 
mx.n=4; 
mx.x=(double*)malloc(mx.n*sizeof(double));
mx.x[0]=-3; mx.x[1]=-1; mx.x[2]=-3; mx.x[3]=-1;
//double *sol=SteepestDescent2(Wood,gWood,mx,myc);
//printf("%lf %lf %lf %lf \n", sol[0],sol[1],sol[2],sol[3]);
//double **train=leerCSV("trainX.csv",50000,784);
//for(int i=0;i<500;i++)
//printf("%1.2lf ",train[0][i]);
printf("Su programa ha terminado\n");
return(0);}
