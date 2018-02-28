#include <stdio.h>
#include "lectura.h"
#include "optimizadores.h"
int main(int argc, char *argv[]){

char inicio[30]; 
Condiciones tol; 
readParams(argc,argv,inicio,&tol.maxiter,&tol.tolg,&tol.tolx,&tol.tolf);

double *mx=(double*)malloc(2*sizeof(double));
mx[0]=3; mx[1]=-2.0;
double *xast; 

xast=SteepestDescent(Rosenbrock,gRosenbrock,mx,2,tol);

printf("%lf,%lf \n",xast[0],xast[1]);

/*Liberacion de memoria*/
free(mx);

printf("Su programa ha terminado\n");
return 0;}
