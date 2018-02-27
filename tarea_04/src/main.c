#include <stdio.h>
#include "lectura.h"
#include "optimizadores.h"
int main(int argc, char *argv[]){

char inicio[30]; 
int miter; 
double tolg, tolx, tolf;
readParams(argc,argv,inicio,&miter,&tolg,&tolx,&tolf);
double *mx=(double*)malloc(2*sizeof(double));
mx[0]=0.5; mx[1]=0.7;

printf("El vector inicial se extraera de: %s.\n",inicio);
printf("Las iteraciones maximas son %d\nLa tol. del gradiente %lf, de la x %lf y de la \
funcion %lf\n",miter, tolg, tolx, tolf);

printf("Prueba rosenbrock f(x)=%lf\n",Rosenbrock(mx,2));

printf("Su programa ha terminado\n");
return 0;}
