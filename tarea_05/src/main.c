#include <stdio.h> 
#include "optimizadores.h"
#include <time.h>
int main(int argc, char *argv[]){
srand(time(NULL));
Condiciones myc; 
myc.tolg=1e-5; myc.maxiter=200;
datos mx; 
mx.n=4; 
mx.x=(double*)malloc(mx.n*sizeof(double));
mx.x[0]=-3; mx.x[1]=-1; mx.x[2]=-3; mx.x[3]=-1;
double *sol=SteepestDescent2(Wood,gWood,mx,myc);
printf("%lf %lf %lf %lf \n", sol[0],sol[1],sol[2],sol[3]);

printf("Su programa ha terminado\n");
return(0);}
