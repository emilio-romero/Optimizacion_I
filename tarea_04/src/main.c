#include <stdio.h>
#include "lectura.h"
#include "optimizadores.h"
int main(int argc, char *argv[]){

char inicio[30]; 
Condiciones tol; 
datos inicial; 
readParams(argc,argv,inicio,&tol.maxiter,&tol.tolg,&tol.tolx,&tol.tolf,tol.msg);
printf("%s\n",tol.msg);
int mn; 
inicial.x=leerVector(inicio,&inicial.n);
//inicial.n=mn;
//inicial.x=(double*)malloc(mn*sizeof(double));
//vector_copiar(mx,mn,inicial.x);
printf("%d\n%lf , %lf\n",inicial.n,inicial.x[0],inicial.x[1]);

//for(int i=0;i<mn;i++) mx[i]=1;
//mx[0]=-1.2; mx[1]=1.0;
//mx[0]=-3.0; mx[1]=-1.0;
//mx[2]=-3.0; mx[3]=-1.0;
double *xast; 
//escribirVector(mx,mn,"e3g.dat");
xast=SteepestDescent(Rosenbrock,gRosenbrock,inicial,tol);
for(int i=0;i<inicial.n;i++)
printf("%lf, ",xast[i]);
printf("\n");
/*Liberacion de memoria*/
free(inicial.x);

printf("Su programa ha terminado\n");
return 0;}
