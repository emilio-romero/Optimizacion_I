#include <stdio.h>
#include "lectura.h"
#include "optimizadores.h"
int main(int argc, char *argv[]){

char inicio[30]; 
Condiciones tol; 
datos inicial; 
double *xast; 
readParams(argc,argv,inicio,&tol.maxiter,&tol.tolg,&tol.tolx,&tol.tolf,tol.msg);
double (*mf)(datos);
int (*mg)(datos,double*);
int (*mh)(datos,double**);

if(strcmp(inicio,"e1g.dat")==0 || strcmp(inicio,"e2g.dat")==0 || strcmp(inicio,"e1m.dat")==0){
  mf=&Rosenbrock; mg=&gRosenbrock; mh=&hRosenbrock;
  inicial.x=leerVector(inicio,&inicial.n);
}
if(strcmp(inicio,"e3g.dat")==0 || strcmp(inicio,"e3m.dat")==0){
  mf=&Wood; mg=&gWood; mh=&hWood;
  inicial.x=leerVector(inicio,&inicial.n);
}
if(strcmp(inicio,"yk.txt")==0){
  mf=&SmoothingModel; mg=&gSModel; mh=&hSModel;
  inicial.x=leeryk(inicio,&inicial.n,1);
  inicial.y=leeryk(inicio,&inicial.n,2);
  if(argc>7){
    inicial.param1=atof(argv[7]);
    printf("%lf\n",inicial.param1);
  } else inicial.param1=100;
}

xast=SteepestDescent(mf,mg,mh,inicial,tol);
//printf("%d\n", inicial.n);
//for(int i=242;i>237;i--) printf("%lf , %lf\n",inicial.x[i],inicial.y[i]);
for(int i=0;i<inicial.n;i++)
printf("%lf, ",xast[i]);
printf("\n");
/*Liberacion de memoria*/
free(inicial.x);
free(xast);
printf("Su programa ha terminado\n");


//escribirVector(mx,mn,"e3g.dat");
return 0;}
