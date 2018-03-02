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
double *xt=(double*)malloc(243*sizeof(double));
if(strcmp(inicio,"e1g.dat")==0 || strcmp(inicio,"e2g.dat")==0 || strcmp(inicio,"e1m.dat")==0 || strcmp(inicio,"e2m.dat")==0){
  mf=&Rosenbrock; mg=&gRosenbrock; mh=&hRosenbrock;
  inicial.x=leerVector(inicio,&inicial.n);
}
if(strcmp(inicio,"e3g.dat")==0 || strcmp(inicio,"e3m.dat")==0){
  mf=&Wood; mg=&gWood; mh=&hWood;
  inicial.x=leerVector(inicio,&inicial.n);
}
if(strcmp(inicio,"yk.txt")==0){
  mf=&SmoothingModel; mg=&gSModel; mh=&hSModel;
  inicial.y=leeryk(inicio,&inicial.n,1);
  inicial.x=leeryk(inicio,&inicial.n,2);
  if(argc>7){
    inicial.param1=atof(argv[7]);
    printf("%lf\n",inicial.param1);
  } else inicial.param1=100;
  vector_copiar(inicial.x,inicial.n,xt);
  for(int i=0;i<inicial.n;i++)
    inicial.x[i]=10*randx();
}

xast=SteepestDescent(mf,mg,mh,inicial,tol);
//printf("%d\n", inicial.n);
//for(int i=242;i>237;i--) printf("%lf , %lf\n",inicial.x[i],inicial.y[i]);
//for(int i=0;i<inicial.n;i++)
//printf("%lf, ",xast[i]);
printf("\n");
/*Liberacion de memoria*/
//Especial para el ejercicio 4 
if(inicial.y!=NULL){
  char archivo[30];
  if(argc>8) strcpy(archivo,argv[8]);
  escribirEjer4(xast,xt,inicial.y,inicial.n,archivo);
}
free(inicial.x);
free(xast);
free(xt);
printf("Su programa ha terminado\n");

//double *mx=(double*)malloc(100*sizeof(double));
//for(int i=0;i<100;i++) mx[i]=-3.0+6*randx();
//escribirVector(mx,100,"e2m.dat");
return 0;}
