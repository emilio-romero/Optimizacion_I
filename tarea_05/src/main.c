#include <stdio.h> 
#include "optimizadores.h"
#include <time.h>
int main(int argc, char *argv[]){
char ejercicio[30]; 
int bandera;
if(argc>1) strcpy(ejercicio,argv[1]);
if(argc>2) bandera=atoi(argv[2]); else bandera=1; 
srand(time(NULL));
Condiciones myc; 
myc.tolg=1e-3; myc.maxiter=15000;
myc.tolf=1e-10; myc.tolx=1e-10;
datos mx; mx.obs=0; 
int nnr; 

if(strcmp("Ej1r2",ejercicio)==0){
  mx.n=2;
  mx.x=(double*)malloc(2*sizeof(double));
  mx.x[0]=-1.2; mx.x[1]=1; 
  double *sol=SteepestDescent2(Rosenbrock,gRosenbrock,mx,myc,&bandera);
  datos xo;
  xo.n=2;
  xo.x=(double*)malloc(2*sizeof(double));
  xo.x[0]=1.0; xo.x[1]=1;
  vector_copiar(sol,2,mx.x);
  printf("|f(x*)-f(xk)|=%g",fabs(Rosenbrock(xo)-Rosenbrock(mx)));
  
  free(xo.x);
  free(mx.x);
  free(sol);
}
if(strcmp("Ej1r100",ejercicio)==0){
  mx.n=100; 
  mx.x=(double*)malloc(100*sizeof(double));
  for(int i=0;i<50;i++){ mx.x[2*i]=-1.2; mx.x[2*i+1]=1.0;} 
  double *sol=SteepestDescent2(Rosenbrock,gRosenbrock,mx,myc,&bandera);
  datos xo;
  xo.n=100;
  xo.x=(double*)malloc(100*sizeof(double));
  for(int i=0;i<100;i++) xo.x[i]=1.0; 
  vector_copiar(sol,100,mx.x);
  printf("|f(x*)-f(xk)|=%g",fabs(Rosenbrock(xo)-Rosenbrock(mx)));
  
  free(xo.x);
  free(mx.x);
  free(sol);
}
if(strcmp("Ej1w",ejercicio)==0){
  mx.n=4;
  mx.x=(double*)malloc(4*sizeof(double));
  mx.x[0]=-3.0; mx.x[1]=-1.0; mx.x[2]=-3.0; mx.x[3]=-1.0; 
  double *sol=SteepestDescent2(Wood,gWood,mx,myc,&bandera);
  datos xo;
  xo.n=4;
  xo.x=(double*)malloc(4*sizeof(double));
  xo.x[0]=1.0; xo.x[1]=1; xo.x[2]=1.0; xo.x[3]=1.0;
  vector_copiar(sol,4,mx.x);
  printf("|f(x*)-f(xk)|=%g",fabs(Wood(xo)-Wood(mx)));
  
  free(xo.x);
  free(mx.x);
  free(sol);
}

//double *sol=SteepestDescent2(Wood,gWood,mx,myc);

if(strcmp("Ej2",ejercicio)==0){
double *desire=(double*)malloc(2*sizeof(double)); desire[0]=0.0; desire[1]=1.0; 
mx.n=785; 
double **trainy=(double**)leerCSV("trainY.csv",50000,1);
int *mind=(int*)filtradoCSV(trainy,50000,desire,&nnr);
mx.obs=nnr;
double **trainx=leerCSV("trainX.csv",50000,784);
double *filty=(double*)aplicarFiltradoY(trainy,mind,nnr);
double **filtx=aplicarFiltradoX(trainx,784,mind,nnr);
int nnf=500; 
mx.dx=(double**)malloc(nnf*sizeof(double*));
for(int i=0;i<nnf;i++) mx.dx[i]=(double*)malloc(mx.n*sizeof(double));
mx.y=(double*)malloc(mx.n*sizeof(double));
reduccionMatriz(filtx,filty,mx.n,mx.obs,nnf,mx.dx,mx.y);
mx.obs=nnf;
/*Liberacion de datos crudos*/
for(int i=0;i<50000;i++){
  free(trainy[i]);
  free(trainx[i]);
} free(trainy); free(trainx);
free(mind);
for(int i=0;i<nnr;i++){
  free(filtx[i]);
} free(filtx); free(filty);
mx.x=(double*)malloc(785*sizeof(double));
for(int i=0;i<785;i++) mx.x[i]=0.02*randx()-0.02;
double *sol=SteepestDescent2(lLogistic,glLogistic,mx,myc,&bandera);//Optimizacion del proceso
double **testy=(double**)leerCSV("testY.csv",10000,1);
mind=(int*)filtradoCSV(testy,10000,desire,&nnr);
double **testx=leerCSV("testX.csv",10000,784);
filty=(double*)aplicarFiltradoY(testy,mind,nnr);
filtx=aplicarFiltradoX(testx,784,mind,nnr);
/*Liberacion de datos crudos*/
for(int i=0;i<10000;i++){
  free(testy[i]);
  free(testx[i]);
} free(testy); free(testx);
free(mind);
printf("El error es: %g \n",errorLogistico(sol,filtx,filty,785,nnr));//Calculo del error
for(int i=0;i<nnr;i++){
  free(filtx[i]);
} free(filtx); free(filty);
}


printf("\n");

//for(int i=0;i<mx.n;i++) printf("%lf ",sol[i]);
//printf("\n");
//printf("%1.2lf ",train[0][i]);

printf("Su programa ha terminado\n");
return(0);}
