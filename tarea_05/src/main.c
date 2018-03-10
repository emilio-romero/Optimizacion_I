#include <stdio.h> 
#include "optimizadores.h"
#include <time.h>
int main(int argc, char *argv[]){
srand(time(NULL));
Condiciones myc; 
myc.tolg=1e-5; myc.maxiter=5000;
datos mx; 
mx.n=785; 
//mx.x=(double*)malloc(mx.n*sizeof(double));
//mx.x[0]=-3; mx.x[1]=-1; mx.x[2]=-3; mx.x[3]=-1;
double *desire=(double*)malloc(2*sizeof(double)); desire[0]=0.0; desire[1]=1.0; 
int nnr; 

//double *sol=SteepestDescent2(Wood,gWood,mx,myc);
//printf("%lf %lf %lf %lf \n", sol[0],sol[1],sol[2],sol[3]);
double **trainy=(double**)leerCSV("trainY.csv",50000,1);
int *mind=(int*)filtradoCSV(trainy,50000,desire,&nnr);
mx.obs=nnr;
double **trainx=leerCSV("trainX.csv",50000,784);
double *filty=(double*)aplicarFiltradoY(trainy,mind,nnr);
double **filtx=aplicarFiltradoX(trainx,784,mind,nnr);
int nnf=50; 
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
for(int i=0;i<785;i++) mx.x[i]=randx();
printf("%lf %lf %lf\n",mx.y[0],mx.y[13],mx.y[60]);

//printf("Funcion evaluada 1 vez: %lf\n",lLogistic(mx));
double *sol=SteepestDescent2(lLogistic,glLogistic,mx,myc);


printf("\n");



//printf("%1.2lf ",train[0][i]);

printf("Su programa ha terminado\n");
return(0);}
