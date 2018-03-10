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
double *desire=(double*)malloc(2*sizeof(double)); desire[0]=0.0; desire[1]=1.0; 
int nnr; 

//double *sol=SteepestDescent2(Wood,gWood,mx,myc);
//printf("%lf %lf %lf %lf \n", sol[0],sol[1],sol[2],sol[3]);
double **trainy=(double**)leerCSV("trainY.csv",50000,1);
int *mind=(int*)filtradoCSV(trainy,50000,desire,&nnr);
//printf("%d %d %d %d\n",mind[1],mind[2],mind[3],nnr);

//double **yfilt=(double**)aplicarFiltradoY(trainy,mind,nnr);

double **trainx=leerCSV("trainX.csv",50000,784);
double **xfilt=(double**)aplicarFiltradoX(trainx,784,mind,nnr);


printf("Hola\n");
printf("%d\n",nnr);
//for(int i=0;i<500;i++)
//printf("%1.2lf ",train[0][i]);
for(int i=0;i<50000;i++){
  free(trainy[i]); 
} free(trainy);
free(mind);

printf("Su programa ha terminado\n");
return(0);}
