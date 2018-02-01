#include <stdio.h>
#include <stdlib.h>
#include "algebralineal.h"
#include "lectura.h"
int main(){
int nr, nc; 
double N1; 
double **A=readMatrix("simetricaDefPos12x12.bin",&nr,&nc);
double **L=(double**)malloc(nr*sizeof(double*));
for(int i=0;i<nr;i++) L[i]=(double*)calloc(nc,sizeof(double));
Chol(A,nr,L);
printf("%dx%d\n",nr,nc);
es_spd("simetricaDefPos12x12.bin");
Norma_1_matriz(A,nr,nc,&N1);
printf("La norma 1 de la matriz es: %f\n",N1);
printf("Su programa ha terminado\n");
return 0;}
