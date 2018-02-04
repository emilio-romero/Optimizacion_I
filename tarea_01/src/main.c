#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "algebralineal.h"
#include "lectura.h"
int main(int argc, char *argv[]){
char carchivo[30]; 
if(argc>1) strcpy(carchivo,argv[1]);
int nr, nc; 
double **A=readMatrix(carchivo,&nr,&nc);
double **L=(double**)malloc(nr*sizeof(double*));
for(int i=0;i<nr;i++) L[i]=(double*)calloc(nc,sizeof(double));
//Chol(A,nr,L);
printf("%dx%d\n",nr,nc);

Cholesky(carchivo,L);

printf("Su programa ha terminado\n");
return 0;}
