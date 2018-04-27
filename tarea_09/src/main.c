#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include "optimizacion.h"
#include "benchmark.h"
int main(int argc, char *argv[]){
  int n,cinicial; 
  int maxi; 
  char metodo[30]; 
  if(argc>1) n=atoi(argv[1]);
  if(argc>2) maxi=atoi(argv[2]);
  if(argc>3) strcpy(metodo,argv[3]);
  double *optimo=crear_vector(n);
  double *x0=crear_vector(n);
  double **h0=crear_matriz(n,n);
  matriz_identidad(n,n,h0);
  matriz_escalar(4.2,h0,n,n,h0);
 if(n==2){
    x0[0]=-1.2; x0[1]=1; 
  }
  if(n==100){
    for(int i=0;i<n;i++) {
      if(i%2==0){
        x0[i]=-1.2;
      }
      else{
        x0[i]=1.0;
      }
    }
  }

    BFGS(Rosenbrock,gRosenbrock,h0,x0,n,maxi,1e-5,optimo);
    printf("%lf %lf\n",optimo[0],optimo[1]);
  free(x0);
  free(optimo); 
  liberar_matriz(h0,n);
  printf("Su programa ha terminado\n");
return 0;}
