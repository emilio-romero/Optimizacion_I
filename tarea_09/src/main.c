#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include "optimizacion.h"
#include "benchmark.h"
int main(int argc, char *argv[]){
  int n,cinicial; 
  int maxi;
  int code; 
  double h; 
  char metodo[30]; 
  if(argc>1) n=atoi(argv[1]);
  if(argc>2) maxi=atoi(argv[2]);
  if(argc>3) code=atoi(argv[3]); else code=1; 
  if(argc>4) h=atof(argv[4]); else h=0.01; 
  if(argc>3) strcpy(metodo,argv[3]);
  double *optimo=crear_vector(n);
  double *x0=crear_vector(n);
  double **h0=crear_matriz(n,n);
  double **h0i=crear_matriz(n,n);
 /*
  matriz_identidad(n,n,h0);
  matriz_escalar(2.2,h0,n,n,h0);
 /**/
  //if(n==100){
    for(int i=0;i<n;i++) {
      if(i%2==0){
        x0[i]=-1.2;
      }
      else{
        x0[i]=1.0;
      }
    }
 // }
 /**/
  if(code==1){
    hRosenbrock(x0,n,h0);
  }
  if(code==2){
    aproximaHessiano(Rosenbrock,x0,n,h,h0);
  }
  matriz_inversa(h0,n,h0i);  
  /**/
  
  BFGS(Rosenbrock,gRosenbrock,h0i,x0,n,maxi,1e-6,optimo);
  
  /*Liberacion*/
  free(x0);
  free(optimo); 
  liberar_matriz(h0,n);
  liberar_matriz(h0i,n);
  printf("Su programa ha terminado\n");
return 0;}
