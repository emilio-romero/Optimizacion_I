#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include "optimizacion.h"
int main(int argc, char *argv[]){
  int cinicial; 
  char metodo[30]; 
  if(argc>1) cinicial=atoi(argv[1]);
  if(argc>2) strcpy(metodo,argv[2]);
  double *x0=crear_vector(3);
  double **a0=crear_matriz(3,3);
  double *xf=crear_vector(3);
  double tau=sqrt(DBL_EPSILON);
  if(cinicial==1){
    x0[0]=x0[1]=x0[2]=0.0;  
  }
  else if(cinicial==2){
    x0[0]=x0[1]=1.1; 
    x0[2]=-1.1; 
  }
  else if(cinicial==3){
    x0[0]=x0[1]=-10.0; 
    x0[2]=10.0; 
  }
  else if(cinicial==4){
    x0[0]=x0[2]=3.0;
    x0[1]=-3.0;
  }
  else {
    printf("Try: ./bin/ejecutable n, where n is an integer 1,2,3 or 4\n");
    return(-1);
  }
  
  if(strcmp(metodo,"newton")==0){
    NLSNewton(funt08,jact08,x0,3,300,tau,xf);
  }
  else if(strcmp(metodo,"broyden")==0){
    jact08(x0,3,a0);
    NLSBroyden(funt08,a0,x0,3,300,tau,xf); 
  }
  else{
    printf("Try: ./bin/ejecutable n metodo\n");
    return(-1);
  }
  //for(int i=0;i<3;++i) printf("%lf ",xf[i]);
  //printf("\n");
  liberar_matriz(a0,3);
  free(x0);
  free(xf);
  //printf("Su programa ha terminado\n");
return 0;}
