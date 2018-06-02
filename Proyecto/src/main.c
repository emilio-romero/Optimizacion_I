#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include "gradienteconjugado.h" 
#include "funciones.h"
int main(int argc, char *argv[]){
int mn=5;
gsl_vector *x=gsl_vector_calloc(mn);
gsl_vector *sol=gsl_vector_calloc(mn);
for(int i=0;i<mn;++i){
  if(i%2==0){
    gsl_vector_set(x,i,-0.04);
  }else{
    gsl_vector_set(x,i,0.04);
  }
}
  GC_NTPA(x,0.00001,0.5,rastrigin,grastrigin,100,sol);
  for(int i=0;i<mn;++i){
    printf("%g ",gsl_vector_get(sol,i));
  }
printf("\n");

gsl_vector_free(sol);
gsl_vector_free(x);
printf("Su optimizacion ha terminado, tenga un buen dia!\n");
return(0);}
