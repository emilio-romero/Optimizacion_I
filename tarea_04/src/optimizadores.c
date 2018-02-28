#include "optimizadores.h"

double *SteepestDescent(double(*f)(double*,int),int(*g)(double*,int,double*),double *x0, int n, Condiciones micond){
  double *xast=(double*)malloc(n*sizeof(double));
  double *xaux=(double*)malloc(n*sizeof(double));
  double *gaux=(double*)malloc(n*sizeof(double));
  double alpha=0.0001;
  double fr,fk,fkm1,gr;
  vector_copiar(x0,n,xaux);
  g(xaux,n,gaux); 
  for(int k=0;k<micond.maxiter;k++){
    vector_escalar(alpha,gaux,n,gaux); 
    vector_resta(xaux,gaux,n,xast);
    g(xast,n,gaux); 
    
    //Comprobaciones de paro de x, gradiente y funcion
    fk=f(xaux,n); fkm1=f(xast,n);
    fr=fabs(fkm1-fk);
    if(1<fk) fr=fr/fabs(fk);
    gr=Norma_2_vector(gaux,n);
    if(fr<micond.tolf || gr<micond.tolg)
      break; 
    vector_copiar(xast,n,xaux);
  }

free(xaux);
free(gaux);
return(xast);}

/*Tamanios de paso*/
double paso_fijo(double a){
  return(a);
}

