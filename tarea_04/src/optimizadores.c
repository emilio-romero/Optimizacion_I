#include "optimizadores.h"

double *SteepestDescent(double(*f)(datos),int(*g)(double*,int,double*),datos mid, Condiciones micond){
  datos aux, ast; 
  aux.n=mid.n; ast.n=mid.n;
  ast.x=(double*)malloc(mid.n*sizeof(double));
  aux.x=(double*)malloc(mid.n*sizeof(double));
  double *gaux=(double*)malloc(mid.n*sizeof(double));
  double *ag=(double*)malloc(mid.n*sizeof(double));
  double **haux=(double**)malloc(mid.n*sizeof(double*));
  for(int i=0;i<mid.n;i++) haux[i]=(double*)calloc(mid.n,sizeof(double));
  double alpha=0.05;
  double fr,fk,fkm1,gr,xr;
  vector_copiar(mid.x,mid.n,aux.x);
  g(aux.x,mid.n,gaux); 
  for(int k=0;k<micond.maxiter;k++){
  //printf("alfa=%lf\n",alpha);
    if(strcmp(micond.msg,"StepFijo")==0) alpha=paso_fijo(0.0001);
    if(strcmp(micond.msg,"StepHess")==0) {
      hRosenbrock(aux.x,mid.n,haux);
      alpha=paso_hessiano(gaux,haux,mid.n);
    }
    vector_escalar(alpha,gaux,mid.n,ag); 
    vector_resta(aux.x,ag,mid.n,ast.x);
    fk=f(aux); fkm1=f(ast);
    if(strcmp(micond.msg,"StepAprox")==0){
      alpha=paso_aproximacion(fkm1,fk,alpha,gaux,mid.n);
    } 
    g(ast.x,mid.n,gaux); 
    //Comprobaciones de paro de x, gradiente y funcion
    fr=fabs(fkm1-fk);
    
    if(1<fk) fr=fr/fabs(fk);
    gr=Norma_2_vector(gaux,mid.n);
    if(fr<micond.tolf || gr<micond.tolg)
      break; 
    vector_copiar(ast.x,mid.n,aux.x);
  }

free(aux.x);
free(gaux);
for(int i=0;i<mid.n;i++) free(haux[i]);
free(haux);
return(ast.x);}

/*Tamanios de paso*/
double paso_fijo(double a){
  return(a);
}

double paso_hessiano(double *g, double **H, int n){
  double aux;
  double *hg=(double*)malloc(n*sizeof(double));
  aux=punto(g,g,n);
  matriz_vector_mul(H,g,n,n,hg);
  aux=aux/punto(g,hg,n);
return(aux);}

double paso_aproximacion(double f, double fk, double ak, double *g,int n){
  double aux,gtg,den;
  gtg=punto(g,g,n);
  aux=ak*ak*gtg;// Utilizo g_k, pero ya esta modificada por el alpha 
  den=2.0*(f-fk+ak*gtg);
  //printf("Numerador=%lf, denomidador=%lf\n",aux,den);
  aux=aux/den;

return(aux);}
