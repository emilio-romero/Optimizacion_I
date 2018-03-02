#include "optimizadores.h"

double *SteepestDescent(double(*f)(datos),int(*g)(datos,double*),int(*h)(datos,double**),datos mid, Condiciones micond){
  datos aux, ast; 
  aux.n=mid.n; ast.n=mid.n;
  ast.x=(double*)malloc(mid.n*sizeof(double));
  aux.x=(double*)malloc(mid.n*sizeof(double));
  //printf("%d\n",mid.n);
  if(mid.y!=NULL){
    ast.y=(double*)malloc(mid.n*sizeof(double));
    aux.y=(double*)malloc(mid.n*sizeof(double));
    vector_copiar(mid.y,mid.n,ast.y);
    vector_copiar(mid.y,mid.n,aux.y);
    ast.param1=mid.param1;
    aux.param1=mid.param1;
  }
  double *gaux=(double*)malloc(mid.n*sizeof(double));
  double *ag=(double*)malloc(mid.n*sizeof(double));
  double **haux=(double**)malloc(mid.n*sizeof(double*));
  for(int i=0;i<mid.n;i++) haux[i]=(double*)calloc(mid.n,sizeof(double));
  double alpha=0.005;
  double fr,fk,fkm1,gr,xr;
  double xk; int miter; 
  double *km2=(double*)malloc(mid.n*sizeof(double));
  double *km3=(double*)malloc(mid.n*sizeof(double));
  double *dx=(double*)malloc(mid.n*sizeof(double));
  vector_copiar(mid.x,mid.n,aux.x);
  g(aux,gaux); 
  printf("\nk & $||x_{k+1}-x_k||$ & $||\\nabla f(x_k) ||$ & $f(x_k)$ \\\\\\hline\n");
  for(int k=0;k<micond.maxiter;k++){
    if(strcmp(micond.msg,"StepFijo")==0) alpha=paso_fijo(0.0005);
    if(strcmp(micond.msg,"StepHess")==0) {
      h(aux,haux);
      alpha=paso_hessiano(gaux,haux,mid.n);
    }
    vector_escalar(alpha,gaux,mid.n,ag); 
    vector_resta(aux.x,ag,mid.n,ast.x);
    fk=f(aux); fkm1=f(ast);
    xk=Norma_2_vector(aux.x,aux.n);
    vector_resta(ast.x,aux.x,aux.n,dx);
    xr=Norma_2_vector(dx,ast.n);
    if(strcmp(micond.msg,"StepAprox")==0){
      alpha=paso_aproximacion(fkm1,fk,alpha,gaux,mid.n);
    } 
    g(ast,gaux); 
    gr=Norma_2_vector(gaux,mid.n);
    //Impresiones para las tablas 
    if(k<3){
      printf("%d & %lf & %lf & %lf \\\\ \n",k,xr,gr,fk);
      //printVector(gaux,mid.n);
    }
    //Comprobaciones de paro de x, gradiente y funcion
    fr=fabs(fkm1-fk); 
    if(1<fabs(fk)) fr=fr/fabs(fk);
    if(1<fabs(xk)) xr=xr/fabs(xk);
    if(fr<micond.tolf || gr<micond.tolg || xr<micond.tolx){ 
      miter=k; 
      //vector_copiar(aux.x,mid.n,km2);
      //vector_copiar(ast.x,mid.n,aux.x);
      break;} 
      vector_copiar(km2,mid.n,km3);
      vector_copiar(aux.x,mid.n,km2);
      vector_copiar(ast.x,mid.n,aux.x);
  miter=k;
  }
  datos am2; am2.n=aux.n; 
  if(mid.y!=NULL){
    am2.y=(double*)malloc(mid.n*sizeof(double));
    //aux.y=(double*)malloc(mid.n*sizeof(double));
    vector_copiar(mid.y,mid.n,am2.y);
    am2.param1=mid.param1;
    //vector_copiar(mid.y,mid.n,aux.y);
  }
  am2.x=(double*)malloc(am2.n*sizeof(double));
  vector_copiar(km2,am2.n,am2.x);
  vector_resta(km2,km3,aux.n,dx); 
  g(am2,gaux);
  printf("%d & %g & %g & %g \\\\ \n",miter-2,Norma_2_vector(dx,mid.n),Norma_2_vector(gaux,mid.n),f(am2));
  vector_resta(aux.x,km2,aux.n,dx); 
  g(aux,gaux);
  printf("%d & %g & %g & %g \\\\ \n",miter-1,Norma_2_vector(dx,mid.n),Norma_2_vector(gaux,mid.n),f(aux));
  vector_resta(ast.x,aux.x,aux.n,dx); 
  g(ast,gaux);
  printf("%d & %g & %g & %g \\\\ \n",miter,Norma_2_vector(dx,mid.n),Norma_2_vector(gaux,mid.n),f(ast));
  printf("%d\n",miter); 

free(aux.x);
free(gaux);
for(int i=0;i<mid.n;i++) free(haux[i]);
free(haux);
free(km2); free(km3);
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
