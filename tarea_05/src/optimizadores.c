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
  double alpha=0.0005;
  double fr,fk,fkm1,gr,xr;
  double xk; int miter; 
  double *km2=(double*)malloc(mid.n*sizeof(double));
  double *km3=(double*)malloc(mid.n*sizeof(double));
  double *dx=(double*)malloc(mid.n*sizeof(double));
  vector_copiar(mid.x,mid.n,aux.x);
  g(aux,gaux); 
  printf("\nk & $||x_{k+1}-x_k||$ & $||\\nabla f(x_k) ||$ & $f(x_k)$ \\\\\\hline\n");
  //FILE *f1=fopen("e1mgf.dat","w");
  for(int k=0;k<micond.maxiter;k++){
    if(strcmp(micond.msg,"StepFijo")==0) alpha=paso_fijo(0.001);
    if(strcmp(micond.msg,"StepHess")==0) {
      h(aux,haux);
      alpha=paso_hessiano(gaux,haux,mid.n);
      if(alpha<0) alpha=fabs(alpha);
    }
    vector_escalar(alpha,gaux,mid.n,ag); 
    vector_resta(aux.x,ag,mid.n,ast.x);
    fk=f(aux); fkm1=f(ast);
    xk=Norma_2_vector(aux.x,aux.n);
    vector_resta(ast.x,aux.x,aux.n,dx);
    xr=Norma_2_vector(dx,ast.n);
    if(strcmp(micond.msg,"StepAprox")==0){
      alpha=paso_aproximacion(fkm1,fk,alpha,gaux,mid.n);
      if(alpha<0) alpha=fabs(alpha);
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
    //fprintf(f1,"%lf %lf %lf\n",ast.x[0],ast.x[1],fkm1);
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
  //fclose(f1);
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


double *SteepestDescent2(double(*f)(datos),int(*g)(datos,double*),datos x, Condiciones mc){
datos aux,ast; 

aux.n=x.n; ast.n=x.n; 
aux.x=(double*)malloc(aux.n*sizeof(double));
ast.x=(double*)malloc(ast.n*sizeof(double));
double paso,fk;
double *gk=(double*)malloc(aux.n*sizeof(double));
double *gaux=(double*)malloc(aux.n*sizeof(double));
vector_copiar(x.x,x.n,aux.x);
double normg,p1=1.0,p2=0.5; 
  fk=f(aux); 
  g(aux,gk); 
  for(int k=0;k<mc.maxiter;k++){
    /*Comprobaciones de tolerancia*/
      normg=Norma_2_vector(gk,aux.n);
      if(normg<=mc.tolg) break; 
  printf("Iteracion <%d> ",k);
    /*Calculo del tamanio de paso*/
    //paso=backtracking(f,aux,fk,gk);
    /*
    * Otros pasos (?)*/
    paso=quadInterpolation(f,aux,fk,gk,1.0);
    //paso=cubicInterpolation(f,aux,fk,gk,1.0,p2);
    printf("Paso: %lf\n",paso);
    /**/
    vector_escalar(paso,gk,aux.n,gaux); //paso por gradiente
    vector_resta(aux.x,gaux,aux.n,ast.x); //calculo del nuevo xk (x_k+1)
    fk=f(ast); //actualizacion del valor de la funcion
    g(ast,gk); //actualizacion del gradiente
    
    vector_copiar(ast.x,ast.n,aux.x); //Actualizacion de xk
  }

free(ast.x);
return(aux.x);}

double backtracking(double(*f)(datos),datos xk,double fk, double *gk){
  double a=1.0;
  double rho=0.5, c1=0.5;
  double nf; 
  double gtg=punto(gk,gk,xk.n);
  double *ng=(double*)malloc(xk.n*sizeof(double));
  double phia0;  
  datos nx; 
  nx.n=xk.n; 
  nx.x=(double*)malloc(nx.n*sizeof(double));
  do{
  a=rho*a; 
  vector_escalar(a,gk,xk.n,ng);
  vector_resta(xk.x,ng,xk.n,nx.x);
  phia0=f(nx);
}while(phia0>=(fk-c1*a*gtg));

return(a);}

double quadInterpolation(double(*f)(datos),datos xk,double fk, double *gk,double a){
  double a1,c1=0.5,a0=a; 
  double gtg=punto(gk,gk,xk.n); 
  double *ng=(double*)malloc(xk.n*sizeof(double));
  double phia1;  
  datos nx; 
  nx.n=xk.n; 
  nx.x=(double*)malloc(nx.n*sizeof(double));
    vector_escalar(a,gk,xk.n,ng);
    vector_resta(xk.x,ng,xk.n,nx.x);
    a1=a0*a0*gtg/(2*(f(nx)+a0*gtg-fk));
  do{
    if(fabs(a1-a0)<1e-5) break;
    a0=a1;
    vector_escalar(a0,gk,xk.n,ng);
    vector_resta(xk.x,ng,xk.n,nx.x);
    a1=a0*a0*gtg/(2*(f(nx)+a0*gtg-fk));
    phia1=f(nx); 
  }while(phia1>=(fk-c1*a1*gtg));

return(a1);}


double cubicInterpolation(double(*f)(datos),datos xk,double fk, double *gk,double pa,double sa){
double a0=pa,a1=sa,a2,c1=0.5;
double gtg=punto(gk,gk,xk.n);
double *ng=(double*)malloc(xk.n*sizeof(double));
double coef[2],*vec1=(double*)malloc(2*sizeof(double)),\
  **mat1=(double**)malloc(2*sizeof(double*));
for(int i=0;i<2;i++) mat1[i]=(double*)malloc(2*sizeof(double));
double phia2,phia1,phia0,den;  
  datos nx; 
  nx.n=xk.n; 
  nx.x=(double*)malloc(nx.n*sizeof(double));
  vector_escalar(a0,gk,xk.n,ng);
  vector_resta(xk.x,ng,xk.n,nx.x);
  phia0=f(nx); 
  vector_escalar(a1,gk,xk.n,ng);
  vector_resta(xk.x,ng,xk.n,nx.x);
  phia1=f(nx);
  mat1[0][0]=a0*a0;mat1[0][1]=-a1*a1;   
  mat1[1][0]=-a0*a0*a0;mat1[1][1]=a1*a1*a1; 
  vec1[0]=phia1+gtg*a1-fk; vec1[1]=phia0+gtg*a0-fk; 
  den=1.0/(a1*a1*a0*a0*(a1-a0));
  matriz_vector_mul(mat1,vec1,2,2,coef);
  vector_escalar(den,coef,2,coef);
  a2=(-coef[1]+sqrt(coef[1]*coef[1]+3*coef[0]*gtg))/(3*coef[0]);
  do{
    if(fabs(a2-a1)<1e-5) break;
    a0=a1; a1=a2; 
    vector_escalar(a0,gk,xk.n,ng);
    vector_resta(xk.x,ng,xk.n,nx.x);
    phia0=f(nx);
    vector_escalar(a1,gk,xk.n,ng);
    vector_resta(xk.x,ng,xk.n,nx.x);
    phia1=f(nx);
    mat1[0][0]=a0*a0;mat1[0][1]=-a1*a1;   
    mat1[1][0]=-a0*a0*a0;mat1[1][1]=a1*a1*a1; 
    vec1[0]=phia1+gtg*a1-fk; vec1[1]=phia0+gtg*a0-fk; 
    den=1.0/(a1*a1*a0*a0*(a1-a0));
    matriz_vector_mul(mat1,vec1,2,2,coef);
    vector_escalar(den,coef,2,coef);
 
   vector_escalar(a2,gk,xk.n,ng);
   vector_resta(xk.x,ng,xk.n,nx.x);
    a2=(-coef[1]+sqrt(coef[1]*coef[1]+3*coef[0]*gtg))/(3*coef[0]);
   phia2=f(nx);
  }while(phia1>=(fk-c1*a1*gtg));
 free(vec1); free(mat1[0]); free(mat1[1]); free(mat1); 
return(a2);}
double backtrackinge(double(*f)(datos),datos xk,double fk, double *gk){
  double a=1.0,a1,a2;
  double rho=randx(), c1=randx();
  double nf; 
  double gtg=punto(gk,gk,xk.n);
  double *ng=(double*)malloc(xk.n*sizeof(double));
  double coef[2],*vec1=(double*)malloc(2*sizeof(double)),\
  **mat1=(double**)malloc(2*sizeof(double*));
  for(int i=0;i<2;i++) mat1[i]=(double*)malloc(2*sizeof(double));
  double phia0, phia1,phia2,den;  
  datos nx; 
  nx.n=xk.n; 
  nx.x=(double*)malloc(nx.n*sizeof(double));
  do{
  a=rho*a; 
  vector_escalar(a,gk,xk.n,ng);
  vector_resta(xk.x,ng,xk.n,nx.x);
  phia0=f(nx);
  if(phia0<=(fk-c1*a*gtg) && c1>1e-4){
    free(ng);
    free(nx.x);
    return(a);
  }else{
    a1=a*a*gtg/(2*(f(nx)+a*gtg-fk));
    vector_escalar(a1,gk,xk.n,ng);
    vector_resta(xk.x,ng,xk.n,nx.x);
    phia1=f(nx);
    if(phia1<=(fk-c1*a1*gtg)){
      free(vec1);
      free(mat1[0]); free(mat1[1]); free(mat1);
      free(ng);
      free(nx.x);
      return(a1);
      } else{
        mat1[0][0]=a*a;mat1[0][1]=-a1*a1;   
        mat1[1][0]=-a*a*a;mat1[1][1]=a1*a1*a1; 
        vec1[0]=phia1+gtg*a1-fk; vec1[1]=phia0+gtg*a-fk; 
        den=1.0/(a1*a1*a*a*(a1-a));
        matriz_vector_mul(mat1,vec1,2,2,coef);
        vector_escalar(den,coef,2,coef);
        a2=(-coef[1]+sqrt(coef[1]*coef[1]+3*coef[0]*gtg))/(3*coef[0]);
        vector_escalar(a1,gk,xk.n,ng);
        vector_resta(xk.x,ng,xk.n,nx.x);
        phia2=f(nx);
        /*if(phia2<=(fk-c1*a1*gtg)){
          free(vec1);
          free(mat1[0]); free(mat1[1]); free(mat1);
          free(ng);
          free(nx.x);*/
          return(a2);
         /*}else{
          a=a2; 
         }*/
      }
  }

  }while(f(nx)<=(fk-c1*a*gtg));
free(vec1);
free(mat1[0]); free(mat1[1]); free(mat1);
free(ng);
free(nx.x);
return(a);}


