#include "optimizacion.h"

int NLSNewton(int(*F)(double *x, int n, double *out),int(*J)(double *x, int n, double **out),
    double *x0, int n, int maxiter, double tol,double *out){
  double *xk=crear_vector(n);
  double *sk=crear_vector(n);
  double *fk=crear_vector(n);
  double **jk=crear_matriz(n,n);
  double nf,kappa; 
  for(int k=0;k<maxiter;++k){
    F(xk,n,fk); 
    nf=Norma_2_vector(fk,n);
    vector_escalar(-1.0,fk,n,fk); 
    if(nf<tol) break; 
    J(xk,n,jk);
    solLU(jk,fk,n,n,sk);
    vector_suma(xk,sk,n,xk);
    kappa=numero_condicion(jk,n);
    printf("%d ",k);
    for(int i=0;i<n;i++) printf("%lf ",xk[i]);
    printf("%g %g\n",nf,kappa);
  }
  vector_copiar(xk,n,out);
  free(sk);
  free(xk); 
  free(fk); 
  liberar_matriz(jk,n);
return(1);}
int NLSBroyden(int(*F)(double *x, int n, double *out),double **A0,
    double *x0, int n, int maxiter, double tol,double *out){
  double *xk=crear_vector(n);
  double *sk=crear_vector(n);
  double *yk=crear_vector(n);
  double *fk=crear_vector(n);
  double *auxv=crear_vector(n);
  double **Ak=crear_matriz(n,n);
  double **auxm=crear_matriz(n,n);
  double nf,kappa; 
  matriz_copiar(A0,n,n,Ak);
  F(xk,n,fk);  
  for(int k=0;k<maxiter;++k){
    kappa=numero_condicion(Ak,n);
    nf=Norma_2_vector(fk,n);
    if(nf<tol) break; 
    vector_escalar(-1.0,fk,n,fk); 
    solLU(Ak,fk,n,n,sk);
    vector_suma(xk,sk,n,xk);
    vector_copiar(fk,n,yk);
    F(xk,n,fk);
    vector_suma(fk,yk,n,yk);
    //Actualizar A_k 
    matriz_vector_mul(Ak,sk,n,n,auxv);
    vector_resta(yk,auxv,n,auxv);
    vector_escalar(punto(sk,sk,n),auxv,n,auxv);
    vector_vector_mul(auxv,sk,n,auxm);
    matriz_suma(Ak,auxm,n,n,Ak);
    printf("%d ",k);
    for(int i=0;i<n;i++) printf("%lf ",xk[i]);
    printf("%lf %lf\n",nf,kappa);
  }
  vector_copiar(xk,n,out);
  free(auxv);
  free(sk);
  free(yk);
  free(xk); 
  free(fk); 
  liberar_matriz(Ak,n);
  liberar_matriz(auxm,n);
return(1);}


int BFGS(double(*f)(double*,int),int(*g)(double*,int,double*),double **H0, double *x0, 
    int n, int maxiter, double tol, double *out){
  double *xk=crear_vector(n);  
  double *xkm1=crear_vector(n);  
  double *gk=crear_vector(n);  
  double *gkm1=crear_vector(n); 
  double *pk=crear_vector(n);
  double **Hk=crear_matriz(n,n);
  double *yk=crear_vector(n);
  double *sk=crear_vector(n);
  double fk, fkm1;
  double normg,alp,rhok,denrho; 
  double **ss=crear_matriz(n,n);
  double **sy=crear_matriz(n,n);
  double **ys=crear_matriz(n,n);
  double **Iys=crear_matriz(n,n);
  double **Isy=crear_matriz(n,n);
  double **Haux=crear_matriz(n,n);
  vector_copiar(x0,n,xk);
  matriz_copiar(H0,n,n,Hk);
  g(xk,n,gk);
  for(int k=0;k<maxiter;++k){
    normg=Norma_2_vector(gk,n);
    if(normg<tol) break; 
    
    fk=f(xk,n);
    printf("%1.4lf ",fk);
    matriz_vector_mul(Hk,gk,n,n,pk);
    vector_escalar(-1.0,pk,n,pk);
    alp=backtracking(f,xk,fk,gk,pk,n);
    vector_zapyx(xk,alp,pk,n,xkm1);

    g(xkm1,n,gkm1);
    vector_resta(xkm1,xk,n,sk);
    vector_resta(gkm1,gk,n,yk);
    denrho=punto(yk,sk,n);
    if(denrho>1e-4){
      rhok=1.0/denrho;
      vector_vector_mul(sk,sk,n,ss);matriz_escalar(rhok,ss,n,n,ss);
      vector_vector_mul(sk,yk,n,sy);matriz_escalar(rhok,ss,n,n,sy);
      vector_vector_mul(yk,sk,n,ys);matriz_escalar(rhok,ss,n,n,ys);
      matriz_identidad(n,n,Iys);
      matriz_resta(Iys,ys,n,n,Iys);
      matriz_identidad(n,n,Isy);
      matriz_resta(Isy,sy,n,n,Isy);
      matriz_mul(Isy,Hk,n,n,n,Haux);
      matriz_mul(Haux,Iys,n,n,n,Hk);
      matriz_suma(Hk,ss,n,n,Hk);
    }
    vector_copiar(xkm1,n,xk);
    vector_copiar(gkm1,n,gk);
  }
  vector_copiar(xk,n,out);
  free(xk); free(xkm1); 
  free(gk); free(gkm1); 
  free(pk); free(yk); free(sk);
  liberar_matriz(Hk,n);
  liberar_matriz(ss,n);
  liberar_matriz(sy,n);
  liberar_matriz(ys,n);
  liberar_matriz(Isy,n);
  liberar_matriz(Iys,n);
  liberar_matriz(Haux,n);
return(1);}

double backtracking(double(*f)(double*,int),double *xk, double fk, double *gk,double *pk, int n){
  double alpha=1.0, c=1e-4, rho=0.5; 
  double fn, aux=c*punto(gk,pk,n); 
  double *nx=crear_vector(n);
  double *np=crear_vector(n);
  int count=0; 
  do{
    count++;
    alpha=rho*alpha; 
    vector_escalar(alpha,pk,n,np);
    vector_suma(xk,np,n,nx);
    fn=f(nx,n);
  }while(fn>(fk+alpha*aux) );
  free(nx); 
  free(np);
return(alpha);}



/*
* Aproximacion del hessiano por medio de diferencias finitas
*/
int aproximaHessiano(double(*f)(double*,int),double *x0, int n, double h, double **out){
  double *xj=crear_vector(n);
  double fx=f(x0,n);
  double fj,fi,fij;
  vector_copiar(x0,n,xj);
  for(int i=0;i<n;++i){
    x0[i]=x0[i]+h;
    fi=f(x0,n);
    for(int j=0;j<n;++j){
      x0[j]=x0[j]+h;
      xj[j]=xj[j]+h; 
      fj=f(xj,n);
      fij=f(x0,n);
      out[i][j]=(fij-fi-fj+fx)/(h*h);
      x0[j]=x0[j]-h;
      xj[j]=xj[j]-h;
    }
    x0[i]=x0[i]-h;
  }
  free(xj);  
return(1);}

//Funciones de prueba
int funt08(double *x, int n, double *f){
  f[0]=3.0*x[0]-cos(x[1]*x[2])-0.5;
  f[1]=x[0]*x[0]-81.0*(x[1]+0.1)*(x[1]+0.1)+sin(x[2])+1.06;
  f[2]=exp(-1.0*x[0]*x[1])+20.0*x[2]+(10.0*M_PI-3.0)/3.0;
}

int jact08(double *x, int n, double **j){
  j[0][0]=3.0;
  j[0][1]=x[2]*sin(x[1]*x[2]);
  j[0][2]=x[1]*sin(x[1]*x[2]);

  j[1][0]=2.0*x[0];
  j[1][1]=-162.0*(x[1]+0.1);
  j[1][2]=cos(x[2]);
  
  j[2][0]=-x[1]*exp(x[0]*x[1]);
  j[2][1]=-x[0]*exp(x[0]*x[1]);
  j[2][2]=20.0;
}

