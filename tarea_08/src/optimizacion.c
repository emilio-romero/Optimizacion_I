#include "optimizacion.h"

int NLSNewton(int(*F)(double *x, int n, double *out),int(*J)(double *x, int n, double **out),
    double *x0, int n, int maxiter, double tol){
  double *xk=crear_vector(n);
  double *sk=crear_vector(n);
  double *fk=crear_vector(n);
  double **jk=crear_matriz(n,n);
  double nf; 
  for(int k=0;k<maxiter;++k){
    F(xk,n,fk); vector_escalar(-1.0,fk,n,fk); 
    nf=Norma_2_vector(fk,n);
    if(nf<tol) break; 
    J(xk,n,jk);
    solLU(jk,fk,n,n,sk);
    vector_suma(xk,sk,n,xk);
  }

  free(sk);
  free(xk); 
  free(fk); 
  liberar_matriz(jk,n);
return(1);}



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

