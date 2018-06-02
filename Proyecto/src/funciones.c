#include "funciones.h"

double rosenbrock(gsl_vector *x){
  double aux=0,xip1, xi;
  int n=(int)x->size; 
  for(int i=0;i<n-1;++i){
    xip1=gsl_vector_get(x,i+1);
    xi=gsl_vector_get(x,1);
    aux=aux+100.0*gsl_pow_2(xip1-xi*xi) + gsl_pow_2(1.0-xi); 
  }
return(aux);}

int grosenbrock(gsl_vector *x, gsl_vector *g){
  int n=(int)x->size;
  double xip1,xi,xim1;
  double t1,t2,t3; 
  xi=gsl_vector_get(x,0);
  xip1=gsl_vector_get(x,1);
  gsl_vector_set(g,0,-400.0*xi*(xip1-xi*xi)-2.0*(1-xi));
  for(int i=1;i<n-1;++i){
    xim1=gsl_vector_get(x,i-1);
    xi=gsl_vector_get(x,i);
    xip1=gsl_vector_get(x,i+1);
    t1=200.0*(xi-xim1*xim1);
    t2=-400.0*xi*(xip1-xi*xi);
    t3=-2.0*(1.0-xi);
    gsl_vector_set(g,i,t1+t2+t3);
  }
  xim1=gsl_vector_get(x,n-2);
  xi=gsl_vector_get(x,n-1);
  gsl_vector_set(g,n-1,200.0*(xi-xim1*xim1));
return(1);}

double sphere(gsl_vector *x){
  double aux=0,xi; 
  int n=(int)x->size; 
  for(int i=0;i<n;++i){
   xi=gsl_vector_get(x,i);
   aux+=xi*xi;
  }
return(aux);}

int gsphere(gsl_vector *x, gsl_vector *g){
int n=(int)x->size; 
double xi; 
for(int i=0;i<n;i++){
  xi=gsl_vector_get(x,i);
  gsl_vector_set(g,i,2.0*xi);
}
return(1);}

double rastrigin(gsl_vector *x){
  double aux=0,xi; 
  int n=(int)x->size; 
  for(int i=0;i<n;++i){
    xi=gsl_vector_get(x,i);
    aux+=xi*xi-10.0*gsl_sf_cos(2.0*xi*M_PI);
  }
return(aux+10*n);}

int grastrigin(gsl_vector *x, gsl_vector *g){
int n=(int)x->size; 
double xi; 
for(int i=0;i<n;i++){
  xi=gsl_vector_get(x,i);
  gsl_vector_set(g,i,2.0*xi+20.0*M_PI*gsl_sf_sin(2.0*xi*M_PI));
}
return(1);}


