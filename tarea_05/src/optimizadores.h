#ifndef OPTIMIZADORES_H 
#define OPTIMIZADORES_H 
#include "algebralineal.h"
#include "funciones.h"
#include <string.h>
/*Estructura para condiciones de paro*/
typedef struct{
  int maxiter; 
  double tolg; 
  double tolx; 
  double tolf;
  char msg[15];
}Condiciones; 



double *SteepestDescent(double(*f)(datos),int(*g)(datos,double*),int(*h)(datos,double**),datos mid, Condiciones micond);

double *SteepestDescent2(double(*f)(datos),int(*g)(datos,double*),datos x,Condiciones mc,int *flag); 

/*Tamanios de paso*/
double paso_fijo(double a);
double paso_hessiano(double *g, double **H, int n);
double paso_aproximacion(double f, double fk, double ak, double *g, int n);
double backtracking(double(*f)(datos),datos xk,double fk,double *gk, int *flag);
double quadInterpolation(double(*f)(datos),datos xk,double fk,double *gk,double a0);
double cubicInterpolation(double(*f)(datos),datos xk,double fk,double *gk,double pa,double sa);
#endif 
