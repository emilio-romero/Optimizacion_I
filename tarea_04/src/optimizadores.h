#ifndef OPTIMIZADORES_H 
#define OPTIMIZADORES_H 
#include "algebralineal.h"
#include "funciones.h"

/*Estructura para condiciones de paro*/
typedef struct{
  int maxiter; 
  double tolg; 
  double tolx; 
  double tolf;
}Condiciones; 



double *SteepestDescent(double(*f)(double*,int),int(*g)(double*,int,double*),double *x0, int n, Condiciones micond);


/*Tamanios de paso*/
double paso_fijo(double a);
double paso_hessiano();
double paso_aproximacion();

#endif 
