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



double *SteepestDescent(double(*f)(datos),int(*g)(double*,int,double*),datos mid, Condiciones micond);


/*Tamanios de paso*/
double paso_fijo(double a);
double paso_hessiano(double *g, double **H, int n);
double paso_aproximacion(double f, double fk, double ak, double *g, int n);

#endif 
