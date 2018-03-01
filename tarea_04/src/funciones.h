#ifndef FUNCIONES_H 
#define FUNCIONES_H
#include "algebralineal.h"

typedef struct{
  int n; 
  double *x; 
  double param1; 
  double param2;
} datos; 


double Rosenbrock(datos md);
int gRosenbrock(double *x, int n, double *out);
int hRosenbrock(double *x, int n, double **out);
double Wood(datos md);
int gWood(double *x, int n, double *out);
int hWood(double *x, int n, double **out);
double SmoothingModel(datos md);

#endif 
