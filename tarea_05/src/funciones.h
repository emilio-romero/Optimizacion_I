#ifndef FUNCIONES_H 
#define FUNCIONES_H
#include "algebralineal.h"

typedef struct{
  int n; 
  double *x; 
  double *y;
  double param1; 
  double param2;
} datos; 


double Rosenbrock(datos md);
int gRosenbrock(datos md, double *out);
int hRosenbrock(datos md, double **out);
double Wood(datos md);
int gWood(datos md, double *out);
int hWood(datos md, double **out);
double SmoothingModel(datos md);
int gSModel(datos md, double *out);
int hSModel(datos md, double **out);
#endif 
