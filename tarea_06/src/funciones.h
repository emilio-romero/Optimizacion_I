#ifndef FUNCIONES_H 
#define FUNCIONES_H
#include "algebralineal.h"

typedef struct{
  int n; 
  int obs; 
  double *x; 
  double *y;
  double **dx;
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

/*Tarea 05*/
double lLogistic(datos x0);
int glLogistic(datos x0,double *g);
double errorLogistico(double *be, double **x,double *y, int nc, int obs);
#endif 
