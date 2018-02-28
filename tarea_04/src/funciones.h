#ifndef FUNCIONES_H 
#define FUNCIONES_H
#include "algebralineal.h"

double Rosenbrock(double *x, int n);
int gRosenbrock(double *x, int n, double *out);

double Wood(double *x, int n);
double SmoothingModel(double *x, int n,double lambda);

#endif 
