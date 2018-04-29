#ifndef BENCHMARK_H 
#define BENCHMARK_H 
#include "algebralineal.h"

double Rosenbrock(double *x, int n); 
int gRosenbrock(double *x, int n, double *out);
int hRosenbrock(double *x, int n, double **out);
#endif 
