#ifndef OPTIMIZACION_H 
#define OPTIMIZACION_H 
#include "algebralineal.h" 
#include <stdlib.h>
#include <stdio.h>
int NLSNewton(int(*F)(double *x, int n, double *out),int(*J)(double *x, int n, double **out),
    double *x0, int n,int maxiter, double tol);

int NLSBroyden(int(*F)(double *x, int n, double *out),double **A0,
    double *x0, int n,int maxiter, double tol);

int funt08(double *x, int n, double *f);
int jact08(double *x, int n, double **j);
#endif 
