#ifndef OPTIMIZACION_H 
#define OPTIMIZACION_H 
#include "algebralineal.h" 
#include <stdlib.h>
#include <stdio.h>
int NLSNewton(int(*F)(double *x, int n, double *out),int(*J)(double *x, int n, double **out),
    double *x0, int n,int maxiter, double tol,double *out);

int NLSBroyden(int(*F)(double *x, int n, double *out),double **A0,
    double *x0, int n,int maxiter, double tol,double *out);
int BFGS(double(*f)(double*,int), int(*g)(double*,int,double*),double **H0,double *x0,
    int n, int maxiter, double tol,double *out);
double backtracking(double(*f)(double*,int), double *xk, double fk, double *gk,double *pk,int n);

int aproximaHessiano(double(*f)(double*,int),double *x0,int n,double h, double **out);

int funt08(double *x, int n, double *f);
int jact08(double *x, int n, double **j);
#endif 
