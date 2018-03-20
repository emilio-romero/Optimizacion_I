#ifndef OPTIMIZADORES06_H
#define OPTIMIZADORES06_H
#include "algebralineal.h"
#include "funciones.h"
#include <string.h>


double *MRC(double(*f)(datos),int(*g)(datos,double*),int(*h)(datos,double**)\
            ,datos x0, Condiciones mc,double dhat, double d0, double eta);

int pCauchy(double **hk, double *gk,double delk,int ndim,double *pc);
int pNewton(double **hk,double *gk,int ndim,double *pb);
int pDogleg(double *pc, double *pb, double dk, int ndim,double *out);

double calidadModelo(double fk,double fkm1, double *pk, double *gk, double **hk, int ndim);
#endif 
