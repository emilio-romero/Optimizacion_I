#ifndef ALGLINEAL_H
#define ALGLINEAL_H 
#include <math.h> 
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "lectura.h"
//=========================
//= Operaciones elementales
//========================

double punto(double *a, double *b, int n);
int matriz_vector_mul(double **A, double *b, int m, int n, double *out);
int matriz_mul(double **A, double **B, int l, int m, int n, double **out);

//=========================
//= Normas 
//========================
int Norma_1_matriz(double **A, int nr, int nc, double *out);

//=========================
//= Soluciondores  
//========================

int Cholesky(double **A, int n, double **out);



//=========================
//= Miscelanea 
//========================

int es_spd(char *cfile); 

#endif 
