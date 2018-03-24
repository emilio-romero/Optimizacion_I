#ifndef ALGLINEAL_H
#define ALGLINEAL_H 
#include <math.h> 
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "lectura.h"
//#include "solucionadores.h"

/*===============================================
 * Creacion y liberacion de matrices y vectores
 *===============================================
 */
double *crear_vector(int n);
double **crear_matriz(int nr, int nc);

//int liberar_vector(double *x);
int liberar_matriz(double **x, int nr);


//=========================
//= Operaciones elementales
//========================
/*
 * Operaciones con vectores 
 */
double punto(double *a, double *b, int n);
int vector_suma(double *x, double *y, int n, double *out);
int vector_resta(double *x, double *y, int n, double *out);
int vector_escalar(double a, double *x,int n,double *out);
int vector_copiar(double *original, int n, double *copia);
/*
 * Operaciones con matrices
 * */
int matriz_vector_mul(double **A, double *b, int m, int n, double *out);
int matriz_mul(double **A, double **B, int l, int m, int n, double **out);
int matriz_suma(double **A, double **B, int nr, int nc, double **out);
int matriz_resta(double **A, double **B, int nr, int nc, double **out);
int matriz_copiar(double **original, int nr, int nc, double **copia);
int matriz_transponer(double **original, int nr, int nc, double **copia);
//=========================
//= Normas 
//========================
int Norma_1_matriz(double **A, int nr, int nc, double *out);

/*Normas vectoriales*/
double Norma_1_vector(double *x, int n);
double Norma_2_vector(double *x, int n);
double Norma_inf_vector(double *x, int n);

//=========================
//= Soluciondores  
//========================

/*
 * Solucionadores basicos 
 */

int sinferior(double **L, double *b, int n, double *out);
int ssuperior(double **U, double *b, int n, double *out);
/*Solucionadores con factoriazaciones */
int factoLU(double **A, int n, double **L, double **U);
int solLU(double **A, double *n, int nr, int nc, double *out);

int Chol(double **A, int n, double **out);
int Cholesky(char *cfile, double**out);
int solLL(double **A,double *b, int nr, int nc, double *out);

//=========================
//= Minimos cuadrados
//========================
//double *aproximaPolinomio(int grado, double **data, int npuntos);



//=========================
//= Miscelanea 
//========================

int reduccionMatriz(double **A,double *a,int nc, int nr, int nnf, double **aux, double *auxv);

int es_spd(char *cfile); 
double randx();
double minimo2(double a, double b);
#endif 
