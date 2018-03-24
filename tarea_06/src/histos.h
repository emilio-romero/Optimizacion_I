#ifndef HISTOS_H 
#define HISTOS_H 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "algebralineal.h"
#include "optimizadores06.h"
void permutar(int n, double **c);
int lee_histo(char *hfile, int bins, double *h);
double fhisto(int bins, double *hc,int nfun, double **c,double *alp, double **mu,double var);
int gahisto(int bins, double *hc,int nfun, double **c,double *alp, double **mu,
    double var,double *ga);

int gmhisto(int bins, double *hc,int nfun, double **c,double *alp, double **mu,
    double var,double **gm);

int hahisto(int bins, double *hc,int nfun, double **c,double *alp, double **mu,
    double var,double **ha);

/*============Optimizacion de la combinacion de la base radial====================*/

int histoMRC(int bins, char *hfile,int nfun, double **c,double *alp, double **mu,
    double var,int maxiter, double tolg, double dhat, double d0, double eta);

int pasoMu(double ***hk, double **gk,double dk, int nfun,double **pmk);
#endif 
