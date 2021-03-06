#ifndef LECTURA_H
#define LECTURA_H 
void readParams(int argc, char *argv[],char *cfile,int *maxiter, double *tg, double *tx, double *tf,char *msg); 
double *readVector(char *cfile, int *nr);
int writeVector(double *vec, int dim, char *cfile);
void printVector(double *vec, int dim); 
double **createMatrix(int nr, int nc); 
double **readMatrix(char *cfile, int *nr, int *nc); 
int writeMatrix(double **mat, int nr, int nc, char *cfile); 
int writeData(double **mat, int nr, int nc, char *cfile); 
void printMatrix(double **mat, int nr, int nc);
void freeMatrix(double **mat);
/*
 * Lectura y escritura de archivos de texto plano
 * (dat, txt ...)
 * */
int escribirVector(double *vec, int dim, char *cfile);
double *leerVector(char *cfile, int *nr);
double *leeryk(char *cfile, int *nr,int col);
int escribirEjer4(double *vec,double *x, double *y, int dim, char *cfile);
#endif 
