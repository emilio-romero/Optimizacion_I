#include "algebralineal.h"
/*
* Producto punto de dos vectores de tamanio n 
*/
double punto(double *a, double *b, int n){
  double suma=0; 
  for(int i=0;i<n;i++){
    suma+=a[i]*b[i]; 
  }
return(suma);
}

/*
* Producto de una matriz por un vector, devuelve un vector 
* m es para filas de la matriz y n para las columnas de la matriz
* n tambien son las filas del vector
*/
int matriz_vector_mul(double **A, double *b, int m, int n, double *out){
  double suma=0; 
  for(int i=0;i<m;i++){
    suma=0;
    for(int j=0;j<n;j++){
      suma+=A[i][j]*b[j]; 
    }
    out[i]=suma;
  }
return 1;}

/*
* Producto de matrices, devuelve una matriz 
* l y m son las filas y columnas de la primera matriz
* m y n son las filas y columnas de la segunda matriz
*/

int matriz_mul(double **A, double **B, int l, int m, int n, double **out){
  double suma=0; 

  for(int i=0;i<l;i++){
    for(int j=0;j<n;j++){
      suma=0; 
      for(int k=0;k<m;k++){
        suma+=A[i][k]*B[k][j];
      }
      out[i][j]=suma; 
    }
  }

return 1;}
/*
* Suma de matrices
*/
int matriz_suma(double **A, double **B, int nr, int nc, double **out){
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      out[i][j]=A[i][j]+B[i][j]; 
    }
  }
return(1);}
/*
* Resta de matrices
*/
int matriz_resta(double **A, double **B, int nr, int nc, double **out){
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      out[i][j]=A[i][j]-B[i][j]; 
    }
  }
return(1);}
/*
* Copia una matriz A a otra B
*/
int matriz_copiar(double **original, int nr, int nc, double **copia){
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      copia[i][j]=original[i][j];
    }
  }
return(1);}
/*
* Copia una matriz A a otra B
*/
int matriz_transponer(double **original, int nr, int nc, double **copia){
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      copia[i][j]=original[j][i];
    }
  }
return(1);}


/*
* Calcula la norma 1 de una matriz 
*/
int Norma_1_matriz(double **A, int nr, int nc, double *out){
  double max=0; 
  double aux=0; 
  for(int i=0;i<nr;i++) 
    max+=fabs(A[i][0]); 
  
  for(int j=1;j<nc;j++){
    for(int i=0;i<nr;i++) aux+=fabs(A[i][j]); 
    if(aux>max){
      max=aux; 
    }
    aux=0; 
  }
*out=max; 
}

int es_spd(char *cfile){
  int nr, nc; 
  double **A=readMatrix(cfile, &nr, &nc);
  double **AA=(double**)malloc(nr*sizeof(double*)); 
  double N1;
  double meps=sqrt(DBL_EPSILON);
  for(int i=0;i<nr;i++) AA[i]=(double*)malloc(nc*sizeof(double));
    
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      AA[i][j]=A[i][j]-A[j][i];
    }
  }
  
  Norma_1_matriz(AA,nr,nc,&N1);
  //printf("%g\n",N1);
  if(N1>meps){
    printf("La matriz dada no es simetrica\n");
    printf("Se debe utilizar la matriz (A+A^T)/2\n");
    for(int i=0;i<nr;i++) free(AA[i]); 
    free(AA);
    freeMatrix(A);
    return(0);
  } else{
    printf("La matriz es simetrica, hurra!\n");
    for(int i=0;i<nr;i++) free(AA[i]); 
    free(AA);
    freeMatrix(A);
    return(1);
  }
  /*Liberacion de memoria*/
return 1;}

/*
* Factorizacion Cholesky
* esta funcion devuelve la matriz L de la factorizacion
*/

int Chol(double **A, int n, double **out){
  double restakj=0; 
  double restaikj=0; 
  double meps=sqrt(DBL_EPSILON);
  out[0][0]=sqrt(A[0][0]);
  if(out[0][0]<meps){
    printf("No se ha podido factorizar\n");
    return(0);
  }
  for(int i=1;i<n;i++){
    out[i][0]=(A[i][0])/out[0][0];
  }
  for(int j=1;j<n;j++){
    for(int k=0;k<j;k++){
      restakj+=out[j][k]*out[j][k]; 
    }
    out[j][j]=sqrt(A[j][j]-restakj);
    if(out[j][j]<meps || (A[j][j]-restakj)<0){
      //printf("No se ha podido factorizar\n");
      return(0);
    }
    restakj=0;
    for(int i=j+1;i<n;i++){
      for(int k=0;k<j;k++){
        restaikj+=out[i][k]*out[j][k]; 
      }
      out[i][j]=(A[i][j]-restaikj)/out[j][j];
      restaikj=0; 
    }
  }

return(1);}

int Cholesky(char *cfile, double **out){
  int simetria=es_spd(cfile); 
  int nr, nc; 
  int df, contador=0;
  double error;
  double **A=readMatrix(cfile, &nr, &nc);
  double **AA=(double**)malloc(nr*sizeof(double*)); 
  double **LL=(double**)malloc(nr*sizeof(double*)); 
  for(int i=0;i<nr;i++){ AA[i]=(double*)malloc(nc*sizeof(double));
     LL[i]=(double*)malloc(nc*sizeof(double));}
  double meps=sqrt(DBL_EPSILON);

  if(simetria==0){
    for(int i=0;i<nr;i++){
      for(int j=0;j<nc;j++){
        AA[i][j]=(A[i][j]+A[j][i])/2.0; 
      }
    }
   matriz_copiar(AA,nr,nc,A);
  }

  df=Chol(A,nr,out);
  //printf("%d\n",df);
  while(df==0){
    for(int i=0;i<nr;i++){
      A[i][i]=A[i][i]+0.01; 
    }
    df=Chol(A,nr,out); 
  //printf("%d\n",df);
    contador++; 
  }
  printMatrix(out,nr,nc);
  matriz_transponer(out,nr,nc,AA);
  matriz_mul(out,AA,nr,nc,nr,LL);
  matriz_resta(A,LL,nr,nc,AA);
  Norma_1_matriz(AA,nr,nc,&error);
  printf("La matriz se perturbo %d veces.\n",contador); 
  printf("El error es %g\n",error);
return(1);}



