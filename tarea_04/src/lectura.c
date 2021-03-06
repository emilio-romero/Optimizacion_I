#include <stdio.h>
#include <stdlib.h> 
#include <string.h> 
void readParams(int argc, char *argv[],char *cfile, int *maxiter, double *tg, double *tx, double *tf,char *msg){
if(argc>1) strcpy(cfile,argv[1]);
if(argc>2) *maxiter=atoi(argv[2]);
if(argc>3) *tg=atof(argv[3]);
if(argc>4) *tx=atof(argv[4]);
if(argc>5) *tf=atof(argv[5]);
if(argc>6) strcpy(msg,argv[6]);
}

//Lectura del vector en el archivo cfile 
//Devuelve NULL si no puede abrir el archivo 
double *readVector(char *cfile, int *nr){
double *vec; 
FILE *f1=fopen(cfile,"rb");
if(!f1) return(NULL); 
fread(nr,sizeof(int),1,f1); 
vec=(double*)malloc((*nr)*sizeof(double));
if(vec==NULL) return(NULL); 
fread(vec,sizeof(double),*nr,f1); 
fclose(f1); 
return(vec); 
}

//Escritura del vector vec en el archivo cfile
//devuelve 0 en caso de exito y 1 si no
int writeVector(double *vec, int dim, char *cfile){
FILE *f1=fopen(cfile,"wb"); 

if(!f1) return(1); 
fwrite(&dim, sizeof(int),1,f1);
fwrite(vec,sizeof(double),dim,f1);
fclose(f1); 
return(0);  
}

//imprime en la consola las entradas del vector 
void printVector(double *vec, int dim){
int i; 
for(i=0;i<dim;++i)
  printf("%8.4f  ",vec[i]);
printf("\n");
}


//Reserva memoria para una matriz de tamanio nr x nc 
double **createMatrix(int nr, int nc){
int i; 
double **mat; 

mat=(double **) malloc((nr)*sizeof(double*));
if(mat==NULL) return(NULL); 
  mat[0]=(double*)malloc(nr*nc*sizeof(double));
if(mat[0]==NULL) return(NULL); 
for(i=1;i<nr;++i)
  mat[i]=mat[i-1]+nc;
return(mat);}

//Liberar memoria para las matricies 
void freeMatrix(double **mat){
free(mat[0]); 
free(mat); 
}

//Lectura de la matriz en el archivo cfile 
//Devuelve NULL si no se pudo abrir el archivo 
double **readMatrix(char *cfile,int *nr, int *nc){
double **mat; 
FILE *f1=fopen(cfile,"rb"); 
if(!f1) return(NULL); 
fread(nr,sizeof(int),1,f1);
fread(nc,sizeof(int),1,f1);
//Reservamos memoria
mat=createMatrix(*nr,*nc);
//Lectura de los datos 
fread(mat[0],sizeof(double),(*nr)*(*nc),f1); 
fclose(f1); 
return(mat); 
}
//Escritura de la matriz en el archivo cfile 
//Devuelve 0 en caso de extio y 1 si no 
int writeMatrix(double **mat, int nr, int nc, char *cfile){
FILE *f1=fopen(cfile,"wb"); 

if(!f1) return(1);
fwrite(&nr,sizeof(int),1,f1);
fwrite(&nc,sizeof(int),1,f1);
fwrite(mat[0],sizeof(double),nr*nc,f1);
fclose(f1);
return(0);
}

void printMatrix(double **mat, int nr, int nc){
int i,j;
for(i=0;i<nr;i++){
  for(j=0;j<nc;j++)
    printf("%6.2f    ",mat[i][j]);
  printf("\n");
}
}

int writeData(double **mat, int nr, int nc, char *cfile){
FILE *f1=fopen(cfile,"w");
if(!f1) return(0);
for(int i=0;i<nr;i++){
 for(int j=0;j<nc;j++){
  fprintf(f1,"%lf ",mat[i][j]); 
 }
 fprintf(f1,"\n");
}
fclose(f1);
return(1);}


/*
 * Funciones para archivos de texto plano
 *
 */
int escribirVector(double *vec, int dim, char *cfile){
FILE *f1=fopen(cfile,"w"); 

if(!f1) return(1); 
fprintf(f1,"%d\n",dim);
for(int i=0;i<dim;i++)
  fprintf(f1,"%lf\n",vec[i]);

fclose(f1); 
return(0);
}

double *leerVector(char *cfile, int *nr){
double *vec; 
FILE *f1=fopen(cfile,"r");
if(!f1) return(NULL); 
fscanf(f1,"%d",nr); //printf("%d\n",*nr); 
vec=(double*)malloc((*nr)*sizeof(double));
//if(vec==NULL) return(NULL); 
for(int i=0;i<(*nr);i++){
  fscanf(f1,"%lf",vec+i);
 // printf("%lf\n",vec[i]);
  }
fclose(f1); 
return(vec);
}

double *leeryk(char *cfile, int *nr, int col){
  double *vec; 
FILE *f1=fopen(cfile,"r");
if(!f1) return(NULL); 
fscanf(f1,"%d",nr); //printf("%d\n",*nr); 
vec=(double*)malloc((*nr)*sizeof(double));


for(int i=0;i<(*nr);i++){
  fscanf(f1,"%lf",vec+i);
}
if(col==1){ 
  fclose(f1);
  return(vec);} 
if(col==2){
  for(int i=0;i<(*nr);i++){
    fscanf(f1,"%lf",vec+i);
  }
} 
fclose(f1);
return(vec);
}


int escribirEjer4(double *vec,double *x, double *y, int dim, char *cfile){
FILE *f1=fopen(cfile,"w"); 

if(!f1) return(1); 
//fprintf(f1,"%d\n",dim);
for(int i=0;i<dim;i++)
  fprintf(f1,"%d %lf %lf %lf\n",i,vec[i],x[i],y[i]);

fclose(f1); 
return(0);
}




