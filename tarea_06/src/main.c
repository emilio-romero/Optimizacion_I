#include <stdio.h>
#include "optimizadores06.h"
#include "histos.h"
#include <time.h>
int main(int argc, char *argv[]){
srand(5);
char histoa[50]; 
if(argc>1) strcpy(histoa,argv[1]);
Condiciones myc; 
myc.maxiter=10000; 
myc.tolg=1e-8; 
datos x0; 
x0.n=100; 
x0.x=(double*)malloc(x0.n*sizeof(double));

for(int i=0;i<x0.n;i++) x0.x[i]=1.0; 
x0.x[0]=-1.2; x0.x[x0.n-2]=-1.2; 
//x0.x[0]=-1.2; x0.x[1]=1.0; 
int nfun=5; 
int nbin=3;
double **mc=crear_matriz(27,3);

permutar(nbin,mc);
//for(int i=0;i<27;i++){
//  printf("%lf %lf %lf\n",mc[i][0],mc[i][1],mc[i][2]);
//}
double **mmu=crear_matriz(nfun,3);
double *ma=crear_vector(nfun);
double *hc=crear_vector(nbin*nbin*nbin);
lee_histo(histoa,nbin,hc);
for(int i=0;i<nfun;i++){
  for(int j=0;j<3;j++){
    mmu[i][j]=2.0*randx(); 
//    printf("%lf ",mmu[i][j]);
  }
//  printf("\n");
  ma[i]=1.0/nfun;
}
//printf("alo\n");
double *gah=crear_vector(nfun);
double **gmh=crear_matriz(nfun,3);
double fh=fhisto(3,hc,nfun,mc,ma,mmu,1.0);
double **hah=crear_matriz(nfun,nfun);
gahisto(3,hc,nfun,mc,ma,mmu,1.0,gah);
gmhisto(3,hc,nfun,mc,ma,mmu,1.0,gmh);
hahisto(3,hc,nfun,mc,ma,mmu,1.0,hah);

printf("Evaluacion de g(a,m): %lf\n",fh);
printf("Gradiente respecto a alfa:\n");
for(int i=0;i<nfun;i++) printf("%lf \n",gah[i]);
printf("\n");
printf("Gradiente respecto a mu:\n");
for(int i=0;i<nfun;i++) printf("%lf %lf %lf \n",gmh[i][0],gmh[i][1],gmh[i][2]);
printf("\n");

printf("hess respecto a alfa:\n");
for(int i=0;i<nfun;i++){
  for(int j=0;j<nfun;j++)
    printf("%lf ",hah[i][j]);
printf("\n");}

histoMRC(nbin,histoa,nfun,mc,ma,mmu,1,1000,1e-6,4.0,1.5,0.17);


/*double *opti=MRC(Rosenbrock,gRosenbrock, hRosenbrock,x0,myc,1.0,0.085,0.22);
for(int i=0;i<10;i++) printf("%lf ",opti[i]);
printf("\n");*/

free(ma);
free(hc);
free(gah);
liberar_matriz(mc,27);
liberar_matriz(mmu,nfun);
liberar_matriz(gmh,nfun);
liberar_matriz(hah,nfun);
free(x0.x);
  printf("Su programa ha terminado\n");
return 0;}
