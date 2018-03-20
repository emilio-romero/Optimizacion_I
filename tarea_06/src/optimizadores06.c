#include "optimizadores06.h"

double *MRC(double(*f)(datos),int(*g)(datos,double*),int(*h)(datos,double**)\
            ,datos x0, Condiciones mc,double dhat,double d0,double eta){
  datos xk; 
  datos xkm1;
  /*Area de copiado para las variables auxiliares*/
   xk.n=x0.n; xkm1.n=x0.n;
   xk.x=(double*)malloc(xk.n*sizeof(double));
   xkm1.x=(double*)malloc(xkm1.n*sizeof(double));
   vector_copiar(x0.x,xk.n,xk.x);
   /*Fin de copiado*/
  double fk,fkm1,dk=d0,dkm1,rhok; 
  double *gk=(double*)malloc(x0.n*sizeof(double));
  double *pck=(double*)malloc(xk.n*sizeof(double));
  double *pbk=(double*)malloc(xk.n*sizeof(double));
  double *pk=(double*)malloc(xk.n*sizeof(double));
  double **hk=(double**)malloc(xk.n*sizeof(double*));
  for(int i=0;i<xk.n;i++) hk[i]=(double*)malloc(xk.n*sizeof(double));
  for(int k=0;k<mc.maxiter;k++){
//printf(" %g %g\n",xk.x[0],xk.x[1]);
   fk=f(xk); 
   g(xk,gk);  
   h(xk,hk);
   //Calculo del avance dada la region de confianza
    pCauchy(hk,gk,dk,xk.n,pck);
    pNewton(hk,gk,xk.n,pbk);
    pDogleg(pck,pbk,dk,xk.n,pk);

   vector_suma(xk.x,pk,xkm1.n,xkm1.x);
   fkm1=f(xkm1);
   //Calculo de la calidad del modelos (yo no lo entiendos xD)
   rhok=calidadModelo(fk,fkm1,pk,gk,hk,xk.n); 
   // Actualizacion de xk
   if(rhok>eta){
    vector_suma(xk.x,pk,xkm1.n,xkm1.x);
   }
   else{
    vector_copiar(xk.x,xkm1.n,xkm1.x);
   }

   //Actualizacion de Delta_k
   if(rhok<0.25){
    dk=0.25*dk;
   }else{
    if(rhok>0.75 && fabs(Norma_2_vector(pk,xk.n)-dk)<1e-6){
      dk=minimo2(2.0*dk,dhat);
    }else{
      dk=dk;
    }
   }

//Criterios de paro
  if(Norma_2_vector(gk,xk.n)<mc.tolg) break; 
  vector_copiar(xkm1.x,xk.n,xk.x);  
  }



//Liberacion de memoria
  free(gk); free(pk);
  for(int i=0;i<xk.n;i++) free(hk[i]);
  free(hk);

  return(xkm1.x);
}

int pCauchy(double **hk,double *gk, double delk,int ndim, double *pc){
  double *bg=(double*)malloc(ndim*sizeof(double));
  double gbg;
  double tauk; 
  double eoco2,ng;
  matriz_vector_mul(hk,gk,ndim,ndim,bg);
  gbg=punto(gk,bg,ndim);
  ng=Norma_2_vector(gk,ndim);
  if(gbg<=0){ //Calculo de \tau_k
    tauk=1.0;
    tauk=-tauk*ng/delk;
    vector_escalar(tauk,gk,ndim,bg);
    vector_copiar(bg,ndim,pc);
  } 
  else{
     eoco2=ng*ng*ng/(delk*gbg);
     if(eoco2<1.0){
      tauk=eoco2;
      tauk=-tauk*ng/delk;
      vector_escalar(tauk,gk,ndim,bg);
      vector_copiar(bg,ndim,pc);
     }else{
      tauk=1.0; 
      tauk=-tauk*ng/delk;
      vector_escalar(tauk,gk,ndim,bg);
      vector_copiar(bg,ndim,pc);
     }
  }


  free(bg);
return(1);}

int pNewton(double **hk, double *gk,int ndim, double *pb){
  double *ng=(double*)malloc(ndim*sizeof(double));
  vector_escalar(-1.0,gk,ndim,ng);
  //Resolver el sistema 
  solLU(hk,ng,ndim,ndim,pb); 
return(1);}

int pDogleg(double *pc, double *pb, double dk, int ndim, double *out){
  double nc, nb, naux; 
  double a,b,c,d,alpha;
  double *aux=(double*)malloc(ndim*sizeof(double));
  nc=Norma_2_vector(pc,ndim);
  nb=Norma_2_vector(pb,ndim);
  if(dk<=nc){
    vector_escalar(dk/nc,pc,ndim,out);
  }
  else if(dk>=nb){
    vector_copiar(pb,ndim,out);
  }
  else{
    a=nc*nc; 
    b=nb*nb;
    vector_resta(pc,pb,ndim,aux);
    naux=Norma_2_vector(aux,ndim);
    naux=naux*naux; 
    d=(a+b-c)/2; 
    alpha=(b-dk*dk)/(b-d+sqrt(d*d-a*b+dk*dk*c));
    vector_escalar(alpha,pc,ndim,out); 
    vector_escalar((1.0-alpha),pb,ndim,aux);
    vector_suma(aux,out,ndim,out);
  }

free(aux);
return(1);}

double calidadModelo(double fk, double fkm1, double *pk,double *gk, double **hk,int ndim){
  double aux; 
  double nume=fk-fkm1; 
  double deno,pbp;
  double *bp=(double*)malloc(ndim*sizeof(double));
  matriz_vector_mul(hk,pk,ndim,ndim,bp);
  pbp=punto(pk,bp,ndim);
  deno=-1.0*(punto(gk,pk,ndim)+0.5*pbp);
  free(bp);
  aux=nume/deno; 
  return(aux);
}
