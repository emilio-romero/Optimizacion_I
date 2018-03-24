#include "histos.h"

void permutar(int bins, double **c){
  for(int i=0, j=0,k=0;i<(bins*bins*bins);i++){
    if(i!=0 && i%(bins*bins)==0){
      k++;
    }
    c[i][0]=(double)k; 
    if(i!=0 && i%bins==0){
      j++; 
      j=j%bins; 
    }
    c[i][1]=(double)j; 
    c[i][2]=(double)(i%bins); 
  }
}


int lee_histo(char *hfile, int bins, double *h){
  int elems=bins*bins*bins; 
  FILE *f1=fopen(hfile,"r");
  char buff[255]; 
  fgets(buff,255,f1);
  for(int i=0;i<elems;i++){
    fscanf(f1,"%lf",&h[i]);
  }
  fclose(f1);
}

double fhisto(int bins, double *hc,int nfun,double **c ,double *alp, double **mu, double var){
  double aux;
  double sumae=0;
  double sumai=0;
  double auxe; 
  int dbins=bins*bins*bins;
  double *vec3=crear_vector(3);
  //double *hc=crear_vector(dbins);
  //lee_histo(hfile,bins,hc);
  for(int om=0;om<dbins;om++){
    //printf("%lf\n",hc[om]); 
    sumai=0;
    for(int i=0;i<nfun;i++){
      vector_resta(c[om],mu[i],3,vec3);
  //printf("obvio\n");
      auxe=punto(vec3,vec3,3);
      //printf("%lf ",auxe);
      auxe=exp(-auxe/(2.0*var));
      sumai+=alp[i]*auxe;
    }
    sumae+=(hc[om]-sumai)*(hc[om]-sumai); 
  }
  free(vec3);
  //free(hc);
  aux=sumae;
  return(aux);}

int gahisto(int bins, double *hc, int nfun, double **c, double *alp, double **mu,
    double var, double *ga){
  int dbins=bins*bins*bins; 
  double sumae=0, sumai=0; 
  double auxe; 
  double *vec3=crear_vector(3);
  for(int om=0;om<dbins;om++){
    sumai=0.0; 
    for(int i=0;i<nfun;i++){
      vector_resta(c[om],mu[i],3,vec3);
      auxe=punto(vec3,vec3,3);
      auxe=exp(-auxe/(2.0*var));
      sumai+=alp[i]*auxe; 
    }
    for(int i=0;i<nfun;i++){
      vector_resta(c[om],mu[i],3,vec3);
      auxe=punto(vec3,vec3,3);
      auxe=exp(-auxe/(2.0*var));
      ga[i]=ga[i]-2.0*(hc[om]-sumai)*auxe; 
    }
  }

  free(vec3);
return(1);}


int gmhisto(int bins, double *hc, int nfun, double **c, double *alp, double **mu,
    double var, double **gm){
  int dbins=bins*bins*bins; 
  double sumai=0; 
  double auxe,escalar; 
  double *vec3=crear_vector(3);
  for(int om=0;om<dbins;om++){
    sumai=0;
    for(int i=0;i<nfun;i++){
      vector_resta(c[om],mu[i],3,vec3);
      auxe=punto(vec3,vec3,3);
      auxe=exp(-auxe/(2.0*var));
      sumai+=alp[i]*auxe; 
    }
    for(int i=0;i<nfun;i++){
      vector_resta(c[om],mu[i],3,vec3);
      auxe=punto(vec3,vec3,3);
      auxe=exp(-auxe/(2.0*var));
      escalar=2.0*(hc[om]-sumai)*alp[i]*auxe/var;
      vector_escalar(escalar,vec3,3,vec3);
      vector_suma(vec3,gm[i],3,gm[i]);
    }
  }

  free(vec3);
return(1);}
int hahisto(int bins, double *hc, int nfun, double **c, double *alp, double **mu,
    double var, double **ha){
  int dbins=bins*bins*bins; 
  double suma=0, auxa,auxb, escalar; 
  double *veca=crear_vector(3);
  double *vecb=crear_vector(3);
  for(int i=0;i<nfun;i++){
    for(int j=0;j<nfun;j++){
      for(int om=0;om<dbins;om++){
        vector_resta(c[om],mu[i],3,veca);  
        vector_resta(c[om],mu[j],3,vecb);
        auxa=punto(veca,veca,3);
        auxb=punto(vecb,vecb,3);
        ha[i][j]=ha[i][j]-2.0*exp(-(auxa+auxb)/(2.0*var)); 
      }
    }
  }
free(veca);
free(vecb);
return(1);}


/*
 *Optimizacion de la funcion
 *
 */
int histoMRC(int bins, char *hfile, int nfun, double **c, double *alp, double **mu,
    double var,int maxiter, double tolg, double dhat, double d0, double eta){
  int dbins=bins*bins*bins; 
  double *hc=crear_vector(dbins); 
  lee_histo(hfile,bins,hc);
  double *gak=crear_vector(nfun);//Gradiente respecto a alfa
  double **hak=crear_matriz(nfun,nfun); //Hessiano respecto a alfa
  double **gmk=crear_matriz(nfun,3);//Gradiente respecto a mu
  double **hmk=crear_matriz(nfun,nfun); //Hessiano respecto a mu
  double nga,ngm,fak,fmk; //Norma gradiente de alfa, norma gradiente de mu, funcion en paso k
  double *auxa=crear_vector(nfun);double *pcak=crear_vector(nfun);double *pbak=crear_vector(nfun);
  double *pak=crear_vector(nfun);
  double dk=d0,fakm1,rhoak,fmkm1; 
  double **auxm=crear_matriz(3,nfun);
  double *pcmk=crear_vector(nfun);double *pbmk=crear_vector(nfun);
  double *pmk=crear_vector(nfun);
  matriz_transponer(mu,3,nfun,auxm);
  for(int k=0;k<maxiter;k++){
    //Procesos con alfa =============================
    fak=fhisto(bins,hc,nfun,c,alp,mu,var); printf("%lf ",fak);
    gahisto(bins,hc,nfun,c,alp,mu,var,gak);
    hahisto(bins,hc,nfun,c,alp,mu,var,hak);
    //Calculo del tamanio de paso para alfa 
    pCauchy(hak,gak,dk,nfun,pcak);pNewton(hak,gak,nfun,pbak);
    pDogleg(pcak,pbak,dk,nfun,pak);
    //printf("%lf %lf ---",pak[0],pak[nfun-1]);
    //==
    vector_suma(alp,pak,nfun,auxa); //reemplazar alp por un auxiliar al 
    fakm1=fhisto(bins,hc,nfun,c,auxa,mu,var); //lo mismo ^
    //Calculo de la calidad del modelo (con alfa)
     rhoak=calidadModelo(fak,fakm1,pak,gak,hak,nfun);
    //Actualizacion de xk
     if(rhoak>eta){
      vector_suma(alp,pak,nfun,auxa);
     }
     else{
      vector_copiar(alp,nfun,auxa);//de nuevo reemplazar por auxiliar al
    }
    vector_copiar(auxa,nfun,alp);//reemplazar por auxiliar al//esto no va :/ ?
    //Fin de procesos con alfa========================
    // Procesos con mu ===============================   
    fmk=fhisto(bins,hc,nfun,c,alp,mu,var);
    gmhisto(bins,hc,nfun,c,alp,mu,var,gmk);// Calculo del gradiente mu
      //Aca debe ir el calculo del hessiano de mu hacerlo 3xnxn 
    //Calculo de los tamanios de paso para mu repetir 3 veces
    //pCauchy(hmk[0],gmkt[0],dk,nfun,pcmk);pNewton(hmk[0],gmkt[0],nfun,pbmk);
    //pDogleg(pcmk,pbmk,dk,nfun,pmk);
    //vector o toda la matriz? 
    // Fin de procesos con mu=========================

    
    //Actualizacion de dk 
    if(rhoak<0.25){
      dk=0.25*dk; 
    }else{
      if(rhoak>0.75 && fabs(Norma_2_vector(pak,nfun)-dk)<1e-6){
        dk=minimo2(2.0*dk,dhat);
      }
      else{
        dk=dk;
      }
    }
    //Criterios de paro
    nga=Norma_2_vector(gak,nfun);
    if(nga<tolg) break;//esto no va :/ ?;
  }

  free(hc); free(auxa); free(auxm);
  free(gak); free(pcak); free(pbak); free(pak);
   free(pcmk); free(pbmk); free(pmk);
  liberar_matriz(hak,nfun);
  liberar_matriz(hmk,nfun); liberar_matriz(auxm,3); liberar_matriz(gmk,nfun);
return(1);}

int pasoMu(double ***hk, double **gk, double dk, int nfun, double **pmk){
  double **gkt=crear_matriz(3,nfun);
  matriz_transponer(gk,3,nfun,gkt);
  double **pkt=crear_matriz(3,nfun);
  double *pc=crear_vector(nfun);
  double *pb=crear_vector(nfun);
  double *pk=crear_vector(nfun);
  pCauchy(hk[0],gkt[0],dk,nfun,pc);pNewton(hk[0],gkt[0],nfun,pb);
  pDogleg(pc,pb,dk,nfun,pk);
  vector_copiar(pk,nfun,pkt[0]);
  pCauchy(hk[1],gkt[1],dk,nfun,pc);pNewton(hk[1],gkt[1],nfun,pb);
  pDogleg(pc,pb,dk,nfun,pk);
  vector_copiar(pk,nfun,pkt[1]);
  pCauchy(hk[2],gkt[2],dk,nfun,pc);pNewton(hk[2],gkt[2],nfun,pb);
  pDogleg(pc,pb,dk,nfun,pk);
  vector_copiar(pk,nfun,pkt[2]);
  matriz_transponer(pkt,nfun,3,pmk);
  free(pc);free(pb);free(pk);
  liberar_matriz(gkt,3);
  liberar_matriz(pkt,3);
return(1);}

