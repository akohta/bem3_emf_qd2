/*
 * m_ipw.c
 *
 *  Created on: Aug 21, 2019
 *      Author: ohta
 */
#include "m_ipw.h"

void read_data_mipw(char *fname,MIPW *obj)
{
  FILE *fp;
  char buf[256]="";
  int i,nn;
  double tmpd,tmpd2;

  strcpy(obj->fn_ipw,fname);

  if((fp=fopen(fname,"rt"))==NULL){
    printf("read_data_mipw(), failed to read %s. Exit...\n",fname);
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("m_ipw.c, read_data_mipw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("m_ipw.c, read_data_mipw(), failed to read the line. exit...\n");
    exit(1);
  }

  if(fscanf(fp,"%d\n",&nn)!=1){
    printf("m_ipw.c, read_data_mipw(), failed to read the nn. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("m_ipw.c, read_data_mipw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(nn==0) {
    printf("read_data_mpiw(), no wave defined. Exit...\n");
    exit(1);
  }
  obj->n_ipw=nn;     obj->bd.ipw=(Ipw *)malloc(sizeof(Ipw)*nn);

  for(i=0;i<nn;i++){
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("m_ipw.c, read_data_mipw(), failed to read the lambda0. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].lambda0=tmpd;
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("m_ipw.c, read_data_mipw(), failed to read the ni. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].ni     =tmpd;
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("m_ipw.c, read_data_mipw(), failed to read the power. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].power  =tmpd;
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("m_ipw.c, read_data_mipw(), failed to read the real(e0x). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("m_ipw.c, read_data_mipw(), failed to read the imag(e0x). exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].e0x    =tmpd+I*tmpd2;
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("m_ipw.c, read_data_mipw(), failed to read the real(e0y). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("m_ipw.c, read_data_mipw(), failed to read the imag(e0y). exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].e0y    =tmpd+I*tmpd2;
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("m_ipw.c, read_data_mipw(), failed to read the fx. exit...\n");
      exit(1);  
    }
    obj->bd.ipw[i].fx     =tmpd;
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("m_ipw.c, read_data_mipw(), failed to read the fy. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].fy     =tmpd;
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("m_ipw.c, read_data_mipw(), failed to read the fz. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].fz     =tmpd;
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("m_ipw.c, read_data_mipw(), failed to read the theta. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].theta  =tmpd;
    if(fscanf(fp,"%lf\n",&tmpd)!=1){
      printf("m_ipw.c, read_data_mipw(), failed to read the phi. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].phi    =tmpd;
  }
  fclose(fp);
}

void print_data_mipw(MIPW *obj)
{
  int i;

  for(i=0;i<obj->n_ipw;i++){
    printf(" \"%s\" No.%02d ",obj->fn_ipw,i);    print_data_ipw(&(obj->bd.ipw[i]));
  }
}

void print_data_mipw_mksa(MIPW *obj)
{
  int i;

  for(i=0;i<obj->n_ipw;i++){
    printf(" \"%s\" No.%02d ",obj->fn_ipw,i);    print_data_ipw_mksa(&(obj->bd.ipw[i]));
  }
}

void check_param_mipw(MIPW *obj)
{
  int i;
  int st=0;

  obj->n_0     =0.0;
  obj->lambda_0=0.0;

  for(i=0;i<obj->n_ipw;i++){
    if(st==0){
      obj->n_0     =obj->bd.ipw[i].ni;
      obj->lambda_0=obj->bd.ipw[i].lambda0;
      st++;
    }
    else {
      if(obj->n_0     !=obj->bd.ipw[i].ni     ){ printf("check_param_mipw(),parameter error ipw ni, exit\n");     exit(1);}
      if(obj->lambda_0!=obj->bd.ipw[i].lambda0){ printf("check_param_mipw(),parameter error ipw lambda, exit\n"); exit(1);}
    }
  }
  
  obj->omega=2.0*M_PI/obj->lambda_0;
}

void setup_mipw(MIPW *obj)
{
  int i;

  for(i=0;i<obj->n_ipw;i++)      setup_ipw(&(obj->bd.ipw[i]));
  
}

void free_mipw(MIPW *obj)
{

  free(obj->bd.ipw);  obj->n_ipw=0;
  strcpy(obj->fn_ipw,"");
}

void calc_mipw_EH(double complex *e,double complex *h,double *x,MIPW *obj)
{
  double complex te[3],th[3];
  int i,j;
  for(i=0;i<3;i++){
    e[i]=0.0;    h[i]=0.0;
  }

  // ipw
  for(i=0;i<obj->n_ipw;i++){
    calc_ipw_EH(te,th,x,&(obj->bd.ipw[i]));
    for(j=0;j<3;j++){
      e[j]+=te[j];      h[j]+=th[j];
    }
  }
}

void calc_mipw_EH_dEHdv(double complex *e,double complex *h,double complex *dedn,double complex *dhdn,
            double *x,double *n,MIPW *obj)
{
  double complex te[3], th[3],tde[3],tdh[3];
  int i, j;
  for (i = 0; i < 3; i++) {
    e[i] = 0.0;    h[i] = 0.0;
    dedn[i]=0.0;    dhdn[i]=0.0;
  }

  for (i = 0; i < obj->n_ipw; i++) {
    calc_ipw_dEHdv(te, th, tde,tdh,x,n, &(obj->bd.ipw[i]));
    for (j = 0; j < 3; j++) {
      e[j] += te[j];      h[j] += th[j];
      dedn[j]+=tde[j];      dhdn[j]+=tdh[j];
    }
  }
}
