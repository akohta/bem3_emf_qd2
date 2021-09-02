/*
 * d3qd2_field.c
 *
 *  Created on: Sep 20, 2019
 *      Author: ohta
 */

#include "bem3_emf_qd2.h"

int mEMP_s(double complex *U,double *rt,int type,DQD2 *qd)
{
  double complex CC[9],kc,ce;
  double F,r0[3],arg;
  int did,s,sd,i,n,l1,l2;

  if(qd->bd.ps==1){
    q0_domain_id_l(&l1,&l2,rt,qd);
    r0[0]=rt[0]-(double)l1*qd->bd.qd.vd1[0]-(double)l2*qd->bd.qd.vd2[0];
    r0[1]=rt[1]-(double)l1*qd->bd.qd.vd1[1]-(double)l2*qd->bd.qd.vd2[1];
    r0[2]=rt[2];
    did=domain_id(r0,qd);
  }
  else {
    l1=0;
    l2=0;
    r0[0]=rt[0];
    r0[1]=rt[1];
    r0[2]=rt[2];
    did=domain_id(r0,qd);
  }

  for(i=0;i<4;i++) U[i]=0.0;
  F=0.0;
  kc=qd->kn[did];

  for(s=1;s<=qd->bd.sb[did].Ne;s++){
    sd=qd->bd.sb[did].sid[s];

    if(qd->bd.ps==1 && did==0) q0_coef_rt(CC,r0,sd,kc,type,&(qd->bd));
    else coef_rt(CC,r0,sd,kc,type,&(qd->bd));

    for(n=0;n<4;n++)
      for(i=0;i<4;i++) U[i]+=CC[n+0]*qd->bd.sb[did].dU[s][n][i]-CC[n+4]*qd->bd.sb[did].U[s][n][i];

    F+=creal(CC[8]);
  }

  if(did==0) for(i=0;i<4;i++) U[i]/=1.0+F;
  else for(i=0;i<4;i++) U[i]/=F;

  if(l1!=0 || l2!=0){
    arg= qd->bd.qd.vk[0]*((double)l1*qd->bd.qd.vd1[0]+(double)l2*qd->bd.qd.vd2[0])
        +qd->bd.qd.vk[1]*((double)l1*qd->bd.qd.vd1[1]+(double)l2*qd->bd.qd.vd2[1]);
    ce=cos(arg)+I*sin(arg);
    for(i=0;i<4;i++) U[i]*=ce;
  }

  return did;
}

int mEMP_t(double complex *U,double *rt,int type,DQD2 *qd)
{
  double complex e[3],h[3];
  int did,i;

  did=mEMP_s(U,rt,type,qd);
  if(did==0){
    calc_mipw_EH(e,h,rt,&(qd->pw));
    for(i=0;i<3;i++) U[i]+=e[i];
  }

  return did;
}

int mEMP_i(double complex *U,double *rt,int type,DQD2 *qd)
{
  double complex e[3],h[3];
  int did,i;

  if(qd->bd.ps==1) did=q0_domain_id(rt,qd);
  else did=domain_id(rt,qd);
  calc_mipw_EH(e,h,rt,&(qd->pw));

  for(i=0;i<3;i++) U[i]=e[i];
  U[3]=0.0;

  return did;
}

int EH_mEMP_s(double complex *E,double complex *H,double *rt,int type,DQD2 *qd)
{
  double complex CC[9],dC[3][9],kc,ch,ce,U[4],dU[3][4];
  double F,i_C,dF[3],r0[3],arg;
  int did,l1,l2,i,s,sd,p,n,d;

  if(qd->bd.ps==1){
    q0_domain_id_l(&l1,&l2,rt,qd);
    r0[0]=rt[0]-(double)l1*qd->bd.qd.vd1[0]-(double)l2*qd->bd.qd.vd2[0];
    r0[1]=rt[1]-(double)l1*qd->bd.qd.vd1[1]-(double)l2*qd->bd.qd.vd2[1];
    r0[2]=rt[2];
    did=domain_id(r0,qd);
  }
  else {
    l1=0;
    l2=0;
    r0[0]=rt[0];
    r0[1]=rt[1];
    r0[2]=rt[2];
    did=domain_id(r0,qd);
  }

  kc=qd->kn[did];
  ch=qd->pw.lambda_0/(2.0*M_PI*I);

  F=0.0;
  for(i=0;i<3;i++){
    dF[i]=0.0;
    for(p=0;p<4;p++) dU[i][p]=0.0;
  }
  for(i=0;i<4;i++){
    U[i]=0.0;
  }

  for(s=1;s<=qd->bd.sb[did].Ne;s++){
    sd=qd->bd.sb[did].sid[s];

    if(did==0 && qd->bd.ps==1) q0_dcoef_rt_grad(CC,dC,r0,sd,kc,type,&(qd->bd));
    else dcoef_rt_grad(CC,dC,r0,sd,kc,type,&(qd->bd));

    for(d=0;d<3;d++){
      for(n=0;n<4;n++){
        for(p=0;p<4;p++){
          if(d==0) U[p]+=CC[n+0]*qd->bd.sb[did].dU[s][n][p]-CC[n+4]*qd->bd.sb[did].U[s][n][p];
          dU[d][p]+=dC[d][n+0]*qd->bd.sb[did].dU[s][n][p]-dC[d][n+4]*qd->bd.sb[did].U[s][n][p];
        }
      }
      dF[d]+=creal(dC[d][8]);
    }
    F+=creal(CC[8]);
  }

  if(fabs(dF[0])>CBD_CDF || fabs(dF[1])>CBD_CDF || fabs(dF[2])>CBD_CDF){
    for(i=0;i<3;i++){
      E[i]=0.0;
      H[i]=0.0;
    }
    if(did==0) return -OPENDID;
    else return -did;
  }

  if(did==0) i_C=1.0/(1.0+F);
  else i_C=1.0/F;
  for(p=0;p<4;p++) U[p]*=i_C;
  for(p=0;p<3;p++)
    for(i=0;i<4;i++) dU[p][i]=(dU[p][i]-U[i]*dF[p])*i_C;

  for(p=0;p<3;p++) E[p]=-dU[p][3]+U[p];
  H[0]=ch*(dU[1][2]-dU[2][1]);
  H[1]=ch*(dU[2][0]-dU[0][2]);
  H[2]=ch*(dU[0][1]-dU[1][0]);

  if(l1!=0 || l2!=0){
    arg= qd->bd.qd.vk[0]*((double)l1*qd->bd.qd.vd1[0]+(double)l2*qd->bd.qd.vd2[0])
        +qd->bd.qd.vk[1]*((double)l1*qd->bd.qd.vd1[1]+(double)l2*qd->bd.qd.vd2[1]);
    ce=cos(arg)+I*sin(arg);
    for(i=0;i<3;i++){
      E[i]*=ce;
      H[i]*=ce;
    }
  }

  return did;
}

int EH_mEMP_t(double complex *E,double complex *H,double *rt,int type,DQD2 *qd)
{
  double complex e[3],h[3];
  int did,i;

  did=EH_mEMP_s(E,H,rt,type,qd);
  if(did==0){
    calc_mipw_EH(e,h,rt,&(qd->pw));
    for(i=0;i<3;i++){
      E[i]+=e[i];
      H[i]+=h[i];
    }
  }

  return did;
}

int EH_mEMP_i(double complex *E,double complex *H,double *rt,int type,DQD2 *qd)
{
  int did;

  if(qd->bd.ps==1) did=q0_domain_id(rt,qd);
  else did=domain_id(rt,qd);
  calc_mipw_EH(E,H,rt,&(qd->pw));

  return did;
}

int EH_s(double complex *E,double complex *H,double *rt,int type,DQD2 *qd)
{
  double complex CC[9],ce;
  double F,r0[3],arg;
  int did,s,sd,i,n,l1,l2;

  if(qd->bd.ps==1){
    q0_domain_id_l(&l1,&l2,rt,qd);
    r0[0]=rt[0]-(double)l1*qd->bd.qd.vd1[0]-(double)l2*qd->bd.qd.vd2[0];
    r0[1]=rt[1]-(double)l1*qd->bd.qd.vd1[1]-(double)l2*qd->bd.qd.vd2[1];
    r0[2]=rt[2];
    did=domain_id(r0,qd);
  }
  else {
    l1=0;
    l2=0;
    r0[0]=rt[0];
    r0[1]=rt[1];
    r0[2]=rt[2];
    did=domain_id(r0,qd);
  }

  for(i=0;i<3;i++){
    E[i]=0.0;
    H[i]=0.0;
  }
  F=0.0;
  for(s=1;s<=qd->bd.sb[did].Ne;s++){
    sd=qd->bd.sb[did].sid[s];

    if(did==0 && qd->bd.ps==1) q0_coef_rt(CC,r0,sd,qd->kn[did],type,&(qd->bd));
    else coef_rt(CC,r0,sd,qd->kn[did],type,&(qd->bd));

    for(n=0;n<4;n++){
      for(i=0;i<3;i++){
        E[i]+=CC[n+0]*qd->bd.sb[did].dE[s][n][i]-CC[n+4]*qd->bd.sb[did].E[s][n][i];
        H[i]+=CC[n+0]*qd->bd.sb[did].dH[s][n][i]-CC[n+4]*qd->bd.sb[did].H[s][n][i];
      }
    }
    F+=creal(CC[8]);
  }

  if(did==0){
    for(i=0;i<3;i++){
      E[i]/=1.0+F;
      H[i]/=1.0+F;
    }
  }
  else{
    for(i=0;i<3;i++){
      E[i]/=F;
      H[i]/=F;
    }
  }

  if(l1!=0 || l2!=0){
    arg= qd->bd.qd.vk[0]*((double)l1*qd->bd.qd.vd1[0]+(double)l2*qd->bd.qd.vd2[0])
        +qd->bd.qd.vk[1]*((double)l1*qd->bd.qd.vd1[1]+(double)l2*qd->bd.qd.vd2[1]);
    ce=cos(arg)+I*sin(arg);
    for(i=0;i<3;i++){
      E[i]*=ce;
      H[i]*=ce;
    }
  }

  return did;
}

int EH_t(double complex *E,double complex *H,double *rt,int type,DQD2 *qd)
{
  double complex e[3],h[3];
  int did,i;

  did=EH_s(E,H,rt,type,qd);
  if(did==0){
    calc_mipw_EH(e,h,rt,&(qd->pw));
    for(i=0;i<3;i++){
      E[i]+=e[i];
      H[i]+=h[i];
    }
  }

  return did;
}

int EH_i(double complex *E,double complex *H,double *rt,int type,DQD2 *qd)
{
  int did;

  if(qd->bd.ps==1) did=q0_domain_id(rt,qd);
  else did=domain_id(rt,qd);
  calc_mipw_EH(E,H,rt,&(qd->pw));

  return did;
}

void EH_s_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DQD2 *qd)
{
  double complex CC[9];
  double F,rt[3];
  int s,sd,l,n,td;

  for(l=0;l<3;l++){
    E[l]=0.0;
    H[l]=0.0;
  }
  F=0.0;

  td=qd->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(qd->bd));

  for(s=1;s<=qd->bd.sb[did].Ne;s++){
    sd=qd->bd.sb[did].sid[s];

    if(did==0 && qd->bd.ps==1) q0_coef_bd(CC,rt,td,zeta_t,eta_t,sd,qd->kn[did],type,&(qd->bd));
    else coef_bd(CC,rt,td,zeta_t,eta_t,sd,qd->kn[did],type,&(qd->bd));

    for(n=0;n<4;n++){
      for(l=0;l<3;l++){
        E[l]+=CC[n+0]*qd->bd.sb[did].dE[s][n][l]-CC[n+4]*qd->bd.sb[did].E[s][n][l];
        H[l]+=CC[n+0]*qd->bd.sb[did].dH[s][n][l]-CC[n+4]*qd->bd.sb[did].H[s][n][l];
      }
    }
    F+=creal(CC[8]);
  }

  if(did==0){
    for(l=0;l<3;l++){
      E[l]/=1.0+F;
      H[l]/=1.0+F;
    }
  }
  else{
    for(l=0;l<3;l++){
      E[l]/=F;
      H[l]/=F;
    }
  }
}

void EH_t_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DQD2 *qd)
{
  double complex CC[9],e[3],h[3];
  double F,rt[3];
  int s,sd,l,n,td;

  for(l=0;l<3;l++){
    E[l]=0.0;
    H[l]=0.0;
  }
  F=0.0;

  td=qd->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(qd->bd));

  for(s=1;s<=qd->bd.sb[did].Ne;s++){
    sd=qd->bd.sb[did].sid[s];

    if(did==0 && qd->bd.ps==1) q0_coef_bd(CC,rt,td,zeta_t,eta_t,sd,qd->kn[did],type,&(qd->bd));
    else coef_bd(CC,rt,td,zeta_t,eta_t,sd,qd->kn[did],type,&(qd->bd));

    for(n=0;n<4;n++){
      for(l=0;l<3;l++){
        E[l]+=CC[n+0]*qd->bd.sb[did].dE[s][n][l]-CC[n+4]*qd->bd.sb[did].E[s][n][l];
        H[l]+=CC[n+0]*qd->bd.sb[did].dH[s][n][l]-CC[n+4]*qd->bd.sb[did].H[s][n][l];
      }
    }
    F+=creal(CC[8]);
  }

  if(did==0){
    calc_mipw_EH(e,h,rt,&(qd->pw));
    for(l=0;l<3;l++){
      E[l]=E[l]/(1.0+F)+e[l];
      H[l]=H[l]/(1.0+F)+h[l];
    }
  }
  else{
    for(l=0;l<3;l++){
      E[l]/=F;
      H[l]/=F;
    }
  }
}

void EH_i_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DQD2 *qd)
{
  double rt[3];
  int td,i;

  td=qd->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(qd->bd));

  if(did==0) calc_mipw_EH(E,H,rt,&(qd->pw));
  else {
    for(i=0;i<3;i++){
      E[i]=0.0;
      H[i]=0.0;
    }
  }
}
