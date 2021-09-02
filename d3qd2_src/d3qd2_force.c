/*
 * d3qd2_force.c
 *
 *  Created on: Oct 1, 2019
 *      Author: ohta
 */

#include "bem3_emf_qd2.h"

int force_FN(double *F,double *N,double *rc,int type,DQD2 *qd)
{
  void bil_force_4p(double *F,double *N,double *rc,int s,DQD2 *qd);
  void bil_force_9p(double *F,double *N,double *rc,int s,DQD2 *qd);
  void lit_force_4p(double *F,double *N,double *rc,int s,DQD2 *qd);
  void lit_force_7p(double *F,double *N,double *rc,int s,DQD2 *qd);

  double tf[3],tn[3],Fx,Fy,Fz,Nx,Ny,Nz;
  int t,td,nc;

  Fx=0.0;    Fy=0.0;    Fz=0.0;
  Nx=0.0;    Ny=0.0;    Nz=0.0;

  nc=0;
  #pragma omp parallel for schedule(dynamic) reduction(+:Fx,Fy,Fz,Nx,Ny,Nz,nc) private(td,tf,tn)
  for(t=1;t<=qd->bd.sb[0].Ne;t++){
    td=qd->bd.sb[0].sid[t];
    nc++;
    if( ELT4==check_element_type(td,&(qd->bd)) ){
      if(type==0) bil_force_4p(tf,tn,rc,t,qd);
      else bil_force_9p(tf,tn,rc,t,qd);
    }
    else {
      if(type==0) lit_force_4p(tf,tn,rc,t,qd);
      else lit_force_7p(tf,tn,rc,t,qd);
    }
    Fx+=tf[0];
    Fy+=tf[1];
    Fz+=tf[2];
    Nx+=tn[0];
    Ny+=tn[1];
    Nz+=tn[2];
  }

  F[0]=Fx;
  F[1]=Fy;
  F[2]=Fz;
  N[0]=Nx;
  N[1]=Ny;
  N[2]=Nz;

  if(nc!=0) return 0;
  else return -1;
}

///////////////////////////////////////////////////////////////////////////////////
void bil_force_4p(double *F,double *N,double *rc,int s,DQD2 *qd)
{
  double complex te[3],th[3],e[3],h[3];
  double epd,aE2,aH2,Tc,T[9],r[3],w[3];
  int i,j,k,asd;

  asd=abs(qd->bd.sb[0].sid[s]);
  epd=creal(qd->n[0])*creal(qd->n[0]);
  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<4;i++){
    for(j=0;j<3;j++){
      e[j]=qd->bd.sb[0].E[s][i][j];
      h[j]=qd->bd.sb[0].H[s][i][j];
      r[j]=qd->bd.ren[asd][i][j];
      w[j]=qd->bd.wen[asd][i][j];
    }
    calc_mipw_EH(te,th,r,&(qd->pw));
    for(j=0;j<3;j++){
      e[j]+=te[j];
      h[j]+=th[j];
    }

    aE2=creal(e[0]*conj(e[0])+e[1]*conj(e[1])+e[2]*conj(e[2]));
    aH2=creal(h[0]*conj(h[0])+h[1]*conj(h[1])+h[2]*conj(h[2]));
    Tc=0.25*(epd*aE2+aH2);

    for(j=0;j<3;j++) for(k=0;k<3;k++) T[3*j+k]=0.5*(epd*creal(e[j]*conj(e[k]))+creal(h[j]*conj(h[k])));
    T[0]-=Tc;
    T[4]-=Tc;
    T[8]-=Tc;
    for(j=0;j<3;j++){
      F[j]+=( T[3*j+0]*w[0]+T[3*j+1]*w[1]+T[3*j+2]*w[2] )*qd->bd.wt_44[i];
      N[j]+=( ((r[1]-rc[1])*T[3*j+2]-(r[2]-rc[2])*T[3*j+1])*w[0]
             +((r[2]-rc[2])*T[3*j+0]-(r[0]-rc[0])*T[3*j+2])*w[1]
             +((r[0]-rc[0])*T[3*j+1]-(r[1]-rc[1])*T[3*j+0])*w[2] )*qd->bd.wt_44[i];
    }
  }

  for(j=0;j<3;j++){
    F[j]*=-1.0;
    N[j]*=-1.0;
  }
}

void bil_force_9p(double *F,double *N,double *rc,int s,DQD2 *qd)
{
  double complex e[3],h[3];
  double epd,aE2,aH2,Tc,T[9],r[3],w[3],cr[3][4],cw[3][3];
  int i,j,k,asd;

  asd=abs(qd->bd.sb[0].sid[s]);
  bil_copy_elem_const_rw(cr,cw,asd,&(qd->bd));

  epd=creal(qd->n[0])*creal(qd->n[0]);
  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<9;i++){
    bil_rw_zeta_eta(r,w,qd->bd.zt_49[i],qd->bd.et_49[i],cr,cw);
    EH_t_bd(e,h,0,s,qd->bd.zt_49[i],qd->bd.et_49[i],PREC_DEF_FN,qd); // integration type

    aE2=creal(e[0]*conj(e[0])+e[1]*conj(e[1])+e[2]*conj(e[2]));
    aH2=creal(h[0]*conj(h[0])+h[1]*conj(h[1])+h[2]*conj(h[2]));
    Tc=0.25*(epd*aE2+aH2);

    for(j=0;j<3;j++) for(k=0;k<3;k++) T[3*j+k]=0.5*(epd*creal(e[j]*conj(e[k]))+creal(h[j]*conj(h[k])));
    T[0]-=Tc;
    T[4]-=Tc;
    T[8]-=Tc;
    for(j=0;j<3;j++){
      F[j]+=( T[3*j+0]*w[0]+T[3*j+1]*w[1]+T[3*j+2]*w[2] )*qd->bd.wt_49[i];
      N[j]+=( ((r[1]-rc[1])*T[3*j+2]-(r[2]-rc[2])*T[3*j+1])*w[0]
             +((r[2]-rc[2])*T[3*j+0]-(r[0]-rc[0])*T[3*j+2])*w[1]
             +((r[0]-rc[0])*T[3*j+1]-(r[1]-rc[1])*T[3*j+0])*w[2] )*qd->bd.wt_49[i];
    }
  }

  for(j=0;j<3;j++){
    F[j]*=-1.0;
    N[j]*=-1.0;
  }
}

void lit_force_4p(double *F,double *N,double *rc,int s,DQD2 *qd)
{
  double complex te[3],th[3],e[3],h[3];
  double epd,aE2,aH2,Tc,T[9],r[3],w[3],cr[3][4],cw[3][3];
  int i,j,k,asd;

  asd=abs(qd->bd.sb[0].sid[s]);
  lit_copy_elem_const_rw(cr,cw,asd,&(qd->bd));

  epd=creal(qd->n[0])*creal(qd->n[0]);
  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<4;i++){
    if(i<3){
      for(j=0;j<3;j++){
        e[j]=qd->bd.sb[0].E[s][i][j];
        h[j]=qd->bd.sb[0].H[s][i][j];
        r[j]=qd->bd.ren[asd][i][j];
        w[j]=qd->bd.wen[asd][i][j];
      }
      calc_mipw_EH(te,th,r,&(qd->pw));
      for(j=0;j<3;j++){
        e[j]+=te[j];
        h[j]+=th[j];
      }
    }
    else {
      lit_rw_zeta_eta(r,w,qd->bd.zt_34[i],qd->bd.et_34[i],cr,cw);
      EH_t_bd(e,h,0,s,qd->bd.zt_34[i],qd->bd.et_34[i],PREC_DEF_FN,qd); //
    }

    aE2=creal(e[0]*conj(e[0])+e[1]*conj(e[1])+e[2]*conj(e[2]));
    aH2=creal(h[0]*conj(h[0])+h[1]*conj(h[1])+h[2]*conj(h[2]));
    Tc=0.25*(epd*aE2+aH2);

    for(j=0;j<3;j++) for(k=0;k<3;k++) T[3*j+k]=0.5*(epd*creal(e[j]*conj(e[k]))+creal(h[j]*conj(h[k])));
    T[0]-=Tc;
    T[4]-=Tc;
    T[8]-=Tc;
    for(j=0;j<3;j++){
      F[j]+=( T[3*j+0]*w[0]+T[3*j+1]*w[1]+T[3*j+2]*w[2] )*qd->bd.wt_34[i];
      N[j]+=( ((r[1]-rc[1])*T[3*j+2]-(r[2]-rc[2])*T[3*j+1])*w[0]
             +((r[2]-rc[2])*T[3*j+0]-(r[0]-rc[0])*T[3*j+2])*w[1]
             +((r[0]-rc[0])*T[3*j+1]-(r[1]-rc[1])*T[3*j+0])*w[2] )*qd->bd.wt_34[i];
    }
  }

  for(j=0;j<3;j++){
    F[j]*=-0.5;
    N[j]*=-0.5;
  }
}

void lit_force_7p(double *F,double *N,double *rc,int s,DQD2 *qd)
{
  double complex e[3],h[3];
  double epd,aE2,aH2,Tc,T[9],r[3],w[3],cr[3][4],cw[3][3];
  int i,j,k,asd;

  asd=abs(qd->bd.sb[0].sid[s]);
  lit_copy_elem_const_rw(cr,cw,asd,&(qd->bd));

  epd=creal(qd->n[0])*creal(qd->n[0]);
  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<7;i++){
    lit_rw_zeta_eta(r,w,qd->bd.zt_37[i],qd->bd.et_37[i],cr,cw);
    EH_t_bd(e,h,0,s,qd->bd.zt_37[i],qd->bd.et_37[i],PREC_DEF_FN,qd); //

    aE2=creal(e[0]*conj(e[0])+e[1]*conj(e[1])+e[2]*conj(e[2]));
    aH2=creal(h[0]*conj(h[0])+h[1]*conj(h[1])+h[2]*conj(h[2]));
    Tc=0.25*(epd*aE2+aH2);

    for(j=0;j<3;j++) for(k=0;k<3;k++) T[3*j+k]=0.5*(epd*creal(e[j]*conj(e[k]))+creal(h[j]*conj(h[k])));
    T[0]-=Tc;
    T[4]-=Tc;
    T[8]-=Tc;
    for(j=0;j<3;j++){
      F[j]+=( T[3*j+0]*w[0]+T[3*j+1]*w[1]+T[3*j+2]*w[2] )*qd->bd.wt_37[i];
      N[j]+=( ((r[1]-rc[1])*T[3*j+2]-(r[2]-rc[2])*T[3*j+1])*w[0]
             +((r[2]-rc[2])*T[3*j+0]-(r[0]-rc[0])*T[3*j+2])*w[1]
             +((r[0]-rc[0])*T[3*j+1]-(r[1]-rc[1])*T[3*j+0])*w[2] )*qd->bd.wt_37[i];
    }
  }

  for(j=0;j<3;j++){
    F[j]*=-0.5;
    N[j]*=-0.5;
  }
}

