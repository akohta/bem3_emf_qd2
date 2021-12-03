/*
 * ipw.c
 *
 *  Created on: Dec 10, 2018
 *      Author: ohta
 */
#include "ipw.h"

void read_data_ipw(char *rfile,Ipw *ipw)
{
  FILE *fp;
  if((fp=fopen(rfile,"rt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  char buf[256]="";  int tmpi;  double tmpd,tmpd2;
  if(fgets(buf,256,fp)==NULL){
    printf("ipw.c, read_data_ipw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("ipw.c, read_data_ipw(), failed to read the line. exit...\n");
    exit(1);
  }

  printf("-- wave parameter --\n");
  if(fscanf(fp,"%d",&tmpi)!=1){
    printf("ipw.c, read_data_ipw(), failed to read the fp. exit...\n");
    exit(1);
  }
  if(tmpi!=0){ printf("beam type %d is not supported\n",tmpi); exit(1);}
  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("ipw.c, read_data_ipw(), failed to read the lambda0. exit...\n");
    exit(1);
  }
  ipw->lambda0=tmpd; 
  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("ipw.c, read_data_ipw(), failed to read the ni. exit...\n");
    exit(1);
  }
  ipw->ni     =tmpd; 
  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("ipw.c, read_data_ipw(), failed to read the power. exit...\n");
    exit(1);
  }
  ipw->power  =tmpd;
  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("ipw.c, read_data_ipw(), failed to read the real(e0x). exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&tmpd2)!=1){
    printf("ipw.c, read_data_ipw(), failed to read the imag(e0x). exit...\n");
    exit(1);
  }
  ipw->e0x    =tmpd+I*tmpd2; 
  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("ipw.c, read_data_ipw(), failed to read the real(e0y). exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&tmpd2)!=1){
    printf("ipw.c, read_data_ipw(), failed to read the imag(e0y). exit...\n");
    exit(1);
  }
  ipw->e0y    =tmpd+I*tmpd2;
  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("ipw.c, read_data_ipw(), failed to read the fx. exit...\n");
    exit(1);
  }
  ipw->fx     =tmpd; 
  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("ipw.c, read_data_ipw(), failed to read the fy. exit...\n");
    exit(1);
  }
  ipw->fy     =tmpd;
  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("ipw.c, read_data_ipw(), failed to read the fz. exit...\n");
    exit(1);
  }
  ipw->fz     =tmpd;
  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("ipw.c, read_data_ipw(), failed to read the theta. exit...\n");
    exit(1);
  }
  ipw->theta  =tmpd; 
  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("ipw.c, read_data_ipw(), failed to read the phi. exit...\n");
    exit(1);
  }
  ipw->phi    =tmpd;

  fclose(fp);
}

void print_data_ipw(Ipw *ipw)
{
  printf("-- plane wave parameter --\n");
  printf("wave length of incident beam in vacuum      : %15.14g\n",ipw->lambda0);
  printf("refractive index of surrounding             : %15.14g\n",ipw->ni);
  printf("incident beam power per unit area           : %15.14g\n",ipw->power);
  printf("x-poralization coefficient                  :%7.6g+%7.6gi\n",creal(ipw->e0x),cimag(ipw->e0x));
  printf("y-poralization coefficient                  :%7.6g+%7.6gi\n",creal(ipw->e0y),cimag(ipw->e0y));
  printf("x-component of translation vector           : %15.14g\n",ipw->fx);
  printf("y-component of translation vector           : %15.14g\n",ipw->fy);
  printf("z-component of translation vector           : %15.14g\n",ipw->fz);
  printf("rotation parameter theta               [rad]: %15.14g\n",ipw->theta);
  printf("rotation parameter phi                 [rad]: %15.14g\n",ipw->phi);
}

void print_data_ipw_mksa(Ipw *ipw)
{
  printf("-- plane wave parameter --\n");
  printf("wave length of incident beam in vacuum   [m]: %15.14g\n",OSUtoMKSA_length(ipw->lambda0));
  printf("refractive index of surrounding             : %15.14g\n",ipw->ni);
  printf("incident beam power per unit area    [W/m^2]: %15.14g\n",OSUtoMKSA_Power_per_unit_area(ipw->power));
  printf("x-poralization coefficient                  :%7.6g+%7.6gi\n",creal(ipw->e0x), cimag(ipw->e0x));
  printf("y-poralization coefficient                  :%7.6g+%7.6gi\n",creal(ipw->e0y), cimag(ipw->e0y));
  printf("x-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(ipw->fx));
  printf("y-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(ipw->fy));
  printf("z-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(ipw->fz));
  printf("rotation parameter theta               [rad]: %15.14g\n",ipw->theta);
  printf("rotation parameter phi                 [rad]: %15.14g\n",ipw->phi);
}

void setup_ipw(Ipw *ipw)
{
  double cs=1.0/sqrt(pow(cabs(ipw->e0x),2)+pow(cabs(ipw->e0y),2));
  double complex ex,ey,hx,hy;
  double sin_t,cos_t,sin_p,cos_p;
  ex= ipw->e0x*cs;
  ey= ipw->e0y*cs;
  hx=-ey*ipw->ni;
  hy= ex*ipw->ni;

  ipw->data.E0=sqrt(ipw->power*2.0/ipw->ni);
  ipw->data.ki=2.0*M_PI*ipw->ni/ipw->lambda0;

  sin_t=sin(ipw->theta);  cos_t=cos(ipw->theta);
  sin_p=sin(ipw->phi  );  cos_p=cos(ipw->phi  );
  ipw->data.sin_t=sin_t;  ipw->data.cos_t=cos_t;
  ipw->data.sin_p=sin_p;  ipw->data.cos_p=cos_p;

  ipw->data.ex=ex*(sin_p*sin_p+cos_p*cos_p*cos_t)+ey*(sin_p*cos_p*cos_t-sin_p*cos_p);
  ipw->data.ey=ex*(sin_p*cos_p*cos_t-sin_p*cos_p)+ey*(sin_p*sin_p*cos_t+cos_p*cos_p);
  ipw->data.ez=-(ex*cos_p+ey*sin_p)*sin_t;
  ipw->data.hx=hx*(sin_p*sin_p+cos_p*cos_p*cos_t)+hy*(sin_p*cos_p*cos_t-sin_p*cos_p);
  ipw->data.hy=hx*(sin_p*cos_p*cos_t-sin_p*cos_p)+hy*(sin_p*sin_p*cos_t+cos_p*cos_p);
  ipw->data.hz=-(hx*cos_p+hy*sin_p)*sin_t;
}

void calc_ipw_E (double complex *e,double *x,Ipw *ipw)
{
  double xc=x[0]-ipw->fx;
  double yc=x[1]-ipw->fy;
  double zc=x[2]-ipw->fz;
  double te=ipw->data.ki*(xc*ipw->data.sin_t*ipw->data.cos_p+yc*ipw->data.sin_t*ipw->data.sin_p+zc*ipw->data.cos_t);
  double complex Ee=(cos(te)+I*sin(te));

  e[0]=ipw->data.ex*Ee*ipw->data.E0;
  e[1]=ipw->data.ey*Ee*ipw->data.E0;
  e[2]=ipw->data.ez*Ee*ipw->data.E0;
}

void calc_ipw_dEdv(double complex *e, double complex *dedn, double *x,double *n, Ipw *ipw)
{
  double xc = x[0] - ipw->fx;
  double yc = x[1] - ipw->fy;
  double zc = x[2] - ipw->fz;
  double te = ipw->data.ki*(  xc * ipw->data.sin_t * ipw->data.cos_p
                            + yc * ipw->data.sin_t * ipw->data.sin_p
                            + zc * ipw->data.cos_t);
  double complex Ee = (cos(te) + I * sin(te));
  double complex gcf = I * ipw->data.ki*(  ipw->data.sin_t * ipw->data.cos_p * n[0]
                                         + ipw->data.sin_t * ipw->data.sin_p * n[1]
                                         + ipw->data.cos_t * n[2]);

  e[0] = ipw->data.ex * Ee * ipw->data.E0;
  e[1] = ipw->data.ey * Ee * ipw->data.E0;
  e[2] = ipw->data.ez * Ee * ipw->data.E0;

  dedn[0] = ipw->data.ex * Ee * ipw->data.E0 * gcf;
  dedn[1] = ipw->data.ey * Ee * ipw->data.E0 * gcf;
  dedn[2] = ipw->data.ez * Ee * ipw->data.E0 * gcf;
}

void calc_ipw_H (double complex *h,double *x,Ipw *ipw)
{
  double xc=x[0]-ipw->fx;
  double yc=x[1]-ipw->fy;
  double zc=x[2]-ipw->fz;
  double te=ipw->data.ki*(xc*ipw->data.sin_t*ipw->data.cos_p+yc*ipw->data.sin_t*ipw->data.sin_p+zc*ipw->data.cos_t);
  double complex Ee=(cos(te)+I*sin(te));

  h[0]=ipw->data.hx*Ee*ipw->data.E0;
  h[1]=ipw->data.hy*Ee*ipw->data.E0;
  h[2]=ipw->data.hz*Ee*ipw->data.E0;
}

void calc_ipw_dHdv(double complex *h,double complex *dhdn,double *x,double *n,Ipw *ipw)
{
  double xc = x[0] - ipw->fx;
  double yc = x[1] - ipw->fy;
  double zc = x[2] - ipw->fz;
  double te = ipw->data.ki*(  xc * ipw->data.sin_t * ipw->data.cos_p
                            + yc * ipw->data.sin_t * ipw->data.sin_p
                            + zc * ipw->data.cos_t);
  double complex Ee = (cos(te) + I * sin(te));
  double complex gcf = I * ipw->data.ki*(  ipw->data.sin_t * ipw->data.cos_p * n[0]
                                         + ipw->data.sin_t * ipw->data.sin_p * n[1]
                                         + ipw->data.cos_t * n[2]);

  h[0] = ipw->data.hx * Ee * ipw->data.E0;
  h[1] = ipw->data.hy * Ee * ipw->data.E0;
  h[2] = ipw->data.hz * Ee * ipw->data.E0;

  dhdn[0] = ipw->data.hx * Ee * ipw->data.E0 * gcf;
  dhdn[1] = ipw->data.hy * Ee * ipw->data.E0 * gcf;
  dhdn[2] = ipw->data.hz * Ee * ipw->data.E0 * gcf;
}

void calc_ipw_EH(double complex *e,double complex *h,double *x,Ipw *ipw)
{
  double xc=x[0]-ipw->fx;
  double yc=x[1]-ipw->fy;
  double zc=x[2]-ipw->fz;
  double te=ipw->data.ki*(xc*ipw->data.sin_t*ipw->data.cos_p+yc*ipw->data.sin_t*ipw->data.sin_p+zc*ipw->data.cos_t);
  double complex Ee=(cos(te)+I*sin(te));

  e[0]=ipw->data.ex*Ee*ipw->data.E0;  e[1]=ipw->data.ey*Ee*ipw->data.E0;  e[2]=ipw->data.ez*Ee*ipw->data.E0;
  h[0]=ipw->data.hx*Ee*ipw->data.E0;  h[1]=ipw->data.hy*Ee*ipw->data.E0;  h[2]=ipw->data.hz*Ee*ipw->data.E0;
}

void calc_ipw_dEHdv(double complex *e,double complex *h,
        double complex *dedn,double complex *dhdn,double *x,double *n,Ipw *ipw)
{
  double xc = x[0] - ipw->fx;
  double yc = x[1] - ipw->fy;
  double zc = x[2] - ipw->fz;
  double te = ipw->data.ki*(  xc * ipw->data.sin_t * ipw->data.cos_p
                            + yc * ipw->data.sin_t * ipw->data.sin_p
                            + zc * ipw->data.cos_t);
  double complex Ee = (cos(te) + I * sin(te));
  double complex gcf = I * ipw->data.ki*(  ipw->data.sin_t * ipw->data.cos_p * n[0]
                                         + ipw->data.sin_t * ipw->data.sin_p * n[1]
                                         + ipw->data.cos_t * n[2]);

  e[0] = ipw->data.ex * Ee * ipw->data.E0;
  e[1] = ipw->data.ey * Ee * ipw->data.E0;
  e[2] = ipw->data.ez * Ee * ipw->data.E0;
  h[0] = ipw->data.hx * Ee * ipw->data.E0;
  h[1] = ipw->data.hy * Ee * ipw->data.E0;
  h[2] = ipw->data.hz * Ee * ipw->data.E0;

  dedn[0] = ipw->data.ex * Ee * ipw->data.E0 * gcf;
  dedn[1] = ipw->data.ey * Ee * ipw->data.E0 * gcf;
  dedn[2] = ipw->data.ez * Ee * ipw->data.E0 * gcf;
  dhdn[0] = ipw->data.hx * Ee * ipw->data.E0 * gcf;
  dhdn[1] = ipw->data.hy * Ee * ipw->data.E0 * gcf;
  dhdn[2] = ipw->data.hz * Ee * ipw->data.E0 * gcf;
}

