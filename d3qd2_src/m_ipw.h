/*
 * m_ipw.h
 *
 *  Created on: Aug 21, 2019
 *      Author: ohta
 */

#ifndef M_IPW_H_
#define M_IPW_H_

#include "ipw.h"
#include <string.h>

typedef struct beam_object{
  Ipw *ipw;
}BOBJ;

typedef struct ipw_data{
  char fn_ipw[64];
  int n_ipw;

  double n_0;
  double lambda_0;

  BOBJ bd;
}MIPW;

// ---- m_ipw.c -----------
void read_data_mipw(char *fname,MIPW *obj);
void print_data_mipw(MIPW *obj);
void print_data_mipw_mksa(MIPW *obj);
void check_param_mipw(MIPW *obj);
void setup_mipw(MIPW *obj);
void free_mipw (MIPW *obj);

void calc_mipw_EH(double complex *e,double complex *h,double *x,MIPW *obj);
void calc_mipw_EH_dEHdv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,MIPW *obj);

#endif /* M_IPW_H_ */
