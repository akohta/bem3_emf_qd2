/*
 * bem3_emf_qd2.h
 *
 *  Created on: Oct 23, 2019
 *      Author: ohta
 */

#ifndef BEM3_EMF_QD2_H_
#define BEM3_EMF_QD2_H_

#include "d3qd2_elem.h"

#define PREC_DEF_BE 2 // type settings for d3qd2_bv_solver.
#define PREC_DEF_FN 0 // type settings of field analysis functions in force_FN().

// -- d3qd2_setup.c --
void read_dqd2(int argc,char **argv,DQD2 *dq1);
void print_dqd2(DQD2 *qd);
void print_dqd2_mksa(DQD2 *qd);
void initialize_dqd2(DQD2 *qd);
void finalize_dqd2(DQD2 *qd);
int domain_id(double *rt,DQD2 *qd); // ret basic domain id of point rt, return the main domain id if on boundary, non periodicity
int q0_domain_id(double *rt,DQD2 *qd);
void q0_domain_id_l(int *l1,int *l2,double *rt,DQD2 *qd); // calc domain id l1,l2
void dat_write(char *filename,DQD2 *qd);
void dat_read (char *filename,DQD2 *qd);


// -- d3qd2_solve_bieq.c --
void solve_bieq(DQD2 *qd);


// -- d3qd2_field.c --
// U[0~3]:U_x~U_z,U[4]=phi. type=0:4pGL,1:9pGL,2:GLNpointGL,3:GHNpoint GL,4:DE. return domain id.
int mEMP_s(double complex *U,double *rt,int type,DQD2 *qd); // scattered or internal field
int mEMP_t(double complex *U,double *rt,int type,DQD2 *qd); // total field
int mEMP_i(double complex *U,double *rt,int type,DQD2 *qd); // incident field

int EH_mEMP_s(double complex *E,double complex *H,double *rt,int type,DQD2 *qd); // using modified electromagnetic potential boundary value
int EH_mEMP_t(double complex *E,double complex *H,double *rt,int type,DQD2 *qd);
int EH_mEMP_i(double complex *E,double complex *H,double *rt,int type,DQD2 *qd);

int EH_s(double complex *E,double complex *H,double *rt,int type,DQD2 *qd); // using electromagnetic field boundary value
int EH_t(double complex *E,double complex *H,double *rt,int type,DQD2 *qd);
int EH_i(double complex *E,double complex *H,double *rt,int type,DQD2 *qd);

void EH_s_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DQD2 *qd);
void EH_t_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DQD2 *qd);
void EH_i_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DQD2 *qd);


// -- d3qd2_force.c --
int force_FN(double *F,double *N,double *rc,int type,DQD2 *md);

#endif /* BEM3_EMF_QD2_H_ */
