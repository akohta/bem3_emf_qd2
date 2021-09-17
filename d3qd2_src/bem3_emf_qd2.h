/*
 * bem3_emf_qd2.h
 *
 *  Created on: Oct 23, 2019
 *      Author: ohta
 */

#ifndef BEM3_EMF_QD2_H_
#define BEM3_EMF_QD2_H_

#include "d3qd2_elem.h"

#define PREC_DEF_BE 2 // type settings for d3qd2_bv_solver
#define PREC_DEF_FN 0 // type settings of field analysis functions in force_FN()

// -- d3qd2_setup.c --
void read_dqd2(int argc,char **argv,DQD2 *dq1);           // read datafile for solver
void print_dqd2(DQD2 *qd);                                // print data 
void print_dqd2_mksa(DQD2 *qd);                           // print data in MKSA system of units
void initialize_dqd2(DQD2 *qd);                           // memory allocation and initialize data for solver
void finalize_dqd2(DQD2 *qd);                             // memory free
int domain_id(double *rt,DQD2 *qd);                       // return domain id of point rt, return the main domain id if on boundary. for non periodic model
int q0_domain_id(double *rt,DQD2 *qd);                    // return domain id of point rt, return tne main domain id if on boundary. for periodic model.
void q0_domain_id_l(int *l1,int *l2,double *rt,DQD2 *qd); // calc domain id l1,l2
void dat_write(char *filename,DQD2 *qd);                  // output analysis result to binary file with specified filename
void dat_read (char *filename,DQD2 *qd);                  // read datafile outputed by dat_write()
void output_node_particles(char *fname,DQD2 *qd);         // outputs the nodes as point cloud data ( .particles file ) 


// -- d3qd2_solve_bieq.c --
void solve_bieq(DQD2 *qd);                                // solve boundary integral equation


// -- d3qd2_field.c --
// analysis modified electromagnetic potential by using boundary integral equations 
int mEMP_s(double complex *U,double *rt,int type,DQD2 *qd); // scattered or internal field
int mEMP_t(double complex *U,double *rt,int type,DQD2 *qd); // total field ( add incident field to scattered field )
int mEMP_i(double complex *U,double *rt,int type,DQD2 *qd); // incident field
// outputs
// U[0]=Ux, U[1]=Uy, U[2]=Uz ( vector potential ), U[3]=phi ( scalar potential ), return object id.
// inputs
// rt[0]=x, rt[1]=y, rt[2]=z, type, pointer of DQD2.
// type=0:4-point GL, type=1:9-pont or 7-point GL, type=2:GLN-point GL, type=3:GHN-point GL, type=4:DE integration.
// the higher the type numbers are, the lower error is, the slower calculation speed is. it is usually set to 0 or 1. 
// the type must be set higher number when analyze near the boundary ( in the order of a element ).

// analysis electromagnetic field by using derivative boundary integral equations. 
int EH_mEMP_s(double complex *E,double complex *H,double *rt,int type,DQD2 *qd); // scattered or internal field
int EH_mEMP_t(double complex *E,double complex *H,double *rt,int type,DQD2 *qd); // total field
int EH_mEMP_i(double complex *E,double complex *H,double *rt,int type,DQD2 *qd); // incidient field
// outputs
// E[0]=Ex, E[1]=Ey, E[2]=Ez, H[0]=Hx, H[1]=Hy, H[2]=Hz, return object id.
// inputs
// rt[0]=x, rt[1]=y, rt[2]=z, type, pointer of DQD2.
// type is the same as mEMP_*()
// these fields are calculated from modified electromagnetic potential using derivative boundary integral equations.
// these are slower than EH_*() functions. the error is smaller than EH_*() functions in far-field.

// analysis electromagnetic field by using boundary integral equations
int EH_s(double complex *E,double complex *H,double *rt,int type,DQD2 *qd); // scattered or internal field
int EH_t(double complex *E,double complex *H,double *rt,int type,DQD2 *qd); // total field
int EH_i(double complex *E,double complex *H,double *rt,int type,DQD2 *qd); // incident field
// outputs and inputs are the same as EH_mEMP_*() functions.

// analysis electromagnetic field on the boundary 
void EH_s_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DQD2 *qd); // scattered or internal field
void EH_t_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DQD2 *qd); // total field
void EH_i_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DQD2 *qd); // incident field
// outputs
// E[0]=Ex, E[1]=Ey, E[2]=Ez, H[0]=Hx, H[1]=Hy, H[2]=Hz.
// inputs
// did:domain id, t:element id defined in each domain, zeta_t,eta_t:point parameter on the element surface ( main domain side coordinate ), type, pointer of DQD2.
// type is the same as EH_*().

// -- d3qd2_force.c --
int force_FN(double *F,double *N,double *rc,int type,DQD2 *md);
// outputs
// F[0]=Fx, F[1]=Fy, F[2]=Fz, N[0]=Nx, N[1]=Ny, N[2]=Nz.
// inputs
// rc[0]=x, rc[1]=y, rc[2]=z ( center of rotation ), type, pointer of DQD2.
// type=0:4-point GL, type!=0:9-point or 7-point GL.


#endif /* BEM3_EMF_QD2_H_ */
