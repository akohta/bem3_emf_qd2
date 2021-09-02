/*
 * ipw.h
 *
 *  Created on: Dec 10, 2018
 *      Author: ohta
 */

#ifndef IPW_H_
#define IPW_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "osu_mksa.h"

typedef struct incident_planewave_data{
  double E0;                // amplitude
  double ki;                // wave number
  double cos_t,sin_t;       // cos(theta),sin(theta), theta and phi are paramter for rotation.
  double cos_p,sin_p;       // cos(phi),sin(phi)
  double complex ex,ey,ez;  // electric field vector
  double complex hx,hy,hz;  // magnetic field vector
}IpwD;
typedef struct incident_planewave{
  double lambda0;          // wavelength of plane wave in vacuum [m]
  double ni;               // refractive index of surrounding
  double power;            // power of planewave per 1 m^2 [W/m^2]
  double complex e0x;      // polarization index of x-direction include phase
  double complex e0y;      // polarization index of y-direction include phase
  double fx,fy,fz;         // translation vector
  double theta;            // theta and phi are parameter for rotation
  double phi;              // unit wave vector k : kx=sin(theta)*cos(phi),ky=sin(theta)*sin(phi),kz=cos(theta)
  IpwD data;               // data
}Ipw;

void read_data_ipw(char *rfile,Ipw *ipw);
void print_data_ipw(Ipw *ipw);
void print_data_ipw_mksa(Ipw *ipw);
void setup_ipw(Ipw *ipw);

void calc_ipw_E (double complex *e,double *x,Ipw *ipw);
void calc_ipw_dEdv(double complex *e,double complex *dedn,double *x,double *v,Ipw *ipw);
void calc_ipw_H (double complex *h,double *x,Ipw *ipw);
void calc_ipw_dHdn(double complex *h,double complex *dhdn,double *x,double *n,Ipw *ipw);
void calc_ipw_EH(double complex *e,double complex *h,double *x,Ipw *ipw);
void calc_ipw_dEHdv(double complex *e,double complex *h,
        double complex *dedn,double complex *dhdn,double *x,double *v,Ipw *ipw);

#endif /* IPW_H_ */
