#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include "defs.h"
#include "transformation.h"
#include "assoclegendre.h"
#include "atmprop.h"
#include "misc.h"



void cart2sph(double x, double y, double z, double *R, double *theta, double *phi) { 
  /* Convert cartesian to spherical coordinates */
  double r;
  r = x*x + y*y;
  *R = sqrt(r + z*z);          
  *phi = atan2(z, sqrt(r));   
  *theta = atan2(y,x);                    
}

void sph2cart(double theta, double phi, double R, double *x, double *y, double *z) {
  /* Convert spherical to cartesian coordinates */
  *x = R*cos(phi)*cos(theta);
  *y = R*cos(phi)*sin(theta);
  *z = R*sin(phi);
}

void cart2kep(double rvec[3], double vvec[3], double kep[6]) {

  /* Cartesian position and velocity (km and km/s) to Keplerian orbit elements (angles in radians) */
  /* Makes temporary assumptions for now. Fix later - Mike */

  int j;
  double hvec[3], h, nvec[3];
  double z[3] = {0.0, 0.0, 1.0};
  double r, v;
  double evec[3], e, E;
  double a, p, i;
  double Om, om;
  double u, nu;

  /* INITIALIZE NVEC AND HVEC */
  for(j=0; j<3; j++) {
    hvec[j] = 0.0;
    nvec[j] = 0.0;
  }

  cross(rvec, vvec, hvec);
  h = norm(hvec);
  cross(z, hvec, nvec);

  r = norm(rvec);
  v = norm(vvec);

  evec[0] = ((v*v - JGM3GMe/r)*rvec[0] - dot(rvec, vvec)*vvec[0])/JGM3GMe;
  evec[1] = ((v*v - JGM3GMe/r)*rvec[1] - dot(rvec, vvec)*vvec[1])/JGM3GMe;
  evec[2] = ((v*v - JGM3GMe/r)*rvec[2] - dot(rvec, vvec)*vvec[2])/JGM3GMe;
  
  e = norm(evec);
  E = v*v/2.0 - JGM3GMe/r;

  /* ASSUME NOT HYPERBOLIC */
  a = -JGM3GMe/(2.0*E);
  p = a*(1.0 - e*e);
  i = acos(hvec[2]/h);

  Om = acos(nvec[0]/norm(nvec));
  if (nvec[1] < 0.0) {
    Om = 2*M_PI - Om;
  }

  /* SPECIAL CASE #1 - CIRCULAR ORBIT */
  if (e < 1.0e-16) {
    om = 0.0; 
  } else {
    om = acos(dot(nvec, evec)/(norm(nvec)*e));
    if (evec[2] < 0.0) { 
      om = 2*M_PI - om;
    }
  }

  /* SPECIAL CASE #2 - CIRCULAR INCLINED */
  if (e < 1.0e-16) {
    u = acos(dot(nvec, rvec)/(norm(nvec)*norm(rvec)));
    if (rvec[2] < 0.0) {
      u = 2*M_PI - u;
    }
    nu = u;
  } else {
    nu = acos(dot(evec, rvec)/(e*r));
    if (dot(rvec, vvec) < 0.0) {
      nu = 2*M_PI - nu;
    }
  }

  kep[0] = a;
  kep[1] = e;
  kep[2] = i;
  kep[3] = Om;
  kep[4] = om;
  kep[5] = nu;

}

void kep2cart(double kep[6], double state[6]) {
    
  /* Convert Keplerian elements to cartesian elements */
  /* input: */
  /* kep = keplerian elements */
  /* units assumed in km and rad/s */
  /* output: */
  /* randv has the cartesian position and velocity in km and km/s (in whatever */
  /* frame the orbit elements are measured in, probably J2000). */
    
  double a, e, i, Om, om, nu;
  double p;
  double cnu, snu;
  double cOm, com, ci, sOm, som, si;
  double r_pqw[3], v_pqw[3];
  double RPQW2IJK[3][3];
  double r_ijk[3], v_ijk[3];

  a = kep[0];   /* semi-major axis, km */
  e = kep[1];   /* eccentricity, unitless */
  i = kep[2];   /* inclination, rad */
  Om = kep[3];  /* right-ascension of the ascending node (RAAN), rad */
  om = kep[4];  /* argument of perigee, rad */
  nu = kep[5];  /* true anomaly, rad */

  p = a*(1.0 - e*e);

  cnu = cos(nu);
  snu = sin(nu);

  r_pqw[0] = p*cnu/(1.0 + e*cnu);
  r_pqw[1] = p*snu/(1.0 + e*cnu);
  r_pqw[2] = 0.0;

  v_pqw[0] = -sqrt(JGM3GMe/p)*snu;
  v_pqw[1] = sqrt(JGM3GMe/p)*(e + cnu);
  v_pqw[2] = 0.0;

  cOm = cos(Om);
  com = cos(om);
  ci = cos(i);
  sOm = sin(Om);
  som = sin(om);
  si = sin(i);

  RPQW2IJK[0][0] =  cOm*com - sOm*som*ci;
  RPQW2IJK[0][1] = -cOm*som - sOm*com*ci;
  RPQW2IJK[0][2] =  sOm*si;

  RPQW2IJK[1][0] =  sOm*com + cOm*som*ci;
  RPQW2IJK[1][1] = -sOm*som + cOm*com*ci;
  RPQW2IJK[1][2] = -cOm*si;

  RPQW2IJK[2][0] =  som*si;
  RPQW2IJK[2][1] =  com*si;
  RPQW2IJK[2][2] =  ci;

  matrix_multiply(RPQW2IJK, r_pqw, r_ijk);
  matrix_multiply(RPQW2IJK, v_pqw, v_ijk);
  
  state[0] = r_ijk[0];
  state[1] = r_ijk[1];
  state[2] = r_ijk[2];
  state[3] = v_ijk[0];
  state[4] = v_ijk[1];
  state[5] = v_ijk[2];
  
}

