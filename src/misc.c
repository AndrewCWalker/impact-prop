#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "defs.h"
#include "misc.h"


/************************** HEIGHT ABOVE ELLIPSOID  **********************************/
double hellipsoid(double rvec[]) {

  /* function to calculate the altitude (km) above an oblate Earth, */
  /* with equatorial radius and polar radius. rvec is the ITRF */
  /* position vector in km. */
  
  /* NOTE: past versions of this function may have used rvec as J2000 instead */
  /* of ITRF. */



  double r, z;
  double Re, Rp;
  double hellp;

  r = norm(rvec);
  z = rvec[2];

  Re = JGM3Re;
  Rp = JGM3Rp;

  hellp = r - Re/(sqrt(1.0 + (z*z/(r*r))*(Re*Re/(Rp*Rp) - 1.0)));

  return(hellp);

}

/************************** FIND NEAREST NEIGHBOR  **********************************/
int findnearest(double et, double dt, double et_doubletime[], int dnp) {

  /* Search through array and find the nearest point */
  int itime = 0;
  double diff1, diff2;

  itime = (et - et_doubletime[0])/(dt/2.0);
  if(itime < 0) {
    itime=0;
  } 

  if (itime >= dnp) {
    itime = dnp - 1;
  } else {
    diff1 = fabs(et-et_doubletime[itime]);
    diff2 = fabs(et-et_doubletime[itime+1]);

    if(diff2 < diff1)
      itime++;
    
    while (et_doubletime[itime]==-1)
      itime--;
  }

  return(itime);

}



/*********************************** DOT PRODUCT *************************************/
double dot(double V[], double W[]) {
   
  /* Computes dot product between vectors V and W */
  /* Inputs: V = Vector #1 [Vx, Vy, Vz] */
  /*         W = Vector #2 [Wx, Wy, Wz] */
  /* Output: a = Dot Product of V and W */

  double a = 0.0;

  a = V[0]*W[0] + V[1]*W[1] + V[2]*W[2];
  return(a);

}

/*********************************** DOT PRODUCT *************************************/
void cross(double V[], double W[], double VEC[]) {
  /* Computes cross product between vectors V and W */

  /* Inputs: V = Vector #1 [Vx, Vy, Vz] */
  /*         W = Vector #2 [Wx, Wy, Wz] */

  /* Output: b = Cross Product of V and W */

  VEC[0] = V[1]*W[2] - V[2]*W[1];
  VEC[1] = V[2]*W[0] - V[0]*W[2];
  VEC[2] = V[0]*W[1] - V[1]*W[0];

}


/*********************************** VECTOR NORM  *************************************/
double norm(double V[]) {

  /* Computes cross product between vectors V and W */
  /* Inputs: V = Vector #1 [Vx, Vy, Vz] */
  /* Output: n = Vector Norm of V*/

  double n;

  n = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);

  return(n);

}


/*********************************** MATRIX MULTIPLY *************************************/
void matrix_multiply(double A[3][3], double V[3], double VEC[3]) {
  /* Computes cross product between vectors V and W */

  /* Inputs: A = 3x3 Matrix */
  /*         V = Vector [Vx, Vy, Vz] */

  /* Output: VEC = A*V */

  VEC[0] = A[0][0]*V[0] + A[0][1]*V[1] + A[0][2]*V[2];
  VEC[1] = A[1][0]*V[0] + A[1][1]*V[1] + A[1][2]*V[2];
  VEC[2] = A[2][0]*V[0] + A[2][1]*V[1] + A[2][2]*V[2];

}



/********************************* RANDOM NUMBER GENERATOR **************************/
double ranf0(gsl_rng *rr)
{
  return gsl_rng_uniform (rr);
}



/********************************* RANDOM NUMBER FROM GAUSSIAN DISTRIBUTION **************************/
double ranfgs(gsl_rng *rr, double sigma)
{

  return gsl_ran_gaussian(rr, sigma);

}






