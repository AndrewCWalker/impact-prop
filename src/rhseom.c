#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include "defs.h"
#include "structs.h"
#include "acceleration.h"
#include "misc.h"

extern int rkcounter;

void rhseom(double state[6], 
	    double et, 
	    double dt, 
	    int n, 
	    double statedot[6], 
	    double doubletime[], 
	    double (*resun)[3], 
	    double (*remoon)[3], 
	    double (*RITJ)[3][3], 
	    double (*CMAT)[MAXGDO],
	    double (*SMAT)[MAXGDO],
	    struct densities_struct *input,
	    int dnp, 
	    struct sat_struct *psatout, 
	    int it,
	    struct msis_struct *pmsis,
	    struct initCd_struct *initCd,
	    double **msismod,
	    int isat,
	    int atmcounter,
	    struct RSM_struct *RSMdata,
	    int nobs) {

  /* this function is the right hand side of the equations of motion */
  /* i.e. statedot = f( state, time) */
  /* t is elapsed seconds from time of et0 */
  /* position and velocity vectors in J2000 frame (in km and km/s) */
  /* et is ephemeris time */
  /* n is the satellite index (usually 1 unless simulating many satellites) */

  int i, j, itime;
  
  double rvec[3], vvec[3], avec[3], nvec[3];
  double rvec_ITRF[3] = {0.0, 0.0, 0.0};
  double vvec_ITRF[3] = {0.0, 0.0, 0.0};

  double RM[3][3], RMT[3][3];

  double cacc[3], sacc[3], dacc[3], nacc[3], tacc[3];

  struct sat_struct *satout;

  satout = psatout + it;

  rkcounter++;

  for(i=0; i<3; i++) {
    rvec[i] = state[i];
    vvec[i] = state[i+3];
    avec[i] = 0.0; /* Acceleration vector */
    cacc[i] = 0.0; /* Central acceleration vector */
    sacc[i] = 0.0; /* SRP acceleration vector */
    dacc[i] = 0.0; /* Drag acceleration vector */
    nacc[i] = 0.0; /* Non-spherical gravity acceleration vector */
    tacc[i] = 0.0; /* Third body acceleration vector */
    nvec[i] = 0.0; /* Temporary Non-spherical vector */
  }
  
  /* Calculate the rotation from earth fixed (use ITRF93 for now) to inertial (use J2000 for now) */
  if (UseSphHarmon || UseDrag) {
        
    itime = findnearest(et, dt, doubletime, dnp);
    for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
	RM[i][j] = RITJ[itime][i][j];
      }
    }

    /* Compute Transpose */
    for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
	RMT[i][j] = RM[j][i];
      }
    }
  
    /* Rotate position into ITRF frame */
    matrix_multiply(RM, rvec, rvec_ITRF);
    matrix_multiply(RM, vvec, vvec_ITRF);
  
  } /* UseSphHarmon || UseDrag */


/********************* Calculate the nonspherical gravity acceleration *******************/
  if (UseSphHarmon) {

    nonsphacc(rvec_ITRF, nacc, CMAT, SMAT);
    matrix_multiply(RMT, nacc, nvec);

    for(i=0; i<3; i++) {
      avec[i] += nvec[i];
      if(rkcounter % RKN == 0) satout->satacc[2][i] = nvec[i];
    }
     
  } /* UseSphHarmon */

/******************* Calculate the Earth point-mass gravity ****************************/
  if (Use2Body) {
        
    centralacc(rvec, cacc);
    for(i=0; i<3; i++) {
      avec[i] += cacc[i];
      if(rkcounter % RKN == 0) satout->satacc[0][i] = cacc[i];
    }

  } /* Use2Body */

/*********************** Calculate the 3rd-body effects from the Moon ******************/
  if (UseMoon) {

    thirdbodyacc(rvec, et, dt, doubletime, resun, remoon, "Moon", tacc, dnp);

    for(i=0; i<3; i++) {
      avec[i] += tacc[i];
      if(rkcounter % RKN == 0) satout->satacc[4][i] = tacc[i];  // FIXME?
    }

  } /* UseMoon */

/*********************** Calculate the 3rd-body effects from the Sun ******************/
  if (UseSun) {

    thirdbodyacc(rvec, et, dt, doubletime, resun, remoon, "Sun", tacc, dnp);

    for(i=0; i<3; i++) {
      avec[i] += tacc[i];
      if(rkcounter % RKN == 0) satout->satacc[5][i] = tacc[i];
    }

  } /* UseSun */   
     
/********* Calculate the SRP acceleration, taking into account Earth eclipse **********/
  if (UseSRP) {

    srpacc(rvec, et, dt, doubletime, resun, sacc, dnp, psatout, it, initCd);
    for(i=0; i<3; i++) {
      avec[i] += sacc[i];
      if(rkcounter % RKN == 0) satout->satacc[3][i] = sacc[i];
    }
   
  } /* UseSRP */
     
/********************** Calculate the atmospheric drag effects *************************/

  if (UseDrag) {
        
    dragacc(state, rvec_ITRF, vvec_ITRF, et, dt, n, input, dacc, psatout, it, RMT, pmsis, initCd, 
	    msismod, isat, atmcounter, RSMdata, nobs);
    for(i=0; i<3; i++) {
      avec[i] += dacc[i];
      if(rkcounter % RKN == 0) satout->satacc[1][i] = dacc[i];
    }

  } /* UseDrag */ 

  statedot[0] = vvec[0];
  statedot[1] = vvec[1];
  statedot[2] = vvec[2];
  statedot[3] = avec[0];
  statedot[4] = avec[1];
  statedot[5] = avec[2];

}




