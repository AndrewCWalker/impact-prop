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
#include "assoclegendre.h"
#include "atmprop.h"
#include "misc.h"
#include "Cd.h"


extern int rkcounter;

void centralacc(double rvec[3], double cacc[3]) {

  /* Calculate the central-body gravitational acceleration, input is position */
  /* vector in inertial frame (km), output is acceleration vector in inertial */
  /* frame (km/s^2). Requires physical constants */

  int i;
  double r;

  r = norm(rvec);

  for(i=0; i<3; i++) {
    cacc[i] = -(JGM3GMe/(r*r*r))*rvec[i];
  }

}

void srpacc(double rvec[3], double et, double dt, double et_doubletime[], double (*r_e_sun_doubletime)[3], double sacc[3], int dnp, struct sat_struct *psatout, int it, struct initCd_struct *initCd) {

  /* Calculate the acceleration vector due to srp in km/s^2, where rvec is */
  /* the position (km) vector of satellite wrt earth center in J2000 */

  double Ls = 3.823e26;                 /*Luminosity of the Sun [W] */
  double rv_es[3], rv_e3[3], rv_sS[3];  /* Position vectors */
  double sval, Psrp; 
  int i;
  int itime;

  struct sat_struct *satout;
  satout = psatout + it;

  for(i=0; i<3; i++) {
    rv_es[i] = rvec[i];                 /* Position vector from earth to satellite */
  } 

  itime = findnearest(et, dt, et_doubletime, dnp);
  for(i=0; i<3; i++) {
    rv_e3[i] = r_e_sun_doubletime[itime][i];
  }

  for(i=0; i<3; i++) {
    rv_sS[i] = rv_e3[i] - rv_es[i];     /* Position vector from satelite to sun */
  }
  
  sval = shadowfunc(rvec, rv_e3); /* Get the shadow value (1 = in sun, 0 = in umbra, 0.0 to 1.0 = penumbra) */
    
  Psrp = Ls/(4.0*M_PI*CSOLV*pow(norm(rv_sS), 2.0))*1.0e-6; /* Account for variation in distance to the sun */
    
  for(i=0; i<3; i++) {
    sacc[i] = - sval * (rv_sS[i]/norm(rv_sS)) * Psrp * initCd->SRP_Cr * initCd->proj_area / initCd->sat_mass;
  }

  if(rkcounter % RKN == 0) {
    satout->Asrp = initCd->proj_area;
    satout->shadowflag = sval;
  }

}


double shadowfunc(double r_J2000[3], double rv_e3[3]) {
    
  /* r_J2000 is position of satellite wrt earth center in J2000 frame in km */
  /* rv_e3 is position vector from earth to sun (in J2000 frame) in km */
    
  /* returns 0.0 for umbra (total eclipse) */
  /* returns 1.0 for sunlight */
  /* returns 0 < value < 1 for penumbra. here = 0.5 */

  /* constant penumbra and umbra angles */
  double alpha_pen = 0.26900424 * DEG2RAD;
  double alpha_umb = 0.26411888 * DEG2RAD;

  double satvert = 0.0;
  double penvert = 0.0;
  double sathoriz = 0.0;
  double umbvert = 0.0;

  double rv_e3_unit[3], r_J2000_unit[3];
  double cos_sun_angle;
  double sval;
  int i;

  /* Get position of sun with respect to the Earth */

  /* positon vector from earth to third body */
  for(i=0; i<3; i++) {
    rv_e3_unit[i] = rv_e3[i]/norm(rv_e3);
    r_J2000_unit[i] = r_J2000[i]/norm(r_J2000);
  }

  if(dot(rv_e3, r_J2000) < 0.0) {

    cos_sun_angle = -1.0*dot(rv_e3_unit, r_J2000_unit);
    sathoriz = norm(r_J2000)*cos_sun_angle;
    satvert = norm(r_J2000)*sin(acos(cos_sun_angle));
    penvert = JGM3Re + tan(alpha_pen)*sathoriz;

    if(satvert <= penvert) { /* In Penumbra or Umbra */
      umbvert = JGM3Re - tan(alpha_umb)*sathoriz; 
      sval = (satvert - umbvert)/(penvert - umbvert);

      if(satvert <= umbvert) { /* In Umbra */              
	sval = 0.0;
      }
    } else {
      sval = 1.0;
    }
  } else {
    sval = 1.0;
  }

  return(sval);
}


void dragacc(double state[6], double r_ITRF[3], double v_ITRF[3], double et, double dt, int n, 
	     struct densities_struct *input, double dacc[3], struct sat_struct *psatout, int it, 
	     double RMT[3][3], struct msis_struct *pmsis, struct initCd_struct *initCd, double **msismod, 
	     int isat, int atmcounter, struct RSM_struct *RSMdata, int nobs) {

  double rvec[3], vvec[3];
  double zhat[3] = {0.0, 0.0, 1.0};
  double rxom[3] = {0.0, 0.0, 0.0};
  double omvecearth[3], vvecrel[3], vJ2000[3];
  double vrel;
  double Cd_T, BC;
  int i;

  double rho = 0.0;
  double temp = 0.0;
  double nden[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  struct sat_struct *satout;
  satout = psatout + it;

  /* calculate the acceleration vector due to drag in km/s^2, where state is */
  /* the position and velocity vector (km and km/s) */

  /* calculate the velocity vector relative to the atmosphere. */
  /* Assume the atmosphere is corotating with the Earth. */
  /* Also, assume the Earth angular velocity vector is simply modeled: */

  for(i=0; i<3; i++) {
    rvec[i] = state[i];
    vvec[i] = state[i+3];
    vvec[i] = v_ITRF[i];
    omvecearth[i] = JGM3ome*zhat[i];
  }

  cross(omvecearth, r_ITRF, rxom);

  for(i=0; i<3; i++) {
    vvecrel[i] = vvec[i] - rxom[i];
  }

  matrix_multiply(RMT, vvecrel, vJ2000);

  vrel = norm(vvecrel);

  /* Calculate the density and other atmospheric properties */  
  atmospheric_properties(r_ITRF, et, input, &rho, &temp, nden, pmsis, msismod, isat, it, atmcounter, nobs);

  if (UseHFCd) {
    Cd_T = CdSetup(vrel, temp, nden, psatout, it, initCd, RSMdata);
  } else { 
    Cd_T =  initCd->constant_Cd;
  } 

       
  BC = Cd_T * initCd->proj_area / initCd->sat_mass;
  
  for(i=0; i<3; i++) {
    dacc[i] = -0.5 * BC * rho * vrel * vJ2000[i];
  }

  if(rkcounter % RKN == 0) {
    satout->rho = rho;
    satout->Cd = Cd_T;
    satout->BC = BC;
    satout->Adrag = initCd->proj_area;

    if (UseHFCd) {
      satout->temp = temp;
      for(i=0; i<6; i++) {
	satout->nden[i] = nden[i];
      }
    } else { 
      satout->temp = 0.0;
      for(i=0; i<6; i++) {
	satout->nden[i] = 0.0;
      }
    } 

  }

}



void nonsphacc(double rvec[3], double nacc[3], double (*CMAT)[MAXGDO], double (*SMAT)[MAXGDO]) {

  /* calculate the nonspherical acceleration in the earth-fixed frame, where */
  /* rvec is the satellite position vector in earth-fixed frame in km */
  /* and N is the degree and order of the gravity field */

  int normflag;
  int n, m;
  int i, j;

  double eps = 1.0e-16;
  double Re, mu;
  double phigc, lambdasat;
  double PMAT[MAXGDO][MAXGDO];
  double r, ri, rj, rk;
  double tansat;
  double dUdr, dUdphi, dUdlam;
  double Pnm, Pn_mplus_1, Cnm, Snm;
  double Rer, Rern;
  double Pi;
  double cml, sml;
  double ai, aj, ak;
  double term1, term2;
  double rijnorm;

  if (EarthGrav == 1) { /* EGM96 */
    Re = EGM96Re;
    mu = EGM96GMe;
    normflag = 1;
  } else { /* JGM3 */
    Re = JGM3Re;
    mu = JGM3GMe;
    normflag = 0;

  } /* EarthGrav */
    
  phigc = atan(rvec[2]/sqrt(rvec[0]*rvec[0] + rvec[1]*rvec[1])); /* Calculate the satellite's geocentric latitude */
  lambdasat = atan2(rvec[1], rvec[0]);                           /* Calculate the satellite's geocentric longitude */

  /* Initialize Legendre Polynomials */
  for(i=0; i<MAXGDO; i++) {
    for(j=0; j<MAXGDO; j++) {
      PMAT[i][j] = 0.0;
    }
  }

  /* Calculate associated Legendre functions based on satellite geocentric latitude */
  /* NOTE: Need to calculate N+1 deg and order because need Pn_mplus1 down below */
  assoclegendre(phigc, PMAT);

  r = norm(rvec);
  ri = rvec[0];
  rj = rvec[1];
  rk = rvec[2];
  Rer = Re/r;

  tansat = tan(phigc);

  /* Calculate the partial derivative of the nonspherical potential (Eq.8-19) w.r.t. radius (Eq.8-25) */
  dUdr = 0.0;
  dUdphi = 0.0;
  dUdlam = 0.0;

  for(n=2; n<GDO+1; n++) {
    for(m=0; m<n+1; m++) {
      Pnm = PMAT[n][m];
      Pn_mplus_1 = PMAT[n][m+1];
      Cnm = CMAT[n][m];
      Snm = SMAT[n][m];

      if(normflag==1) {
	/* Gravity coefficients are normalized, so must also normalize the associated Legendre polynomials */
	Pi = normcoef(n,m);
	Cnm = Cnm/Pi;
	Snm = Snm/Pi;
      }

      Rern = pow(Rer, n);
      
      /* Note: If I wanted to speed things up even more, I could use recursive relations on these trig functions */
      cml = cos(m*lambdasat);
      sml = sin(m*lambdasat);
            
      dUdr =   dUdr   + Rern * (n+1.0) * Pnm * (Cnm*cml + Snm*sml);         
      dUdphi = dUdphi + Rern * (Pn_mplus_1 - m*tansat*Pnm)*(Cnm*cml + Snm*sml);       
      dUdlam = dUdlam + Rern * m * Pnm * (Snm*cml - Cnm*sml);

    }
  }
 
  dUdr = (-mu/(r*r))*dUdr;
  dUdphi = (mu/r)*dUdphi;
  dUdlam = (mu/r)*dUdlam;

  rijnorm = sqrt(ri*ri + rj*rj);

  /* Note: I don't think Vallado warns about the case where you are exactly at */
  /* phi = 90 or -90, and this rijnorm term is zero, hence 1/0... */
  if(rijnorm < eps) {
    term1 = 0.0;
    term2 = 0.0;
  } else {
    term1 = (dUdr/r - dUdphi*rk/(r*r*rijnorm));
    term2 = (dUdlam/(ri*ri + rj*rj));
  }

  ai = term1*ri - term2*rj;
  aj = term1*rj + term2*ri;
  ak = dUdr*rk/r + dUdphi*rijnorm/(r*r);

  nacc[0] = ai;
  nacc[1] = aj;
  nacc[2] = ak;

}


void thirdbodyacc(double rvec[3], double et, double dt, double et_doubletime[], double (*r_e_sun_doubletime)[3], 
		  double (*r_e_moon_doubletime)[3], char thirdbodyname[], double tacc[3], int dnp) {
    
  /* Note: uses simple but potentially numerically unstable formula (8-35 in Vallado).  */
  /* Do not use for non-Earth-orbiting satellite. */
    
  /* rvec is satellite position vector in J2000 frame (km) */
  /* et is ephemeris time */
  /* thirdbodyname is a string recognized by JPL SPICE for third body (e.g. 'MOON' or 'SUN'). */

  double rv_es[3], rv_e3[3], rv_s3[3];
  double norm_s3, norm_e3;
  double mu3 = 0.0;
  int i;  
  int itime;

  for(i=0; i<3; i++) {
    rv_es[i] = rvec[i]; /* Position vector from earth to satellite */
  }
  
  itime = findnearest(et, dt, et_doubletime, dnp); /* Positon vector from earth to third body */

  if(strcmp(thirdbodyname,"Sun")==0) {
    mu3 = SunGM;
    for(i=0; i<3; i++) {
      rv_e3[i] = r_e_sun_doubletime[itime][i];
    }
  } else if(strcmp(thirdbodyname,"Moon")==0) {
    mu3 = MoonGM;
    for(i=0; i<3; i++) {
      rv_e3[i] = r_e_moon_doubletime[itime][i];
    }
  } else {
    printf("Unknown 3rd Body\n");
  }

  for(i=0; i<3; i++) {
    rv_s3[i] = rv_e3[i] - rv_es[i]; /* Position vector from satelite to third body */
  }

  norm_s3 = norm(rv_s3);
  norm_e3 = norm(rv_e3);

  for(i=0; i<3; i++) {
    tacc[i] = mu3*(rv_s3[i]/pow(norm_s3, 3) - rv_e3[i]/pow(norm_e3, 3));
  }

}
























