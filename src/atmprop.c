#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include "defs.h"
#include "structs.h"
#include "atmprop.h"
#include "misc.h"
#include "interpolation.h"
#include "nrlmsise-00.h"
#if !DIORAMA
#include "SpiceUsr.h"
#endif
#include "transformation.h"
#include "msisinputs.h"

void atmospheric_properties(double r_ITRF[3], double et, struct densities_struct *input, double *rho, 
			    double *temp, double nden[6], struct msis_struct *pmsis, double **msismod, 
			    int isat, int it, int atmcounter, int nobs) {

  /* Calculate atmospheric properties at satellite position */
    
  /* Inputs: */
  /* r_ITRF is the satellite position in ITRF frame in km   */
  /* et is ephemeris time                                   */
  /* density rho is returned in kg/km^3                     */

  double mass[MAXSPECIES];     /* [O, O2, N, N2, He, H] */
  double theta = 0.0;
  double phi = 0.0;
  double R = 0.0;
  double Re;
  double geocentric_Long, geocentric_Lat;
  double geodetic_Long, geodetic_Lat;
  double rdsat;
  double alpha, delta;
  double old_Lat = 1.0e10;;
  double Cv = 0.0;
  double alt;

  /* Define atomic/molecular masses of each species */
  mass[0] = 2.657e-26; /* atomic oxygen */
  mass[1] = 5.314e-26; /* diatomic oxygen */
  mass[2] = 2.326e-26; /* atomic nitrogen */
  mass[3] = 4.652e-26; /* diatomic nitrogen */
  mass[4] = 6.646e-27; /* helium */
  mass[5] = 1.674e-27; /* hydrogen */
     
  Re = JGM3Re;

  /* Convert to spherical coordinates*/
  cart2sph(r_ITRF[0], r_ITRF[1], r_ITRF[2], &R, &theta, &phi);

  /* Convert to convention used to define density profiles: */
  geocentric_Long = theta;
  if(geocentric_Long < 0) {
    geocentric_Long = geocentric_Long + 2*M_PI;
  }

  geocentric_Long *= RAD2DEG;
  geocentric_Lat = phi*RAD2DEG;
  
  /* Because GITM altitude is above a spherical Earth */
  geodetic_Long = geocentric_Long;

  rdsat = sqrt(pow(r_ITRF[0], 2.0) + pow(r_ITRF[1], 2.0));
  alpha = asin(r_ITRF[1]/rdsat);
  delta = atan(r_ITRF[2]/rdsat);

  geodetic_Lat = delta;

  while(fabs(geodetic_Lat-old_Lat) > 1e-3) {
    Cv = Re/sqrt(1.0 - pow(EECC*sin(geodetic_Lat), 2));
    old_Lat = geodetic_Lat;
    geodetic_Lat = atan((r_ITRF[2] + Cv*pow(EECC, 2.0)*sin(old_Lat))/rdsat);
  }
        
  alt = rdsat/cos(geodetic_Lat)-Cv;
  geodetic_Lat *= RAD2DEG;

  if (DensityModel==1 && DynamicMSIS==1) { /* MSIS on the fly */

    densitydynamicMSIS(alt, et, geodetic_Lat, geodetic_Long, rho, temp, nden, pmsis, msismod, isat, nobs);

  } else { /* Not MSIS on the fly */
    int k=0;
    int i;

    if (DensityModel==1 || DensityModel==2) { /* MSIS or GITM */

      if (DensityTimeLookup == 0) { /* Find nearest atmospheric data set */

	/* Compute time difference to nearest data points */
	double td1 = fabs(et - input->etlist[atmcounter-1]);
	double td2 = fabs(et - input->etlist[atmcounter]);

	/* Find closest data point */
	if(td2 > td1) {
	  k = 0;
	} else {
	  k = 1;
	}
                  
      } else { /* Linear interpolation in time  - FIXME BROKEN!!!!!!*/

	/* Compute data point weights */
	double w1 = fabs(et - input->etlist[atmcounter])/(input->etlist[1]-input->etlist[0]);
	double w2 = fabs(et - input->etlist[atmcounter+1])/(input->etlist[1]-input->etlist[0]);

	/* Find closest data point */
	if(w2 > w1) {
	  k = 0;
	} else {
	  k = 1;
	}
             
      } /* DensityTimeLookUp */
    } /* DensityModel */

    if (DensityModel == 0) { /* CIRA72 */
        
      double hellp = hellipsoid(r_ITRF);
      *rho = densitycira72(hellp);

    } else if (DensityModel == 3) { /* US STANDARD ATMOSPHERE */
#if 0
      double hellp = hellipsoid(r_ITRF);
      *rho = densityussa76(hellp, input->AltKm, input->Rho);
#endif
    } else { /* DensityModel */

         /* Check if data is beyond upper domain boundary */
      if(alt > input->AltKm[NVERT-1]) {

	if (UseHFCd) {

	  if(DensityModel==4) { /* HASDM */	    
	    dynamicMSISforHASDM(alt, et, geodetic_Lat, geodetic_Long, temp, nden, pmsis, msismod, isat, nobs);  
	  } else { /* MSIS or GITM */
	    extrapolate_atm_properties(k, geodetic_Lat, geodetic_Long, alt, input->LatDeg, input->LonDeg, 
				       input->AltKm, input->Rho, 
				       input->Temperature, input->ndenvec, &rho, &temp, nden, mass, it);

	  }
      
	} else { /* UseHFCd */
                
	  double hellp = hellipsoid(r_ITRF);
	  *rho = densitycira72(hellp);

	} /* UseHFCd */

      } else {

	/* Interpolate density in the latitudinal and longitudinal directions */
	*rho = interpLatLong(geodetic_Lat, geodetic_Long, alt, input->Rho[k], input->LatDeg, 
			     input->LonDeg, input->AltKm, it);

	if (UseHFCd) {

	  if(DensityModel==4) { /* HASDM */
	    
	    dynamicMSISforHASDM(alt, et, geodetic_Lat, geodetic_Long, temp, nden, pmsis, msismod, isat, nobs);
	    
	  } else { /* GITM or Non-dynamic MSIS */
	    
	    *temp = interpLatLong(geodetic_Lat, geodetic_Long, alt, input->Temperature[k], input->LatDeg, 
				  input->LonDeg, input->AltKm, it);
	    for(i=0; i<NSPECIES; i++) {
	      nden[i] = interpLatLong(geodetic_Lat, geodetic_Long, alt, input->ndenvec[i][k], 
				      input->LatDeg, input->LonDeg, input->AltKm, it);
	    }
	    
	  }
         
	} /* UseHFCd */     

      } /* Above / Below Upper Domain Boundary */ 
            
    } /* DensityModel */ 

  } /* MSIS on the fly */
        
}


double densitycira72(double h) {

  /* Uses table from pp. 537 (Table 8-4) in Vallado for simple exponential */
  /* atmospheric density model. */
  
  /* h is the height above ellipsoid in km. */
  /* density p is returned in kg/km^3 */
  /* H is the scale height in km. */
  
  /* these numbers for density units of kg/m^3, conversion to kg/km^3 is at */
  /* the end of the function. */

  double h0, p0, H, p;

  if(h >= 0 && h < 25) {
    h0 = 0.0;
    p0 = 1.225;
    H = 7.249;
  } else if(h >= 25 && h < 30) {
    h0 = 25.0;
    p0 = 3.899e-2;
    H = 6.349;
  } else if(h >= 30 && h < 40) {
    h0 = 30.0;
    p0 = 1.774e-2;
    H = 6.682;
  } else if(h >= 40 && h < 50) {
    h0 = 40.0;
    p0 = 3.972e-3;
    H = 7.554;
  } else if(h >= 50 && h < 60) {
    h0 = 50.0;
    p0 = 1.057e-3;
    H = 8.382;
  } else if(h >= 60 && h < 70) {
    h0 = 60.0;
    p0 = 3.206e-4;
    H = 7.714;
  } else if(h >= 70 && h < 80) {
    h0 = 70.0;
    p0 = 8.770e-5;
    H = 6.549;
  } else if(h >= 80 && h < 90) {
    h0 = 80.0;
    p0 = 1.905e-5;
    H = 5.799;
  } else if(h >= 90 && h < 100) {
    h0 = 90.0;
    p0 = 3.396e-6;
    H = 5.382;
  } else if(h >= 100 && h < 110) {
    h0 = 100.0;
    p0 = 5.297e-7;
    H = 5.877;
  } else if(h >= 110 && h < 120) {
    h0 = 110;
    p0 = 9.661e-8;
    H = 7.263;
  } else if(h >= 120 && h < 130) {
    h0 = 120.0;
    p0 = 2.438e-8;
    H = 9.473;
  } else if(h >= 130 && h < 140) {
    h0 = 130.0;
    p0 = 8.484e-9;
    H = 12.636;
  } else if(h >= 140 && h < 150) {
    h0 = 140.0;
    p0 = 3.845e-9;
    H = 16.149;
  } else if(h >= 150 && h < 180) {
    h0 = 150.0;
    p0 = 2.070e-9;
    H = 22.523;
  } else if(h >= 180 && h < 200) {
    h0 = 180.0;
    p0 = 5.464e-10;
    H = 29.740;
  } else if(h >= 200 && h < 250) {
    h0 = 200.0;
    p0 = 2.789e-10;
    H = 37.105;
  } else if(h >= 250 && h < 300) {
    h0 = 250.0;
    p0 = 7.248e-11;
    H = 45.546;
  } else if(h >= 300 && h < 350) {
    h0 = 300.0;
    p0 = 2.418e-11;
    H = 53.628;
  } else if(h >= 350 && h < 400) {
    h0 = 350.0;
    p0 = 9.518e-12;
    H = 53.298;
  } else if(h >= 400 && h < 450) {
    h0 = 400.0;
    p0 = 3.725e-12;
    H = 58.515;
  } else if(h >= 450 && h < 500) {
    h0 = 450.0;
    p0 = 1.585e-12;
    H = 60.828;
  } else if(h >= 500 && h < 600) {
    h0 = 500.0;
    p0 = 6.967e-13;
    H = 63.822;
  } else if(h >= 600 && h < 700) {
    h0 = 600.0;
    p0 = 1.454e-13;
    H = 71.835;
  } else if(h >= 700 && h < 800) {
    h0 = 700.0;
    p0 = 3.614e-14;
    H = 88.667;
  } else if(h >= 800 && h < 900) {
    h0 = 800.0;
    p0 = 1.170e-14;
    H = 124.64;
  } else if(h >= 900 && h < 1000) {
    h0 = 900.0;
    p0 = 5.245e-15;
    H = 181.05;
  } else {
    h0 = 1000.0;
    p0 = 3.019e-15;
    H = 268.00;
  }
          
  p = p0 * exp((h0-h)/H);
    
  /* Convert density to kg/km^3 */
  p *= pow(1000.0, 3);

  return(p);

}


double densityussa76(double h, double altitude[MAXVERT],
		     double density[MAXLON][MAXLAT][MAXVERT]) {

  /* Uses 1976 US Standard Atmosphere Look-up Table */
  
  /* h is the height above ellipsoid in km. */
  /* density p is returned in kg/km^3 */
  
  /* these numbers for density units of kg/m^3, conversion to kg/km^3 is at */
  /* the end of the function. */

  int iAlt;
  double irho = 0.0;
  double y1, y2, x1, x2;
  double H;

  iAlt = 0;

  /* Find correct altitude index */
  while(h > altitude[iAlt]) {
    iAlt++;
    /* If altitude is above 1000 km, use last two cells to calculate effective scale height */
    if(iAlt==NVERT) {
      iAlt--;
    }
  }
    
  /* Assign neighboring densities and altitudes to dummy variables */
  y1 = density[0][0][iAlt-1];
  y2 = density[0][0][iAlt];
  x1 = altitude[iAlt-1];
  x2 = altitude[iAlt];
  
  /* Calculate effective scale height */
  H = -(x2-x1)/log(y2/y1);
  
  /* Compute density at satellite altitude */
  irho = y1*exp(-(h-x1)/H);
 
  /* Convert density to kg/km^3 */
  irho *= pow(1.0e3, 3);

  return(irho);

}

void densitydynamicMSIS(double alt, double et, double geodetic_Lat, double geodetic_Long, double *rho, double *temp, double nden[6], struct msis_struct *pmsis, double **msismod, int isat, int nobs) {

  struct nrlmsise_output MSISoutput;
  struct nrlmsise_input MSISinput;
  struct nrlmsise_flags MSISflags;
  struct ap_array MSISaph;
  struct msis_struct *msis_inputs;
  char utcmsis[24];
  char msisyear[5], msisdoy[4], msishour[3], msismin[3], msissec[7];
  double year, doy, hour, min, sec, UTCsec;
  double hr2sec = 3600.0;
  double min2sec = 60.0;
  double msis_ap_array[7];

  int i;
  int imsis;

  for(i=0; i<7; i++) {
    msis_ap_array[i] = 0.0;
  }

#if !DIORAMA
  /* Compute MSIS time parameters from et */
  et2utc_c(et, "D", 10, 24, utcmsis);
#endif

  /* Copy pieces of utcmsis to individual strings */
  memcpy(msisyear, utcmsis,   4*sizeof(char));
  memcpy(msisdoy,  utcmsis+5, 3*sizeof(char));
  memcpy(msishour,  utcmsis+12,  2*sizeof(char));
  memcpy(msismin,   utcmsis+15, 2*sizeof(char));
  memcpy(msissec,   utcmsis+18, 6*sizeof(char));

  /* Terminate string since memcpy does not */
  msisyear[4] = '\0';
  msisdoy[3] = '\0';
  msishour[2] = '\0';
  msismin[2] = '\0';
  msissec[6] = '\0';

  /* Convert strings to floats */
  year = atof(msisyear);
  doy = atof(msisdoy);
  hour = atof(msishour);
  min = atof(msismin);
  sec = atof(msissec);

  /* Convert hour, min, sec to UTC seconds */
  UTCsec = hour*hr2sec + min*min2sec + sec;

  /* Find msis index */
  imsis = find_msis_time(nobs, year, doy);
  msis_inputs = pmsis + imsis;

  /* Compute proper ap index format for msis */
  get_msis_ap_array(pmsis, imsis, hour, msis_ap_array);

  /* Define MSIS inputs */
  MSISinput.year = (int)year;
  MSISinput.doy = (int)doy;
  MSISinput.sec = UTCsec;
  MSISinput.alt = alt;
  MSISinput.g_lat = geodetic_Lat;
  MSISinput.g_long = geodetic_Long;
  MSISinput.lst = UTCsec/hr2sec + MSISinput.g_long/15.0;
  if(StateEnsembles && !ReadSW) {
    MSISinput.f107 = msismod[isat][0];
    MSISinput.f107A = msismod[isat][1];
    MSISinput.ap = msismod[isat][2];
    for(i=0; i<7; i++) {
      MSISaph.a[i] = msismod[isat][2];
    }
  } else {
    MSISinput.f107 = msis_inputs->F107;
    MSISinput.f107A = msis_inputs->F107A;
    MSISinput.ap = msis_inputs->ap_index[8];
    for(i=0; i<7; i++) {
      if(ReadSW) {
	MSISaph.a[i] = msis_ap_array[i];
      } else {
	MSISaph.a[i] = msis_inputs->ap_index[8];
      }
    }
  }
  MSISinput.ap_a = &MSISaph;

  /* Define MSIS switches - See nrlmsise-00.h */
  for(i=0; i<24; i++) {
    MSISflags.switches[i]=1;
  }
  
  if(ReadSW) {
    MSISflags.switches[9] = -1; /* Set to -1 to read ap_array */
  }

  /* Calculate MSIS density at satellite location */
  gtd7(&MSISinput, &MSISflags, &MSISoutput);

  /* Convert density to kg/km^3 */
  *rho = MSISoutput.d[5]*1.0e9;
  *temp = MSISoutput.t[0];

  /* Compute temperature and composition for high-fidelity drag coefficient */
  if (UseHFCd) {
    *temp = MSISoutput.t[0];
    nden[0] = MSISoutput.d[1];
    nden[1] = MSISoutput.d[3];
    nden[2] = MSISoutput.d[7];
    nden[3] = MSISoutput.d[2];
    nden[4] = MSISoutput.d[0];
    nden[5] = MSISoutput.d[6];
  } /* UseHFCd */

}





void dynamicMSISforHASDM(double alt, double et, double geodetic_Lat, double geodetic_Long, double *temp, double nden[6], struct msis_struct *pmsis, double **msismod, int isat, int nobs) {

  struct nrlmsise_output MSISoutput;
  struct nrlmsise_input MSISinput;
  struct nrlmsise_flags MSISflags;
  struct ap_array MSISaph;
  struct msis_struct *msis_inputs;
  char utcmsis[24];
  char msisyear[5], msisdoy[4], msishour[3], msismin[3], msissec[7];
  double year, doy, hour, min, sec, UTCsec;
  double hr2sec = 3600.0;
  double min2sec = 60.0;
  double msis_ap_array[7];

  int i;
  int imsis;

  for(i=0; i<7; i++) {
    msis_ap_array[i] = 0.0;
  }

#if !DIORAMA
  /* Compute MSIS time parameters from et */
  et2utc_c(et, "D", 10, 24, utcmsis);
#endif

  /* Copy pieces of utcmsis to individual strings */
  memcpy(msisyear, utcmsis,   4*sizeof(char));
  memcpy(msisdoy,  utcmsis+5, 3*sizeof(char));
  memcpy(msishour,  utcmsis+12,  2*sizeof(char));
  memcpy(msismin,   utcmsis+15, 2*sizeof(char));
  memcpy(msissec,   utcmsis+18, 6*sizeof(char));

  /* Terminate string since memcpy does not */
  msisyear[4] = '\0';
  msisdoy[3] = '\0';
  msishour[2] = '\0';
  msismin[2] = '\0';
  msissec[6] = '\0';

  /* Convert strings to floats */
  year = atof(msisyear);
  doy = atof(msisdoy);
  hour = atof(msishour);
  min = atof(msismin);
  sec = atof(msissec);

  /* Convert hour, min, sec to UTC seconds */
  UTCsec = hour*hr2sec + min*min2sec + sec;

  /* Find msis index */
  imsis = find_msis_time(nobs, year, doy);
  msis_inputs = pmsis + imsis;

  /* Compute proper ap index format for msis */
  get_msis_ap_array(pmsis, imsis, hour, msis_ap_array);

  /* Define MSIS inputs */
  MSISinput.year = (int)year;
  MSISinput.doy = (int)doy;
  MSISinput.sec = UTCsec;
  MSISinput.alt = alt;
  MSISinput.g_lat = geodetic_Lat;
  MSISinput.g_long = geodetic_Long;
  MSISinput.lst = UTCsec/hr2sec + MSISinput.g_long/15.0;
  if(StateEnsembles && !ReadSW) {
    MSISinput.f107 = msismod[isat][0];
    MSISinput.f107A = msismod[isat][1];
    MSISinput.ap = msismod[isat][2];
    for(i=0; i<7; i++) {
      MSISaph.a[i] = msismod[isat][2];
    }
  } else {
    MSISinput.f107 = msis_inputs->F107;
    MSISinput.f107A = msis_inputs->F107A;
    MSISinput.ap = msis_inputs->ap_index[8];
    for(i=0; i<7; i++) {
      MSISaph.a[i] = msis_ap_array[i];
    }
  }
  MSISinput.ap_a = &MSISaph;

  /* Define MSIS switches - See nrlmsise-00.h */
  for(i=0; i<24; i++) {
    MSISflags.switches[i]=1;
  }
  MSISflags.switches[9] = -1; /* Set to -1 to read ap_array */

  /* Calculate MSIS density at satellite location */
  gtd7(&MSISinput, &MSISflags, &MSISoutput);

  /* Compute temperature and composition for high-fidelity drag coefficient */
  *temp = MSISoutput.t[0];
  nden[0] = MSISoutput.d[1];
  nden[1] = MSISoutput.d[3];
  nden[2] = MSISoutput.d[7];
  nden[3] = MSISoutput.d[2];
  nden[4] = MSISoutput.d[0];
  nden[5] = MSISoutput.d[6];

}
