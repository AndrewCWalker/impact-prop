#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include "defs.h"
#include "interpolation.h"


double interpLatLong(double lati, double loni, double alt, double parameter[MAXLON][MAXLAT][MAXVERT], 
		     double LatDeg[MAXLAT], double LonDeg[MAXLON], double AltKm[MAXVERT], int it) {

  /* Interpolate atmospheric properties in latitude and longitude */
  /* lat and lon in degrees */
  /* Constant contains latitude, longitude, and altitude vectors */
  /* parameter is the variable that is interpolated e.g. rho, T, etc. */
  /* DensityModel is GITM or MSIS */

  int LatIndxMax, LonIndxMax;
  int iLat, iLon, i;

  double delta = LonDeg[1] - LonDeg[0]; /* Latitude and Longitude have 5 or 10 degree increments */
  double vi = 0.0;
  double LatMin, LonMin;
  double LatMax, LonMax;
  double param[4];

  for(i=0; i<4; i++) {
    param[i] = 0.0;
  }

  if (DensityModel == 2 || DensityModel == 4) { /* GITM or HASDM */

    /* Minimum Latitude and Longitude data points */
    LatMin = LatDeg[0];
    LonMin = LonDeg[0];
        
    /* Maximum Latitude and Longitude index points */
    LatIndxMax = NLAT-1;
    LonIndxMax = NLON-1;

    LatMax = LatDeg[LatIndxMax];
    LonMax = LonDeg[LonIndxMax];
        
    iLat = (int)((lati - LatMin)/delta);
    iLon = (int)((loni - LonMin)/delta);

    if((loni < LonMax && loni > LonMin) && (lati < LatMax && lati > LatMin)) {
	
      param[0] = interpVert(iLat, iLon, alt, AltKm, parameter, it);
      param[1] = interpVert(iLat+1, iLon, alt, AltKm, parameter, it);
      param[2] = interpVert(iLat, iLon+1, alt, AltKm, parameter, it);
      param[3] = interpVert(iLat+1, iLon+1, alt, AltKm, parameter, it);
      
      vi = BilinearInterp(LatDeg, LonDeg, lati, loni, param, iLat, iLon, it);
  
    } else if(lati < LatMin) { /* Average all polar points if past polar end points */
      for(i=0; i<NLON; i++) {
	vi += interpVert(0, i, alt, AltKm, parameter, it);  
      }
      vi = vi/(double)NLON;

    } else if(lati > LatMax) { /* Average all polar points if past polar end points */
      for(i=0; i<NLON; i++) {
	vi += interpVert(LatIndxMax, i, alt, AltKm, parameter, it);  
      }
      vi = vi/(double)NLON;

    } else if (loni < LonMin || loni > LonMax) {   /* Average first and last longitude points */
      param[0] = interpVert(iLat, 0, alt, AltKm, parameter, it);   
      param[1] = interpVert(iLat+1, 0, alt, AltKm, parameter, it);  
      param[2] = interpVert(iLat, LonIndxMax, alt, AltKm, parameter, it); 
      param[3] = interpVert(iLat+1, LonIndxMax, alt, AltKm, parameter, it);  
      
      double x1 = LonDeg[0];
      double x2 = LonDeg[LonIndxMax];
      double y1 = LatDeg[iLat];
      double y2 = LatDeg[iLat+1];
      
      double xp = loni;
      double yp = lati;
      
      double f11 = param[0];
      double f12 = param[1];
      double f21 = param[2];
      double f22 = param[3];
      
      vi = (1.0/((x2-x1)*(y2-y1)))*(f11*(x2-xp)*(y2-yp)+f21*(xp-x1)*(y2-yp)+f12*(x2-xp)*(yp-y1)+f22*(xp-x1)*(yp-y1));
    }

  } else { /* DensityModel == 1 */

    LatMin = LatDeg[0]; 
    LonMin = LonDeg[0];

    LatIndxMax = length(LatDeg)-1;
    LonIndxMax = length(LonDeg)-1;

    iLat = (int)((lati - LatMin)/delta);
    iLon = (int)((loni - LonMin)/delta);
  
    param[0] = interpVert(iLat, iLon, alt, AltKm, parameter, it);
    param[1] = interpVert(iLat+1, iLon, alt, AltKm, parameter, it);
    param[2] = interpVert(iLat, iLon+1, alt, AltKm, parameter, it);
    param[3] = interpVert(iLat+1, iLon+1, alt, AltKm, parameter, it);

    vi = BilinearInterp(LatDeg, LonDeg, lati, loni, param, iLat, iLon, it);

  } /* DensityModel */

  return(vi);

}


double BilinearInterp(double LatDeg[MAXLAT], double LonDeg[MAXLON], double lati, double loni, double param[], int iLat, int iLon, int it) {

  double x1, x2, y1, y2;
  double xp, yp;
  double f11, f12, f21, f22;
  double vi;

  x1 = LonDeg[iLon];
  x2 = LonDeg[iLon+1];
  y1 = LatDeg[iLat];
  y2 = LatDeg[iLat+1];
        
  xp = loni;
  yp = lati;

  f11 = param[0];
  f12 = param[1];
  f21 = param[2];
  f22 = param[3];
        
  vi = (1.0/((x2-x1)*(y2-y1)))*(f11*(x2-xp)*(y2-yp)+f21*(xp-x1)*(y2-yp)+f12*(x2-xp)*(yp-y1)+f22*(xp-x1)*(yp-y1));

  return(vi);

}



double interpVert(int iLat, int iLon, double alt, double AltKm[MAXVERT], double parameter[MAXLON][MAXLAT][MAXVERT], int it) {

  /* Interpolate atmospheric properties in vertical direction */
  /* iLat and iLon determine satellite position in latitude and longitude */
  /* alt is the satellite altitude */
  /* Constant contains latitude, longitude, and altitude vectors */
  /* parameter is the variable to be interpolated e.g. rho, T, etc. */
  
  int iAlt;
  double irho = 0.0;
  double y1, y2, x1, x2;
  double H;
  
  if(alt < AltKm[NVERT-1]) {
    
    iAlt = 1;
    while(alt > AltKm[iAlt]) {
      iAlt++;
    }

    y1 = parameter[iLon][iLat][iAlt-1];
    y2 = parameter[iLon][iLat][iAlt];
    x1 = AltKm[iAlt-1];
    x2 = AltKm[iAlt];
    
    H = -(x2-x1)/log(y2/y1);
    
    irho = y1*exp(-(alt-x1)/H);
    
  }
  
  return(irho);
  
}




void extrapolate_atm_properties(int k, double lati, double loni, double alt, double LatDeg[MAXLAT], 
				double LonDeg[MAXLON], double AltKm[MAXVERT], double Rho[2][MAXLON][MAXLAT][MAXVERT], 
				double Temperature[2][MAXLON][MAXLAT][MAXVERT], 
				double ndenvec[NSPECIES][2][MAXLON][MAXLAT][MAXVERT], double **rho, double **temp, 
				double nden[6], double mass[NSPECIES], int it) {


  double kB = 1.3806488e-3;  /* Boltzmann constant */
  double hdelta = LonDeg[1] - LonDeg[0];
  /* Minimum Latitude and Longitude */
  double LatMin = LatDeg[0];
  double LonMin = LonDeg[0];
  /* Maximum Latitude and Longitude index points */
  int LatIndxMax = NLAT-1;
  int LonIndxMax = NLON-1;
  /* Maximum Latitude and Longitude */
  double LatMax = LatDeg[LatIndxMax];
  double LonMax = LonDeg[LonIndxMax];
  int iLat = (int)((lati - LatMin)/hdelta);
  int iLon = (int)((loni - LonMin)/hdelta);  
  double rho_bilinear[4];
  double nden_bilinear[6][4];
  
  /* Assume temperature is constant above top of the domain */
  double temperature = Temperature[k][iLon][iLat][NVERT-1];
  double gz = 9.80065*pow(JGM3Re/(alt+JGM3Re), 2);     /* Gravity at altitude */
  double e_alt = (alt - AltKm[NVERT-1])*1.0e3; /* Altitude above highest cell in meters */
  double Hz = gz*e_alt/(kB*temperature);

  double Ntot = 0.0;
  double m_avg = 0.0;

  int i, plat, ispec, ilon;

  if(DensityModel==2) { /* GITM */
    
    if((loni < LonMax && loni > LonMin) && (lati < LatMax && lati > LatMin)) {

      for(i=0; i<NSPECIES; i++) {
	
	/* Assume scale height also remains constant */
	nden_bilinear[i][0] = ndenvec[i][k][iLon][iLat][NVERT-1]*exp(-mass[i]*Hz);
	nden_bilinear[i][1] = ndenvec[i][k][iLon][iLat+1][NVERT-1]*exp(-mass[i]*Hz);
	nden_bilinear[i][2] = ndenvec[i][k][iLon+1][iLat][NVERT-1]*exp(-mass[i]*Hz);
	nden_bilinear[i][3] = ndenvec[i][k][iLon+1][iLat+1][NVERT-1]*exp(-mass[i]*Hz);
	
	nden[i] = BilinearInterp(LatDeg, LonDeg, lati, loni, nden_bilinear[i], iLat, iLon, it);
	
      }
      
      /* Calculate the total number density */
      for(i=0; i<NSPECIES; i++) {
	Ntot += nden[i];
      }
      
      /* Calculate average mass */ 
      for(i=0; i<NSPECIES; i++) {
	m_avg += nden[i]*mass[i]/Ntot;
      }
      
      rho_bilinear[0] = Rho[k][iLon][iLat][NVERT-1]*exp(-m_avg*Hz);
      rho_bilinear[1] = Rho[k][iLon][iLat+1][NVERT-1]*exp(-m_avg*Hz);
      rho_bilinear[2] = Rho[k][iLon+1][iLat][NVERT-1]*exp(-m_avg*Hz);
      rho_bilinear[3] = Rho[k][iLon+1][iLat+1][NVERT-1]*exp(-m_avg*Hz);
      
      **rho = BilinearInterp(LatDeg, LonDeg, lati, loni, rho_bilinear, iLat, iLon, it);
      
      **temp = temperature;
      
    } else if(lati < LatMin || lati > LatMax) { /* Average all polar points if past polar end points */

      if(lati<LatMin)
	plat = 0;
      if(lati>LatMax)
	plat = NLAT-1;
      
      for(ilon=0; ilon<NLON; ilon++) {
	
	for(ispec=0; ispec<NSPECIES; ispec++) {
	  nden[ispec] += ndenvec[ispec][k][ilon][plat][NVERT-1]*exp(-mass[ispec]*Hz);
	}

	/* Calculate the total number density */
	for(ispec=0; ispec<NSPECIES; ispec++) {
	  Ntot += nden[ispec];
	}
	
	/* Calculate average mass */
	for(ispec=0; ispec<NSPECIES; ispec++) {
	  //printf("i = %d m_avg = %e nden[i] = %e mass[i] = %e Ntot = %e\n", i, m_avg, nden[ispec], mass[ispec], Ntot);
	  m_avg += nden[ispec]*mass[ispec]/Ntot;
	}
	
	**rho += Rho[k][ilon][plat][NVERT-1]*exp(-m_avg*Hz);
	
      }
      
      /* Divide by number of longitude bins */
      for(ispec=0; ispec<NSPECIES; ispec++) {
	nden[ispec] /= (double)NLON;
      }

      **rho /= (double)NLON;
      **temp = temperature;
      
    } else if (loni < LonMin || loni > LonMax) {   /* Average first and last longitude points */

      double x1 = LonDeg[0];
      double x2 = LonDeg[LonIndxMax];
      double y1 = LatDeg[iLat];
      double y2 = LatDeg[iLat+1];
      
      double xp = loni;
      double yp = lati;
      
      double f11, f12, f21, f22;
      
      for(ispec=0; ispec<NSPECIES; ispec++) {
	
	f11 = ndenvec[ispec][k][0][iLat][NVERT-1]*exp(-mass[ispec]*Hz);
	f12 = ndenvec[ispec][k][0][iLat+1][NVERT-1]*exp(-mass[ispec]*Hz);
	f21 = ndenvec[ispec][k][LonIndxMax][iLat][NVERT-1]*exp(-mass[ispec]*Hz);
	f22 = ndenvec[ispec][k][LonIndxMax][iLat+1][NVERT-1]*exp(-mass[ispec]*Hz);

	nden[ispec] = (1.0/((x2-x1)*(y2-y1)))*(f11*(x2-xp)*(y2-yp)+f21*(xp-x1)*(y2-yp)+f12*(x2-xp)*(yp-y1)+f22*(xp-x1)*(yp-y1));

      }

      /* Calculate the total number density */
      for(i=0; i<NSPECIES; i++) {
	Ntot += nden[i];
      }
      
      /* Calculate average mass */
      for(i=0; i<NSPECIES; i++) {
	m_avg += nden[i]*mass[i]/Ntot;
      }
      
      f11 = Rho[k][iLon][iLat][NVERT-1]*exp(-m_avg*Hz);
      f12 = Rho[k][iLon][iLat+1][NVERT-1]*exp(-m_avg*Hz);
      f21 = Rho[k][iLon+1][iLat][NVERT-1]*exp(-m_avg*Hz);
      f22 = Rho[k][iLon+1][iLat+1][NVERT-1]*exp(-m_avg*Hz);
           
      **rho = (1.0/((x2-x1)*(y2-y1)))*(f11*(x2-xp)*(y2-yp)+f21*(xp-x1)*(y2-yp)+f12*(x2-xp)*(yp-y1)+f22*(xp-x1)*(yp-y1));
      **temp = temperature;

      }
    
  }

  if(DensityModel==1) {
    
    for(i=0; i<NSPECIES; i++) {
      
      /* Assume scale height also remains constant */
      nden_bilinear[i][0] = ndenvec[i][k][iLon][iLat][NVERT-1]*exp(-mass[i]*Hz);
      nden_bilinear[i][1] = ndenvec[i][k][iLon][iLat+1][NVERT-1]*exp(-mass[i]*Hz);
      nden_bilinear[i][2] = ndenvec[i][k][iLon+1][iLat][NVERT-1]*exp(-mass[i]*Hz);
      nden_bilinear[i][3] = ndenvec[i][k][iLon+1][iLat+1][NVERT-1]*exp(-mass[i]*Hz);
      
      nden[i] = BilinearInterp(LatDeg, LonDeg, lati, loni, nden_bilinear[i], iLat, iLon, it);
      
    }
    
    /* Calculate the total number density */
    for(i=0; i<NSPECIES; i++) {
      Ntot += nden[i];
    }
    
    /* Calculate average mass */
    for(i=0; i<NSPECIES; i++) {
      m_avg += nden[i]*mass[i]/Ntot;
    }
    
    rho_bilinear[0] = Rho[k][iLon][iLat][NVERT-1]*exp(-m_avg*Hz);
    rho_bilinear[1] = Rho[k][iLon][iLat+1][NVERT-1]*exp(-m_avg*Hz);
    rho_bilinear[2] = Rho[k][iLon+1][iLat][NVERT-1]*exp(-m_avg*Hz);
    rho_bilinear[3] = Rho[k][iLon+1][iLat+1][NVERT-1]*exp(-m_avg*Hz);
    
    **rho = BilinearInterp(LatDeg, LonDeg, lati, loni, rho_bilinear, iLat, iLon, it);
    
    **temp = temperature;
    
  }

}	    
