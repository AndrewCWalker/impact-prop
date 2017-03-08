#ifndef INTERPOLATION_H
#define INTERPOLATION_H


double interpLatLong(double lati, double loni, double alt, double parameter[MAXLON][MAXLAT][MAXVERT], double LatDeg[MAXLAT], 
		     double LonDeg[MAXLON], double AltKm[MAXVERT], int it);

double BilinearInterp(double LatDeg[MAXLAT], double LonDeg[MAXLON], double lati, double loni, double param[], int iLat, int iLon, int it);

double interpVert(int iLat, int iLon, double alt, double AltKm[MAXVERT], double parameter[MAXLON][MAXLAT][MAXVERT], int it);

void extrapolate_atm_properties(int k, double lati, double loni, double alt, double LatDeg[MAXLAT], 
				double LonDeg[MAXLON], double AltKm[MAXVERT], double Rho[2][MAXLON][MAXLAT][MAXVERT], 
				double Temperature[2][MAXLON][MAXLAT][MAXVERT], 
				double ndenvec[NSPECIES][2][MAXLON][MAXLAT][MAXVERT], double **rho, double **temp, 
				double nden[6], double mass[NSPECIES], int it);


#endif
