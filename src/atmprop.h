#ifndef ATMPROP_H
#define ATMPROP_H

#include "structs.h"

void atmospheric_properties(double r_ITRF[3], double et, struct densities_struct *input, double *rho, 
			    double *temp, double nden[6], struct msis_struct *pmsis, double **msismod, 
			    int isat, int it, int atmcounter, int nobs);

double densitycira72(double h);

double densityussa76(double h, double altitude[MAXVERT],
		     double density[MAXLON][MAXLAT][MAXVERT]);

void densitydynamicMSIS(double alt, double et, double geodetic_Lat, double geodetic_Long, double *rho, double *temp, double nden[6], struct msis_struct *pmsis, double **msismod, int isat, int nobs);

void dynamicMSISforHASDM(double alt, double et, double geodetic_Lat, double geodetic_Long, double *temp, double nden[6], struct msis_struct *pmsis, double **msismod, int isat, int nobs);

#endif
