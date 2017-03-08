#ifndef ACCELERATION_H
#define ACCELERATION_H

#include "structs.h"

void centralacc(double rvec[3], double cacc[3]);

void srpacc(double rvec[3], double et, double dt, double et_doubletime[], double (*r_e_sun_doubletime)[3], 
	    double sacc[3], int dnp, struct sat_struct *psatout, int it, struct initCd_struct *initCd);

double shadowfunc(double r_J2000[3], double rv_e3[3]);

void dragacc(double state[6], double r_ITRF[3], double v_ITRF[3], double et, double dt, int n, 
	     struct densities_struct *input, double dacc[3], struct sat_struct *psatout, int it, 
	     double RMT[3][3], struct msis_struct *pmsis, struct initCd_struct *initCd, double **msismod, 
	     int isat, int atmcounter, struct RSM_struct *RSMdata, int nobs);

void nonsphacc(double rvec[3],double nacc[3], double (*CMAT)[MAXGDO], double (*SMAT)[MAXGDO]);

void thirdbodyacc(double rvec[3], double et, double dt, double et_doubletime[], double (*r_e_sun_doubletime)[3], 
		  double (*r_e_moon_doubletime)[3], char thirdbodyname[], double tacc[3], int dnp);


#endif
