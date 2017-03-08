#ifndef RK_H
#define RK_H

#include "rhseom.h"

void rk4(rhseom_t rhseom, 
	 double dt, double x[6], double t, int n, 
	 double nstate[6], double doubletime[], 
	 double (*resun)[3], double (*remoon)[3], double (*RITJ)[3][3], 
	 double (*CMAT)[MAXGDO], double (*SMAT)[MAXGDO], struct densities_struct *input,
	 //double etlist[MAXFILES], int Nfiles, double LatDeg[MAXLAT], double LonDeg[MAXLON], 
	 //double AltKm[MAXVERT], double Rho[2][MAXLON][MAXLAT][MAXVERT], 
	 //double ndenvec[MAXSPECIES][2][MAXLON][MAXLAT][MAXVERT], 
	 //double Temperature[2][MAXLON][MAXLAT][MAXVERT], 
	 int dnp, struct sat_struct *psatout, int it, struct msis_struct *pmsis,
	 struct initCd_struct *initCd, double **msismod, int isat, int atmcounter,
	 struct RSM_struct *RSMdata, int nobs);

void rk8(rhseom_t rhseom, 
	 double dt, double x[6], double t, int n, 
	 double nstate[6], double doubletime[], 
	 double (*resun)[3], double (*remoon)[3], double (*RITJ)[3][3], 
	 double (*CMAT)[MAXGDO], double (*SMAT)[MAXGDO], struct densities_struct *input,
	 //double etlist[MAXFILES], int Nfiles, double LatDeg[MAXLAT], double LonDeg[MAXLON], 
	 //double AltKm[MAXVERT], double Rho[2][MAXLON][MAXLAT][MAXVERT], 
	 //double ndenvec[MAXSPECIES][2][MAXLON][MAXLAT][MAXVERT], 
	 //double Temperature[2][MAXLON][MAXLAT][MAXVERT], 
	 int dnp, struct sat_struct *psatout, int it, struct msis_struct *pmsis,
	 struct initCd_struct *initCd, double **msismod, int isat, int atmcounter,
	 struct RSM_struct *RSMdata, int nobs);

#endif
