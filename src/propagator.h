#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include "structs.h"
#include <gsl/gsl_rng.h>

typedef void (*data_read)(char *, const char *, struct config_struct,
			  int, struct densities_struct *, double *, int *, 
			  int, double);

typedef void (*data_store)(char *, struct constants_struct *,
			   int, double *, struct sat_struct *);

void propagate(char *atm_dir, const char *data_dir,
	       struct config_struct config,
	       data_read get_data, data_store record_data,
	       double initSV[6], double DioramaDeltaT, double finalSV[6]);

void propagator_init(char *atm_dir, const char *data_dir, struct config_struct config, data_read get_data, 
		     data_store record_data, double DioramaDeltaT,
		     int *natm, int *n, int *dnp, int *numpoints, int *scn, int *nobs, 
		     char kernelpath[NKER][PATH_MAX], 
		     double **t, double **doubletime, double *deltat, double state0[6], double *psig, double *vsig, 
		     double ***state_cart, double ***state_elem, double ***initstate, double ***msismod, 
		     double (**RITJ)[3][3], double (**resun)[3], double (**remoon)[3], double (**CMAT)[MAXGDO], 
		     double (**SMAT)[MAXGDO], struct initCd_struct **initCd, struct RSM_struct **RSMdata, 
		     struct sat_struct **psatout, struct constants_struct **psatconstant, struct msis_struct **pmsis, 
		     struct densities_struct **input, gsl_rng **rr);

void propagator_core(char *atm_dir, const char *data_dir, struct config_struct config, data_read get_data, 
		     double initSV[6], double deltaT, double finalSV[6], double psig, double vsig, 
		     double **state_cart, double **state_elem, double state0[6], 
		     double **initstate, double *t, double deltat, double *doubletime,  
		     double (*RITJ)[3][3], double (*resun)[3], double (*remoon)[3], double (*CMAT)[MAXGDO], 
		     double (*SMAT)[MAXGDO], double **msismod, int numpoints, int *natm, int isat, int n, 
		     int dnp, int nobs, struct densities_struct *input, struct initCd_struct *initCd, 
		     struct sat_struct *psatout, struct msis_struct *pmsis, struct RSM_struct *RSMdata, 
		     gsl_rng *rr);

void propagator_write(char *atm_dir, const char *data_dir, data_store record_data, int numpoints, 
		      struct constants_struct *psatconstant, double *t, struct sat_struct *psatout,
		      int scn, int isat, double **state_cart);

void propagator_close(struct initCd_struct *initCd, struct RSM_struct *RSMdata, double *t, double *doubletime,
		      double (*CMAT)[MAXGDO], double (*SMAT)[MAXGDO], struct densities_struct *input, 
		      double (*RITJ)[3][3], double (*resun)[3], double (*remoon)[3], 
		      struct constants_struct *psatconstant, struct sat_struct *psatout,
		      double **state_cart, double **state_elem, double **initstate, int *natm, int dnp, int n, 
		      double **msismod, struct msis_struct *pmsis, gsl_rng *rr, char kernelpath[NKER][PATH_MAX]);


#endif
