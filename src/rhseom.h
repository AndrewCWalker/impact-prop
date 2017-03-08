#ifndef RHSEOM_H
#define RHSEOM_H

typedef void (*rhseom_t)(double [6], double, double, int, double [6], double *, double (*)[3], double (*)[3], double (*)[3][3], double [MAXGDO][MAXGDO], double [MAXGDO][MAXGDO], struct densities_struct *, int, struct sat_struct *, int, struct msis_struct *, struct initCd_struct *, double **, int, int, struct RSM_struct *, int);

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
	    int nobs);

#endif
