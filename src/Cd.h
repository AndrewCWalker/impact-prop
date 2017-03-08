#ifndef CD_H
#define CD_H


double CdSetup(double U, double Ta, double n[NSPECIES], struct sat_struct *psatout, int it, 
	       struct initCd_struct *initCd, struct RSM_struct *RSMdata);

double compute_Cd_CLL(int Geometry, double U, double Ts, double Ta, double n, double sigmat, double alphan, double X[6], 
		      double theta, double phi, double R, double L, double H, double W);

double compute_Cd_diffuse(int Geometry, double U, double Ts, double Ta, double n, double alpha, double X[NSPECIES], 
			  double theta, double phi, double R, double L, double H, double W);



#endif
