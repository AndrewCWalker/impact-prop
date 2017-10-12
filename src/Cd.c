#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_bessel.h>
#include <time.h>
#include "defs.h"
#include "structs.h"
#include "Cd.h"
#include "emu.h"

extern int rkcounter;

double CdSetup(double U, double Ta, double n[NSPECIES], struct sat_struct *psatout, int it, 
	       struct initCd_struct *initCd, struct RSM_struct *RSMdata) {

  /* Geometries: 1 = sphere; 2 = cylinder; 3 = flat plate; 4 = cube; 5 = cuboid */

  /* Setup all the input properties to calculate the drag coefficient */
  /* #INPUTS: */
  /* U is the satellite velocity relative to the atmosphere in m/s */
  /* Ta is the atmospheric translational temperature in K */
  /* n is a list of the species number densities in the order [O, O2, N, N2, He, H] */
  /* theta is the orientation angle of the object in the x-y plane */
  /* phi is the orientation angle in the x-z plane */
  /* Sphere: phi and theta are irrelevant */
  /* Cylinder: theta = 0 represents cylinder axis perpendicular to velocity vector, U,  */
  /* while theta = Pi/2 represents the cylinder axis */
  /* parallel to the velocity vector. Phi doesn't matter by symmetry. */
  /* Flat Plate: theta = 0, phi = 0 represents the surface normal of the plate parallel to the velocity vector, U  */
  /* (e.g. broadside) */
  /* Cube: theta = 0, phi = 0 represents the surface normal of the front of the cube parallel  */
  /* to the velocity vector, U */
  /* Cuboid: theta = 0, phi = 0 represents the LxH face surface normal parallel to the velocity vector, U.  */
  /* At, theta = Pi/2 and phi = 0 */
  /* the WxH face's surface normal is parallel to the velocity vector, and at theta = 0, phi = Pi/2,  */
  /* the LxW face is normal to the velocity vector */
  /* R is the radius [m] and is only relevant to the cylinder, currently. */
  /* L is the length [m] and is only relevant to the cylinder and the cuboid, currently. */
  /* H is the height [m] and is only relevant to the cuboid, currently. */
  /* W is the width [m] and is only relevant to the cuboid, currently. */

  int i;
  int Geometry = initCd->geometry;

  /* Define geometric parameters */
  double theta = initCd->yaw;
  double phi = initCd->pitch;
  double R = initCd->radius;
  double L = initCd->length;
  double H = initCd->height;
  double W = initCd->width;

  double mass[NSPECIES], X[NSPECIES], Xemu[NSPECIES];
  double ads_input[7], srf_input[7];
  double n_TOT = 0.0;
  double m_avg = 0.0;

  /* Define mass of satellite surface */
  double m_surf = initCd->surf_mass;

  double Ts = 300.0;    /* Constant satellite surface temperature of 300K */
  double sigmat = 1.0;  /* Constant tangential momentum accomodation of 1.0 based on Comsa (1980)/Porodnov (1974) */
  double kB = 1.38e-23; /* Boltzmann's constant */

  double kL_CLL = 2.89e6;    /* CLL - Define Langmuir Adsorpate Constant */
  double aF_CLL = 0.089;    /* CLL - Define Freundlich Alpha Constant */
  double kF_CLL = 2.35;      /* CLL - Define Freundlich K constant */

  double kL_DRIA = 1.44e6;   /* DRIA - Define Langmuir Adsorpate Constant */
  double aF_DRIA = 0.17;     /* DRIA - Define Freundlich Alpha Constant */
  double kF_DRIA = 5.7 ;     /* DRIA - Define Freundlich K constant */

  /* alpha for atomic oxygen covered surface is unity */
  double alpha_ads = 1.0;
  double alphan_ads = 1.0;

  double mu;
  double fsc = 0.0;
  double alpha_sub, alphan_sub;
  double Cd_ads, Cd_srf;
  double Po;
  double Cd_TOTAL;

  struct sat_struct *satout;
  satout = psatout + it;

  U *= 1.0e3;   /* Convert speed to m/s */

  /* Define atomic/molecular masses of each species */
  mass[0] = 2.657e-26;  /* atomic oxygen */
  mass[1] = 5.314e-26;  /* diatomic oxygen */
  mass[2] = 2.326e-26;  /* atomic nitrogen */
  mass[3] = 4.652e-26;  /* diatomic nitrogen */
  mass[4] = 6.646e-27;  /* helium */
  mass[5] = 1.674e-27;  /* hydrogen */

  /* Compute total number density */
  for(i=0; i<NSPECIES; i++) {
    n_TOT += n[i];
    X[i] = 0.0;
  }
  
  /* Compute mole fractions for each species */
  for(i=0; i<NSPECIES; i++) {
    X[i] = n[i]/n_TOT;
    satout->X[i] = X[i];
  }

  /* Swap order of species for response surface emulator */
  Xemu[0] = X[5];
  Xemu[1] = X[4];
  Xemu[2] = X[2];
  Xemu[3] = X[3];
  Xemu[4] = X[0];
  Xemu[5] = X[1];
  
  /* Compute the average mass */
  for(i=0; i<NSPECIES; i++) {
    m_avg += X[i]*mass[i];
  }
    
  /* alpha for clean surface is based on Goodman's (1966) empirical formula */
  mu = m_avg/m_surf;
  //alpha_sub = 3.6*mu/pow(1+mu, 2.0);
  alpha_sub = 2.4*mu/pow(1+mu, 2.0);

  alphan_sub = 2.0*alpha_sub-1.0;
  if(alphan_sub < 0.0) {
    alphan_sub = 0.0;
  }

  if(initCd->gsi_model) { /* CLL */

    /* Define response surface emulator inputs */
    ads_input[0] = srf_input[0] = U;
    ads_input[1] = srf_input[1] = Ts;
    ads_input[2] = srf_input[2] = Ta;
    ads_input[4] = srf_input[4] = sigmat;
    ads_input[5] = srf_input[5] = theta;
    ads_input[6] = srf_input[6] = phi;

    ads_input[3] = alphan_ads;
    srf_input[3] = alphan_sub;

    //printf("alphan_ads = %e alphan_sub = %e\n", alphan_ads, alphan_sub);

    if(RSM) {
      emu(RSMdata, ads_input, Xemu, &Cd_ads);
      emu(RSMdata, srf_input, Xemu, &Cd_srf);
    } else {
      Cd_ads = compute_Cd_CLL(Geometry, U, Ts, Ta, n_TOT, sigmat, alphan_ads, X, theta, phi, R, L, H, W);
      Cd_srf = compute_Cd_CLL(Geometry, U, Ts, Ta, n_TOT, sigmat, alphan_sub, X, theta, phi, R, L, H, W);
    }
  } else { /* Diffuse */
    Cd_ads = compute_Cd_diffuse(Geometry, U, Ts, Ta, n_TOT, alpha_ads, X, theta, phi, R, L, H, W);
    Cd_srf = compute_Cd_diffuse(Geometry, U, Ts, Ta, n_TOT, alpha_sub, X, theta, phi, R, L, H, W) ; 
  }     

  /* Compute the partial pressure of atomic oxygen */
  Po = n[0]*kB*Ta;

  /* Compute fsc: fractional coverage of surface by atomic oxygen */
  if(initCd->gsi_model) {                                     /* CLL */
    if(initCd->ads_model) {                                   /* Freundlich */
      fsc = kF_CLL*pow(Po, aF_CLL);
      if(fsc > 1.0) {
	fsc = 1.0;
      }
    } else {                                          /* Langmuir */
      fsc = (kL_CLL*Po)/(1+kL_CLL*Po);
    }
  } else {                                            /* Diffuse Reflection with Incomplete Accommodation (DRIA) */
    if(initCd->ads_model) {                                   /* Freundlich */
      fsc = kF_DRIA*pow(Po, aF_DRIA);
      if(fsc > 1.0) {
	fsc = 1.0;
      } else {                                        /* Langmuir */
	fsc = (kL_DRIA*Po)/(1+kL_DRIA*Po);
      }
    }
  }

  Cd_TOTAL = fsc*Cd_ads + (1.0-fsc)*Cd_srf;

  //if(it>22 && it<28) { 
  //  printf("it = %d Cd_TOTAL = %e\n", it, Cd_TOTAL);
  //  printf("Cd_TOTAL = %e X[0] = %e X[1] = %e X[2] = %e X[3] = %e X[4] = %e X[5] = %e\n", Cd_TOTAL, X[0], X[1], X[2], X[3], X[4], X[5]);
  //  printf("U = %e Ts = %e Ta = %e sigmat = %e alphan_ads = %e theta = %e phi = %e\n", U, Ts, Ta, sigmat, alphan_ads, theta, phi);
  //}

  if(rkcounter % RKN ==0) {
    satout->alpha_n = fsc*alphan_ads + (1.0-fsc)*alphan_sub;
    satout->theta = theta;
    satout->phi = phi;
  }

  return(Cd_TOTAL);

}


double compute_Cd_CLL(int Geometry, double U, double Ts, double Ta, double n, double sigmat, double alphan, double X[6], 
		      double theta, double phi, double R, double L, double H, double W) {

  /* Compute the drag coefficient based on satellite and atmospheric properties */
  int i;

  double mass[NSPECIES];
  double m_avg = 0.0;
  
  double alpha[3][NSPECIES];
  double beta[3][NSPECIES];
  double gamma[3][NSPECIES];
  double delta[3][NSPECIES];
  double sigman, sigman2;
  double Cd[NSPECIES];
  double Vmp, s;
  double s2, s3, s4;
  double Cd_TOTAL = 0.0;

  double cp = cos(phi*RAD2DEG);
  double ct = cos(theta*RAD2DEG);
  double sp = sin(phi*RAD2DEG);
  double st = sin(theta*RAD2DEG);

  double sqrtPI = sqrt(M_PI);

  if( n > 1.0e16) {
    printf("Atmosphere is no longer free molecular - Cd may be wrong\n");
  }

  /* Define masses of each chemical species */
  mass[0] = 2.657e-26; /* atomic oxygen */
  mass[1] = 5.314e-26; /* diatomic oxygen */
  mass[2] = 2.326e-26; /* atomic nitrogen */
  mass[3] = 4.652e-26; /* diatomic nitrogen */
  mass[4] = 6.646e-27; /* atomic helium */
  mass[5] = 1.674e-27; /* atomic hydrogen */

  /* Compute the average mass */
  for(i=0; i<NSPECIES; i++) {
    m_avg += X[i]*mass[i];
  }

  /* Define least-squares coefficients for curve fit of sigman2 */   
  /* ######################## Flat Plate ########################## */

  /* Atomic Oxygen */
  alpha[0][0] = 5.85;
  beta[0][0] = 0.20;
  gamma[0][0] = 0.48;
  delta[0][0] = 31.0;

  /* Diatomic Oxygen */
  alpha[0][1] = 6.3;
  beta[0][1] = 0.26;
  gamma[0][1] = 0.42;
  delta[0][1] = 20.5;

  /* Atomic Nitrogen */
  alpha[0][2] = 4.9;
  beta[0][2] = 0.32;
  gamma[0][2] = 0.42;
  delta[0][2] = 8.0;

  /* Diatomic Nitrogen */
  alpha[0][3] = 6.6;
  beta[0][3] = 0.22;
  gamma[0][3] = 0.48;
  delta[0][3] = 35.0;

  /* Helium */
  alpha[0][4] = 4.5;
  beta[0][4] = 0.38;
  gamma[0][4] = 0.51;
  delta[0][4] = 5.8;

  /* Atomic Hydrogen */
  alpha[0][5] = 3.6;
  beta[0][5] = 0.48;
  gamma[0][5] = 0.52;
  delta[0][5] = 2.8;


  /* ######################## Sphere ########################## */

  /* Atomic Oxygen */
  alpha[1][0] = 4.4;
  beta[1][0] = 0.32;
  gamma[1][0] = 0.48;
  delta[1][0] = 11.0;
  
  /* Diatomic Oxygen */
  alpha[1][1] = 5.45;
  beta[1][1] = 0.18;
  gamma[1][1] = 0.50;
  delta[1][1] = 49.0;

  /* Atomic Nitrogen */
  alpha[1][2] = 4.75;
  beta[1][2] = 0.24;
  gamma[1][2] = 0.50;
  delta[1][2] = 20.0;

  /* Diatomic Nitrogen */
  alpha[1][3] = 5.5;
  beta[1][3] = 0.18;
  gamma[1][3] = 0.50;
  delta[1][3] = 51.0;

  /* Helium */
  alpha[1][4] = 4.15;
  beta[1][4] = 0.35;
  gamma[1][4] = 0.52;
  delta[1][4] = 8.0;

  /* Atomic Hydrogen */
  alpha[1][5] = 3.4;
  beta[1][5] = 0.54;
  gamma[1][5] = 0.54;
  delta[1][5] = 2.8;

  /* ######################## Cylinder ########################## */

  /* Atomic Oxygen */
  alpha[2][0] = 5.35;
  beta[2][0] = 0.36;
  gamma[2][0] = 0.56;
  delta[2][0] = 11.0;

  /* Diatomic Oxygen */
  alpha[2][1] = 6.0;
  beta[2][1] = 0.34;
  gamma[2][1] = 0.58;
  delta[2][1] = 15.0;

  /* Atomic Nitrogen */
  alpha[2][2] = 5.15;
  beta[2][2] = 0.30;
  gamma[2][2] = 0.48;
  delta[2][2] = 13.0;

  /* Diatomic Nitrogen */
  alpha[2][3] = 6.3;
  beta[2][3] = 0.24;
  gamma[2][3] = 0.54;
  delta[2][3] = 33.0;

  /* Helium */
  alpha[2][4] = 4.55;
  beta[2][4] = 0.34;
  gamma[2][4] = 0.48;
  delta[2][4] = 8.0;

  /* Atomic Hydrogen */
  alpha[2][5] = 3.7;
  beta[2][5] = 0.44;
  gamma[2][5] = 0.56;
  delta[2][5] = 3.7;

  sigman = 1.0 - sqrt(1.0 - alphan);

  /* Initialize Species Specific Cd */
  for(i=0; i<NSPECIES; i++) {
    Cd[i] = 0.0;
  }

  /* Loop over each species */
  for(i=0; i<NSPECIES; i++) {

    /* Calculate speed ratio, s, and reflected kinetic temperature, Tkr */
    Vmp = sqrt(2.0*1.38e-23*Ta/mass[i]);
    s = U/Vmp;
    s2 = pow(s, 2);
    s3 = pow(s, 3);
    s4 = pow(s, 4);
    
    /*************************************************************************************/
    /***************************************SPHERE****************************************/
    /*************************************************************************************/
    if(Geometry==1) {
      
      sigman2 = exp(-alpha[1][i]*pow(1.0-alphan,beta[1][i]))*pow(Ts/Ta, gamma[1][i])*delta[1][i]/s;
        
      if(alphan > 0.99) {
	sigman2 = 1.0 - pow(1.0 - alphan, 0.04);
      }

      Cd[i] = ((2.0-sigman+sigmat)/2.0)*(((2.0*s2+1.0)/(sqrtPI*s3))*exp(-s2)+
					 ((4.0*s4+4.0*s2-1.0)/(2.0*s4))*erf(s))+
	                                 (2.0*sigman2*sqrtPI/(3.0*s))*sqrt(Ts/Ta);
     
    }

    /*************************************************************************************/
    /***************************************CYLINDER**************************************/
    /*************************************************************************************/
    if(Geometry==2) {

      sigman2 = exp(-alpha[2][i]*pow(1.0-alphan,beta[2][i]))*pow(Ts/Ta, gamma[2][i])*delta[2][i]/s;
      if(alphan > 0.99) {
	sigman2 = 1.0 - pow(1.0 - alphan, 0.04);
      }
         
      Cd[i] = ((
	    /* Flat Plate Axial */
	    ((M_PI*R*R)/(2.0*R*L))*sigmat*(2.0*ct*st*erf(s*st)+
	    (2.0*ct)/(sqrtPI*s)*exp(-s2*pow(st,2.0)))+

	    /* Cylinder perpendicular to flow */
	    (((sqrtPI/3.0)*(4.0 - 2.0*sigman + sigmat)*ct*exp(-(s2*pow(ct,2.0))/2.0))/s)*
	    ((s2*pow(ct,2.0)+3.0/2.0)*gsl_sf_bessel_I0((s2*pow(ct,2.0))/2.0)+
	     ((s2*pow(ct,2.0))+1.0/2.0)*gsl_sf_bessel_I1((s2*pow(ct,2.0))/2.0))+
	    (pow(M_PI, 3.0/2.0)/(4.0*s))*sqrt(Ts/Ta)*sigman2*ct

	    /* Multiply Sum by cos(theta) */
	    )*ct+

	    /* Flat Plate Normal */
	    (((M_PI*R*R)/(2.0*R*L))*((2.0*(2.0-sigman)*st/(sqrtPI*s))*
	    exp(-s2*pow(st, 2.0))+(2.0-sigman)*(2.0*pow(st, 2.0)+
	    1.0/s2)*erf(s*st)+(sqrtPI/s)*sqrt(Ts/Ta)*sigman2*st)+

	    /* Cylinder parallel to flow */
	    ((sqrtPI*sigmat*st)/s)*exp(-s2*pow(ct,2.0)/2.0)*
	    ((1.0+s2*pow(ct,2.0))*
	     gsl_sf_bessel_I0(s2*pow(ct,2.0)/2.0)+s2*pow(ct,2.0)*
	    gsl_sf_bessel_I1((s2*pow(ct,2.0)/2.0))))

	    /* Multiply Sum by sin(theta) */
	    *st);

    }

    /***************************************************************************************/
    /***************************************FLAT PLATE**************************************/
    /***************************************************************************************/
    if(Geometry==3) {
    
      sigman2 = exp(-alpha[0][i]*pow(1.0-alphan,beta[0][i]))*pow(Ts/Ta, gamma[0][i])*delta[0][i]/s;
      if(alphan > 0.99) {
	sigman2 = 1.0 - pow(1.0 - alphan, 0.04);
      }

      Cd[i] = ((((2.0*(2.0-sigman)*ct*cp)/(sqrtPI*s))*
	       exp(-s2*pow(ct*cp, 2.0))+
	       (2.0-sigman)*(2.0*pow(ct*cp, 2.0)+1.0/s2)*erf(s*ct*cp)+
	       (sqrtPI*ct*cp/s)*sqrt(Ts/Ta)*sigman2)*ct*cp+                  

	       sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(ct*cp,2.0))+
	       2.0*ct*cp*erf(s*ct*cp))*pow((st*cp),2.0)+

	       sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(ct*cp,2.0))+
	       2.0*ct*cp*erf(s*ct*cp))*pow(sp,2.0));

    }

    /***************************************************************************************/
    /********************************************CUBE***************************************/
    /***************************************************************************************/


    if (Geometry==4) {
     
      sigman2 = exp(-alpha[0][i]*pow(1.0-alphan,beta[0][i]))*pow(Ts/Ta, gamma[0][i])*delta[0][i]/s;
      if(alphan > 0.99) {
	sigman2 = 1.0 - pow(1.0 - alphan, 0.04);
      }

      Cd[i] = (
	       /* Front Surface (Z-direction */
	       (((2.0*(2.0-sigman)*ct*cp)/(sqrtPI*s))*
		exp(-s2*pow(ct*cp, 2.0))+
		(2.0-sigman)*(2.0*pow(ct*cp, 2.0)+1.0/s2)*erf(s*ct*cp)+
		(sqrtPI*ct*cp/s)*sqrt(Ts/Ta)*sigman2)*ct*cp+
	       
	       /* Side Surface (Y-direction) */
	       (((2.0*(2.0-sigman)*st*cp)/(sqrtPI*s))*
		exp(-s2*pow(st*cp, 2.0))+
		(2.0-sigman)*(2.0*pow(st*cp, 2.0)+1.0/s2)*erf(s*st*cp)+
		(sqrtPI*st*cp/s)*sqrt(Ts/Ta)*sigman2)*st*cp+

	       /* Top Surface (Z-direction) */
	       (((2.0*(2.0-sigman)*sp)/(sqrtPI*s))*exp(-s2*pow(sp, 2.0))+
	       (2.0-sigman)*(2.0*pow(sp, 2.0)+1.0/s2)*erf(s*sp)+
	       (sqrtPI*sp/s)*sqrt(Ts/Ta)*sigman2)*sp+  

	       /* Y-direction Shear on X-face */
	       sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(ct*
	       cp,2.0))+2*ct*cp*erf(s*ct*
	       cp))*pow((st*cp),2.0)+

	       /* Z-direction Shear on X-face */
	       sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(ct*cp,2.0))+
	       2*ct*cp*erf(s*ct*cp))*pow(sp,2.0)+

	       /* X-direction Shear on Y-face */
	       sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(st*cp,2.0))+
	       2*st*cp*erf(s*st*cp))*pow((ct*cp),2.0)+

	       /* Z-direction Shear on Y-face */
	       sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(st*cp,2.0))+
	       2*st*cp*erf(s*st*cp))*pow(sp,2.0)+

	       /* X-direction Shear on Z-face */
	       sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(sp,2.0))+
	       2*sp*erf(s*sp))*pow((ct*cp),2.0)+

	       /* Y-direction Shear on Z-face */
	       sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(sp,2.0))+
	       2*sp*erf(s*sp))*pow((st*cp),2.0));

    }

    /***************************************************************************************/
    /********************************************CUBOID*************************************/
    /***************************************************************************************/

    if (Geometry==5) {
    
      sigman2 = exp(-alpha[0][i]*pow(1.0-alphan,beta[0][i]))*pow(Ts/Ta, gamma[0][i])*delta[0][i]/s;
      if(alphan > 0.99) {
	sigman2 = 1.0 - pow(1.0 - alphan, 0.04);
      }

      Cd[i] = (
	       /* Front Surface (Z-direction) */
	       (((2.0*(2.0-sigman)*ct*cp)/(sqrtPI*s))*
	       exp(-s2*pow(ct*cp, 2.0))+
	       (2.0-sigman)*(2.0*pow(ct*cp, 2.0)+1.0/s2)*erf(s*ct*cp)+
	       (sqrtPI*ct*cp/s)*sqrt(Ts/Ta)*sigman2)*ct*cp+

	       /* Side Surface (Y-direction) */
	       ((H*W)/(L*H))*(((2.0*(2.0-sigman)*st*cp)/(sqrtPI*s))*
	       exp(-s2*pow(st*cp, 2.0))+
	       (2.0-sigman)*(2.0*pow(st*cp, 2.0)+1.0/s2)*erf(s*st*cp)+
	       (sqrtPI*st*cp/s)*sqrt(Ts/Ta)*sigman2)*st*cp+

	       /* Top Surface (Z-direction) */
	       ((L*W)/(L*H))*(((2.0*(2.0-sigman)*sp)/(sqrtPI*s))*exp(-s2*pow(sp, 2.0))+
               (2.0-sigman)*(2.0*pow(sp, 2.0)+1.0/s2)*erf(s*sp)+
	       (sqrtPI*sp/s)*sqrt(Ts/Ta)*sigman2)*sp+  

	       /* Y-direction Shear on X-face */
	       sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(ct*cp,2.0))+
	       2*ct*cp*erf(s*ct*cp))*pow((st*cp),2.0)+

	       /* Z-direction Shear on X-face */
	       sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(ct*cp,2.0))+
	       2*ct*cp*erf(s*ct*cp))*pow(sp,2.0)+

	       /* X-direction Shear on Y-face */
	       ((H*W)/(L*H))*sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(st*
	       cp,2.0))+2*st*cp*erf(s*st*
	       cp))*pow((ct*cp),2.0)+

	       /* Z-direction Shear on Y-face */
	       ((H*W)/(L*H))*sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(st*
	       cp,2.0))+2*st*cp*erf(s*st*cp))*pow(sp,2.0)+

	       /* X-direction Shear on Z-face */
	       ((L*W)/(L*H))*sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(sp,2.0))+
	       2*sp*erf(s*sp))*pow((ct*cp),2.0)+

	       /* Y-direction Shear on Z-face */
	       ((L*W)/(L*H))*sigmat*(2.0/(sqrtPI*s)*exp(-s2*pow(sp,2.0))+
	       2*sp*erf(s*sp))*pow((st*cp),2.0));

    }
  }

  

  /* Perform weighted sum over all species */
  for(i=0; i<NSPECIES; i++) {
    Cd_TOTAL += X[i]*mass[i]*Cd[i];
  }

  Cd_TOTAL /= m_avg;

  return(Cd_TOTAL);

}
    

double compute_Cd_diffuse(int Geometry, double U, double Ts, double Ta, double n, double alpha, double X[NSPECIES], 
			  double theta, double phi, double R, double L, double H, double W) {

  int i;
  
  double mass[NSPECIES];
  double m_avg = 0.0;

  double kB = 1.38e-23; /* Define Boltzmann Constant */
  double Cd[NSPECIES];
  double s, Vmp, Tkr;
  double s2, s3, s4;

  double Cd_TOTAL = 0.0;

  double cp = cos(phi*RAD2DEG);
  double ct = cos(theta*RAD2DEG);
  double sp = sin(phi*RAD2DEG);
  double st = sin(theta*RAD2DEG);

  double sqrtPI = sqrt(M_PI);

  /* Check that atmospheric number density, n, is in the free molecular limit */
  if(n>1.0e16) {
    printf("Atmosphere is no longer free molecular - Cd may be wrong\n");
  }

  /* Define masses of each chemical species */
  mass[0] = 2.657e-26; /* atomic oxygen */
  mass[1] = 5.314e-26; /* diatomic oxygen */
  mass[2] = 2.326e-26; /* atomic nitrogen */
  mass[3] = 4.652e-26; /* diatomic nitrogen */
  mass[4] = 6.646e-27; /* atomic helium */
  mass[5] = 1.674e-27; /* atomic hydrogen */

  /* Compute the average mass */
  for(i=0; i<NSPECIES; i++) {
    m_avg += X[i]*mass[i];
    Cd[i] = 0.0;
  }

  /* Loop over each species */
  for(i=0; i<NSPECIES; i++) {
    /* Calculate speed ratio, s, and reflected kinetic temperature, Tkr */
    Vmp = sqrt(2.0*1.38e-23*Ta/mass[i]);
    s = U/Vmp;
    s2 = pow(s, 2);
    s3 = pow(s, 3);
    s4 = pow(s, 4);
    Tkr = (mass[i]*pow(U,2)/(3.0*kB))*(1-alpha)+alpha*Ts;

    /***************************************************************************************/
    /********************************************SPHERE*************************************/
    /***************************************************************************************/
   
    if(Geometry==1) {
    
      Cd[i] = (((2.0*s2+1.0)/(sqrtPI*s3))*exp(-s2)+
	       ((4.0*s4+4.0*s2-1.0)/(2.0*s4))*erf(s))+(2.0*sqrtPI/(3.0*s))*sqrt(Tkr/Ta);
    }

    /***************************************************************************************/
    /******************************************CYLINDER*************************************/
    /***************************************************************************************/


    if (Geometry==2) {
     
      Cd[i] = ((
		/* Flat Plate Axial */
		((M_PI*R*R)/(2.0*R*L))*(2.0*ct*st*erf(s*st)+
		(2.0*ct)/(sqrtPI*s)*exp(-s2*pow(st,2.0)))+
                                
		/* Cylinder perpendicular to flow */
		(((sqrtPI/3.0)*(3.0)*ct*exp(-(s2*pow(ct,2.0))/2.0))/s)*
		((s2*pow(ct,2.0)+3.0/2.0)*gsl_sf_bessel_I0((s2*pow(ct,2.0))/2.0)+
		((s2*pow(ct,2.0))+1.0/2.0)*gsl_sf_bessel_I1((s2*pow(ct,2.0))/2.0))+
		(pow(M_PI, 3.0/2.0)/(4.0*s))*sqrt(Tkr/Ta)*ct
                                
		/* Multiply Sum by cos(theta) */
		)*ct+
                                
	       /* Flat Plate Normal */
	       (((M_PI*R*R)/(2.0*R*L))*((2.0*st/(sqrtPI*s))*exp(-s2*pow(st, 2.0))+
	       (2.0*pow(st, 2.0)+1.0/s2)*erf(s*st)+(sqrtPI/s)*sqrt(Tkr/Ta)*st)+
                                
		/* Cylinder parallel to flow */
		((sqrtPI*st)/s)*exp(-s2*pow(ct,2.0)/2.0)*
		((1.0+s2*pow(ct,2.0))*gsl_sf_bessel_I0(s2*pow(ct,2.0)/2.0)+
		 s2*pow(ct,2.0)*gsl_sf_bessel_I1((s2*pow(ct,2.0)/2.0))))*st);

    }

    /***************************************************************************************/
    /****************************************FLAT PLATE*************************************/
    /***************************************************************************************/


    if (Geometry==3) {
    
      Cd[i] = ((((2.0*ct*cp)/(sqrtPI*s))*exp(-s2*pow(ct*cp, 2.0))+
		(2.0*pow(ct*cp, 2.0)+1.0/s2)*erf(s*ct*cp)+
		(sqrtPI*ct*cp/s)*sqrt(Tkr/Ta))*ct*cp+                         

	       (2.0/(sqrtPI*s)*exp(-s2*pow(ct*cp,2.0))+
		2.0*ct*cp*erf(s*ct*cp))*pow((st*cp),2.0)+

	       (2.0/(sqrtPI*s)*exp(-s2*pow(ct*cp,2.0))+
		2.0*ct*cp*erf(s*ct*cp))*pow(sp,2.0));

    }


    /***************************************************************************************/
    /*******************************************CUBE****************************************/
    /***************************************************************************************/

    if (Geometry==4) {
   
      Cd[i] = (
	       /* Front Surface (Z-direction */
	       (((2.0*ct*cp)/(sqrtPI*s))*exp(-s2*pow(ct*cp, 2.0))+
	       (2.0*pow(ct*cp, 2.0)+1.0/s2)*erf(s*ct*cp)+
	       (sqrtPI*ct*cp/s)*sqrt(Tkr/Ta))*ct*cp+

	       /* Side Surface (Y-direction) */
	       (((2.0*st*cp)/(sqrtPI*s))*exp(-s2*pow(st*cp, 2.0))+
	       (2.0*pow(st*cp, 2.0)+1.0/s2)*erf(s*st*cp)+
	       (sqrtPI*st*cp/s)*sqrt(Tkr/Ta))*st*cp+

	       /* Top Surface (Z-direction) */
	       (((2.0*sp)/(sqrtPI*s))*exp(-s2*pow(sp, 2.0))+
		(2.0*pow(sp, 2.0)+1.0/s2)*erf(s*sp)+
		(sqrtPI*sp/s)*sqrt(Tkr/Ta))*sp+  

	       /* Y-direction Shear on X-face */
	       (2.0/(sqrtPI*s)*exp(-s2*pow(ct*cp,2.0))+
		2*ct*cp*erf(s*ct*cp))*pow((st*cp),2.0)+

	       /* Z-direction Shear on X-face */
	       (2.0/(sqrtPI*s)*exp(-s2*pow(ct*cp,2.0))+
		2*ct*cp*erf(s*ct*cp))*pow(sp,2.0)+

	       /* X-direction Shear on Y-face */
	       (2.0/(sqrtPI*s)*exp(-s2*pow(st*cp,2.0))+
		2*st*cp*erf(s*st*cp))*pow((ct*cp),2.0)+

	       /* Z-direction Shear on Y-face */
	       (2.0/(sqrtPI*s)*exp(-s2*pow(st*cp,2.0))+
		2*st*cp*erf(s*st*cp))*pow(sp,2.0)+

	       /* X-direction Shear on Z-face */
	       (2.0/(sqrtPI*s)*exp(-s2*pow(sp,2.0))+
		2*sp*erf(s*sp))*pow((ct*cp),2.0)+

	       /* Y-direction Shear on Z-face */
	       (2.0/(sqrtPI*s)*exp(-s2*pow(sp,2.0))+
		2*sp*erf(s*sp))*pow((st*cp),2.0));

    }


    /***************************************************************************************/
    /*******************************************CUBE****************************************/
    /***************************************************************************************/

    if (Geometry==5) {

      Cd[i] = (
	       /* Front Surface (Z-direction */
	       (((2.0*ct*cp)/(sqrtPI*s))*exp(-s2*pow(ct*cp, 2.0))+
	       (2.0*pow(ct*cp, 2.0)+1.0/s2)*erf(s*ct*cp)+
	       (sqrtPI*ct*cp/s)*sqrt(Tkr/Ta))*ct*cp+

	       /* Side Surface (Y-direction) */
	       ((H*W)/(L*H))*(((2.0*st*cp)/(sqrtPI*s))*
	       exp(-s2*pow(st*cp, 2.0))+
	       (2.0*pow(st*cp, 2.0)+1.0/s2)*erf(s*st*cp)+
	       (sqrtPI*st*cp/s)*sqrt(Tkr/Ta))*st*cp+

	       /* Top Surface (Z-direction) */
	       ((L*W)/(L*H))*(((2.0*sp)/(sqrtPI*s))*exp(-s2*pow(sp, 2.0))+
	       (2.0*pow(sp, 2.0)+1.0/s2)*erf(s*sp)+
	       (sqrtPI*sp/s)*sqrt(Tkr/Ta))*sp+  

	       /* Y-direction Shear on X-face */
	       (2.0/(sqrtPI*s)*exp(-s2*pow(ct*cp,2.0))+
	       2*ct*cp*erf(s*ct*cp))*pow((st*cp),2.0)+

	       /* Z-direction Shear on X-face */
	       (2.0/(sqrtPI*s)*exp(-s2*pow(ct*cp,2.0))+
	       2*ct*cp*erf(s*ct*cp))*pow(sp,2.0)+

	       /* X-direction Shear on Y-face */
	       ((H*W)/(L*H))*(2.0/(sqrtPI*s)*exp(-s2*pow(st*cp,2.0))+
	       2*st*cp*erf(s*st*cp))*pow((ct*cp),2.0)+

	       /* Z-direction Shear on Y-face */
	       ((H*W)/(L*H))*(2.0/(sqrtPI*s)*exp(-s2*pow(st*cp,2.0))+
	       2*st*cp*erf(s*st*cp))*pow(sp,2.0)+

	       /* X-direction Shear on Z-face */
	       ((L*W)/(L*H))*(2.0/(sqrtPI*s)*exp(-s2*pow(sp,2.0))+
	       2*sp*erf(s*sp))*pow((ct*cp),2.0)+

	       /* Y-direction Shear on Z-face */
	       ((L*W)/(L*H))*(2.0/(sqrtPI*s)*exp(-s2*pow(sp,2.0))+
	       2*sp*erf(s*sp))*pow((st*cp),2.0));

    }

  }

  /* Perform weighted sum over all species */
  for(i=0; i<NSPECIES; i++) {
    Cd_TOTAL += X[i]*mass[i]*Cd[i];
  }

  /* Find total Cd for specifici geometry */
  Cd_TOTAL /= m_avg;

  return(Cd_TOTAL);

}
