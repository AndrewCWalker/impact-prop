#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <errno.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "defs.h"
#include "propagator.h"
#if !DIORAMA
#include "SpiceUsr.h"
#endif
#include "misc.h"
#include "transformation.h"
#include "rk.h"
#include "rhseom.h"
#include "io.h"
#include "msisinputs.h"
#if PARALLEL
#include "mpi.h"
#endif /* PARALLEL */

#if PARALLEL
int rank, nProcs;
MPI_Comm io_comm;
#endif /* PARALLEL */


void propagate(char *atm_dir, const char *data_dir, struct config_struct config, data_read get_data, 
	       data_store record_data, double initSV[6], double DioramaDeltaT, double finalSV[6]) {

  gsl_rng *rr = NULL;                /* Pointer for GSL random number generator */

  char kernelpath[NKER][PATH_MAX];   /* Path for the cspice kernels        [Units: unitless] */

  int *natm = NULL;   /* Array of atmospheric ensemble index                       [Units: unitless] */
  int nobs = 0;       /* Number of observations for MSIS data array                [Units: unitless] */
  int n = 0;          /* Number of satellites to simulate                          [Units: unitless] */
  int dnp = 0;        /* Number of points for half timestep calculation            [Units: unitless] */
  int numpoints = 0;  /* Number of points for full timestep calculation            [Units: unitless] */
  int scn = 0;        /* Satellite catalog number                                  [Units: unitless] */
  int isat = 0;       /* Satellite index number for loop                           [Units: unitless] */
  int satloopmin = 0; /* Starting point for loop over satellites                   [Units: unitless] */
  int satloopmax = 0; /* Ending point for loop over satellites                     [Units: unitless] */

  double satperproc = 0.0;       /* Satellites computeed per processor      [Units: unitless] */
  double psig = 0.0;             /* Position standard deviation             [Units: km]       */
  double vsig = 0.0;             /* Velocity standard deviation             [Units: km/s]     */
  double deltat = 0.0;           /* Timestep of simulation time             [Units: seconds]  */
  double *t = NULL;              /* Ephemeris time array (since J2000)      [Units: seconds]  */
  double *doubletime = NULL;     /* Ephemeris half time array (since J2000) [Units: seconds]  */
  double **state_elem = NULL;    /* State vector in keplerian elements      [Units: {km, unitless, rad, rad, rad, rad}] */
  double **state_cart = NULL;    /* State vector in cartesian elements      [Units: {km, km, km, km/s, km/s, km/s}] */
  double **initstate = NULL;     /* Initial J2000 state vector for each sat [Units: {km, km, km, km/s, km/s, km/s}] */
  double **msismod = NULL;       /* Modification variables for MSIS         [Units: {sfu, sfu, unitless} */
  double (*RITJ)[3][3] = NULL;   /* ITRF93 to J2000 3x3 rotation matrix at all times           [Units: ???]  */
  double (*resun)[3] = NULL;     /* Earth to Sun position vector at all times                  [Units: {km, km, km}] */
  double (*remoon)[3] = NULL;    /* Earth to Moon position vector at all times                 [Units: {km, km, km}] */
  double (*CMAT)[MAXGDO] = NULL; /* Constants for calculating spherical harmonic gravity field [Units: unitless] */
  double (*SMAT)[MAXGDO] = NULL; /* Constants for calculating spherical harmonic gravity field [Units: unitless] */
  double state0[6];              /* Initial state vector from Input File     [Units: {km, km, km, km/s, km/s, km/s] */

  struct initCd_struct *initCd = NULL;          /* Structure for drag variables */
  struct RSM_struct *RSMdata = NULL;            /* Structure for Response Surface Model data */
  struct sat_struct *psatout=NULL;              /* Structure for satellite output variables */
  struct constants_struct *psatconstant=NULL;   /* Structure for satellite output constants */
  struct msis_struct *pmsis=NULL;               /* Structure for MSIS variables */
  struct densities_struct *input=NULL;          /* Structure for atmospheric properties */

  /* Initialize all of the propagator variables into memory */
  propagator_init(atm_dir, data_dir, config, get_data, record_data, DioramaDeltaT, natm, &n, &dnp, 
		  &numpoints, &scn, &nobs, kernelpath, &t, &doubletime, &deltat, state0, &psig, &vsig, 
		  &state_cart, &state_elem, &initstate, &msismod, &RITJ, &resun, &remoon, &CMAT, &SMAT, &initCd, 
		  &RSMdata, &psatout, &psatconstant, &pmsis, &input, &rr);

#if PARALLEL
  satperproc = (double)n/(double)nProcs;
  satloopmin = rank*satperproc;
  satloopmax = (rank+1)*satperproc;
#else /* PARALLEL */
  satperproc = 0;
  satloopmin = 0;
  satloopmax = n;
#endif /* PARALLEL */

  /* Loop over satellites */
  for(isat=satloopmin; isat<satloopmax; isat++) {

    /* Perform the propagation integration */
    propagator_core(atm_dir, data_dir, config, get_data, initSV, DioramaDeltaT, finalSV, psig, vsig, 
    		    state_cart, state_elem, state0, initstate, t, deltat, doubletime, RITJ, resun, remoon, 
		    CMAT, SMAT, msismod, numpoints, natm, isat, n, dnp, nobs, input, initCd, psatout, 
		    pmsis, RSMdata, rr);

    /* Write out the orbital trajectory */
    propagator_write(atm_dir, data_dir, record_data, numpoints, psatconstant, t, psatout, scn, isat, state_cart);

  }

  /* Free allocated memory and unload kernels */
  propagator_close(initCd, RSMdata, t, doubletime, CMAT, SMAT, input, RITJ, resun, remoon, 
		   psatconstant, psatout, state_cart, state_elem, initstate, natm, dnp, n,
		   msismod, pmsis, rr, kernelpath);
 
}

/* Initialize the propagator input data */
void propagator_init(char *atm_dir, const char *data_dir, struct config_struct config, data_read get_data, 
		     data_store record_data, double DioramaDeltaT,
		     int *natm, int *n, int *dnp, int *numpoints, int *scn, int *nobs, 
		     char kernelpath[NKER][PATH_MAX], 
		     double **t, double **doubletime, double *deltat, double state0[6], double *psig, double *vsig, 
		     double ***state_cart, double ***state_elem, double ***initstate, double ***msismod, 
		     double (**RITJ)[3][3], double (**resun)[3], double (**remoon)[3], double (**CMAT)[MAXGDO], 
		     double (**SMAT)[MAXGDO], struct initCd_struct **initCd, struct RSM_struct **RSMdata, 
		     struct sat_struct **psatout, struct constants_struct **psatconstant, struct msis_struct **pmsis, 
		     struct densities_struct **input, gsl_rng **rr) {

  Use2Body          = config.Use2Body;            /* Use Two-body Gravity - Should Always Be ON */
  UseSphHarmon      = config.UseSphHarmon;        /* Use Spherical Harmonic Gravity Field */
  GDO               = config.GDO;                 /* Degree and Order of Gravity Field */
  EarthGrav         = config.EarthGrav;           /* Earth Gravity Model - JGM3 or EGM96 */
  UseSRP            = config.UseSRP;              /* Use Solar Radiation Pressure */
  UseDrag           = config.UseDrag;             /* Use Atmospheric Drag */
  DensityModel      = config.DensityModel;        /* Atmospheric Density Model - CIRA72, MSIS, GITM, USSA1976, HASDM */ 
  DynamicMSIS       = config.DynamicMSIS;         /* Compute MSIS properties on the fly */
  ReadSW            = config.ReadSW;              /* Read in tabular space weather indices */
  RSM               = config.RSM;                 /* Use Response Surface Model for Satellite Drag - GRACE, CHAMP */
  UseHFCd           = config.UseHFCd;             /* Use closed-form solution for Cd for simple satellite geometries */
  UseSun            = config.UseSun;              /* Use Third-body Perturbation of the Sun */
  UseMoon           = config.UseMoon;             /* Use Third-body Perturbation of the Moon */
  DensityTimeLookup = config.DensityTimeLookup;   /* Method for Density Interpolation - Nearest Neighbor, Linear */
  RungeKutta        = config.RungeKutta;          /* Runge-Kutta Integration Order - 4th, 8th */
  StateEnsembles    = config.StateEnsembles;      /* Read Initial Satellite State from Text File */
  Keplerian         = config.Keplerian;           /* Compute Keplerian Elements */

  if (get_data == NULL) {  // overrride
    DensityModel = config.DensityModel= 1;  // MSIS
    DynamicMSIS  = config.DynamicMSIS = 1;
  }

  /* Assign number of vertical (radial), latitudinal, and longitudinal cells for each atmospheric model */
  if (DensityModel==4) {        /* HASDM */
    NVERT = 30;
    NLAT = 36;
    NLON = 72;
  } else if (DensityModel==3) { /* US STANDARD ATMOSPHERE, 1976 */
    NVERT = 429;
    NLAT  = 1;
    NLON  = 1;
  } else if (DensityModel==2) { /* GITM */
    NVERT = 50;
    NLAT  = 36;
    NLON  = 72;
  } else if (DensityModel<=1) { /* MSIS & CIRA 72*/
    NVERT = 161;
    NLAT  = 37;
    NLON  = 73; 
  }

  if (RungeKutta==4)
    RKN = 4;
  else if (RungeKutta==8)
    RKN = 13;

#if PARALLEL
  /***************************** MPI Initialization *******************************/
  /*Create communicator for the I/O functions */
  MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
  /* Find my ID number */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  /* Find the total number of procs */
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs); 
#endif /* PARALLEL */

  /* Setup Random Number Generator */
  long int seed;                     /* Seed for GSL random number generator */
  const gsl_rng_type * T;            /* Type for GSL random number generator */
  gsl_rng_env_setup();
  T = gsl_rng_default;
  *rr = gsl_rng_alloc (T);
#if PARALLEL
  seed = time (NULL) * (rank+1);
#else /* PARALLEL */
  seed = time (NULL);
#endif /* DYNAMIC SEED */
  gsl_rng_set(*rr, seed);

  int i;                           /* General index */
  int maxt;                        /* Max elapsed time of simulation in seconds */
  char *utc0 = NULL;               /* ISO formatted UTC string */
  double a0, e0, i0, om0, Om0, f0; /* Keplerian elements */
  double keplerian0[6];            /* Kepelerian state vector */
  double mean, sigma;              /* Mean and standard deviation for satellite ballistic coefficient distribution */

  /* Constant Space Weather inputs for MSIS */
  double cf107;                    /* Constant daily F10.7 solar flux index */
  double cf107A;                   /* Constant 81-day average F10.7 solar flux index */
  double cAp;                      /* Constant daily averaged Ap geomagnetic index */

  /* Allocate memory for drag coefficient structures */
  *initCd =  (struct initCd_struct *) calloc(1, sizeof(struct initCd_struct));
  *RSMdata = (struct RSM_struct *)    calloc(1, sizeof(struct RSM_struct));

  read_input_file(data_dir, &a0, &e0, &i0, &Om0, &om0, &f0, &utc0, &maxt, deltat, scn, n, &sigma, &mean, psig, vsig, &cf107, &cf107A, &cAp, *initCd);

#if DIORAMA
  double et0;
  *n = 1; /* DIORAMA will only run one satellite at a time */
#else
  int j, k; 
  double lt;
  SpiceDouble et0;
#endif

  *numpoints = fabs(maxt/(*deltat)+1); /* Number of time points for full timestep */
  *dnp = 2*(*numpoints)-1;             /* Number of time points for half timestep */

  double moonstate[6];                 /* Earth-Moon position vector    [Units: {km, km, km}] */
  double sunstate[6];                  /* Eart-Sun position vector      [Units: {km, km, km}] */

  for(i=0; i<6; i++) {
    moonstate[i] = 0.0;
    sunstate[i] = 0.0;
  }

  /* Initialize RITJ, resun, and remoon arrays for each half timestep */
  *RITJ = calloc(*dnp , sizeof **RITJ);     /* ITRF93_TO_J2000 3x3 rotation matrix at all times */
  *resun = calloc(*dnp, sizeof **resun);    /* Earth to Sun position vector at all times */
  *remoon = calloc(*dnp, sizeof **remoon);  /* Earth to Moon position vector at all times */

  *CMAT = calloc( (MAXGDO), sizeof **CMAT); /* Gravity constants for Spherical Harmonic Gravity Field */
  *SMAT = calloc( (MAXGDO), sizeof **SMAT); /* Gravity constants for Spherical Harmonic Gravity Field */
  
#if !DIORAMA
  SpiceDouble mat[3][3]; /* ITRF93_TO_J2000 3x3 rotation matrix */
#endif

  *t = (double *) calloc(*numpoints, sizeof(double));
  *doubletime = (double *) calloc(*dnp+2, sizeof(double));

  /* Initialize state_cart and state_elem vectors for each full timestep */
  /* state_cart is cartesian position and velocity vectors (units of km and km/s) */
  /* state_elem is classical orbit elements (a,e,i,Om,om,f) (units of km and rad) */
  *state_cart = (double **) calloc(*dnp, sizeof(double));
  *state_elem = (double **) calloc(*dnp, sizeof(double));
  for(i=0; i<*dnp; i++) {
    (*state_cart)[i] = (double *) calloc(6, sizeof(double));
    (*state_elem)[i] = (double *) calloc(6, sizeof(double));
  }

  /* Initialize array for atmospheric ensemble number */
  natm = (int *) calloc(*n, sizeof(int));

  *msismod = (double **) calloc(*n, sizeof(double));
  for(i=0; i<*n; i++) {
    (*msismod)[i] = (double *) calloc(3, sizeof(double));
  }
  
  /* Initialize array for initial state */
  *initstate = (double **) calloc(*n, sizeof(double));
  for(i=0; i<*n; i++) {
    (*initstate)[i] = (double *) calloc(6, sizeof(double));
  }
  
  *nobs = 0;
   
  /* Initialize satellite structure to save all output */
  *psatout = (struct sat_struct *) calloc(*numpoints, sizeof(struct sat_struct));
  *psatconstant = (struct constants_struct *) calloc(1, sizeof(struct constants_struct));

  if(ReadSW) {
    *nobs = read_sw_lines(data_dir);
  } else {
    *nobs = 1;
  }
  *pmsis = (struct msis_struct *) calloc(*nobs, sizeof(struct msis_struct));

#if !DIORAMA

  /* Define kernels to load */
  char kernels[NKER][PATH_MAX] = { "spice_kernels/naif0010.tls",
				   "spice_kernels/earth_000101_180101_171010.bpc",
				   "spice_kernels/earth_720101_070426.bpc",
				   "spice_kernels/de421.bsp" };

  if( UseSphHarmon || UseSun || UseMoon || UseSRP || UseDrag) {

 
    for(i=0; i<NKER; i++) {   
      strcpy(kernelpath[i], data_dir);
      strcat(kernelpath[i], "/");
      strcat(kernelpath[i], kernels[i]);
    }
    
    /* Load kernels */
    for(i=0; i<NKER; i++) {
      furnsh_c(kernelpath[i]);
    }

  }
#endif

  /* Define constants structure for future output */
  create_constants_structure(*psatconstant, *initCd);

  if((DensityModel==1 && DynamicMSIS) || (DensityModel==4 && ReadSW)) {

    read_space_weather(data_dir, *pmsis, *nobs, cf107, cf107A, cAp);

  } /* DensityModel && DynamicMSIS */

  if(StateEnsembles) {

  /* Read in initial state for satellite */
    read_state_ensemble(data_dir, *n, *scn, *initstate, natm, *msismod);

  } else {

    for(i=0; i<*n; i++) {
      natm[i] = 1;
    }

  }

  if (RSM) {

    /* Read in RSM data */
    read_RSM(data_dir, *RSMdata);

  }

  if (UseSphHarmon) {

    /* Read in Earth gravity field coefficients */
    read_earth_grav(data_dir, *CMAT, *SMAT);

  } /* UseSphHarmon */

  /* Convert time to ET */
  /* ET = "ephemeris time", aka solar-system barycentric dynamical time, which */
  /* is used by JPL functions, which is essentially seconds since 01/01/2000 */
  
  if( UseSphHarmon || UseSun || UseMoon || UseSRP || UseDrag ) {

#if !DIORAMA
    str2et_c(utc0, &et0);
#endif

  } else { 
    et0 = 0.0;
  }

  /* Define the propagation time */
  for(i=0; i<*numpoints; i++) {
    (*t)[i] = et0 + i*(*deltat);
  }

  /* Allocate data for atmospheric data structure */
  if (!DynamicMSIS) {
    *input = calloc(1, sizeof(struct densities_struct));
  }

  
  if(!StateEnsembles) {

    /* Define Keplerian Element Vector*/
    keplerian0[0] = a0;
    keplerian0[1] = e0;
    keplerian0[2] = i0 * DEG2RAD;
    keplerian0[3] = Om0* DEG2RAD;
    keplerian0[4] = om0* DEG2RAD;
    keplerian0[5] = f0 * DEG2RAD;
    
    /* Convert to cartesian state */
    kep2cart(keplerian0, state0);

  } else {

    for(i=0; i<6; i++) {
      keplerian0[i] = 0.0;
      state0[i] = 0.0;
    }

  }

  for(i=0; i<*dnp; i++) {

    (*doubletime)[i] = et0 + i*(*deltat)/2.0;   /* First make time vector again but with twice as many points */

#if !DIORAMA
    if( UseSphHarmon || UseSun || UseMoon || UseSRP || UseDrag ) {

      pxform_c("J2000", "ITRF93", (*doubletime)[i], mat);   /* Get the earth-orientation rotation matrix */
      
      spkezr_c("SUN", (*doubletime)[i], "J2000", "NONE", "EARTH", sunstate, &lt);  /* Get the Earth-Sun position vector */
      spkezr_c("MOON", (*doubletime)[i], "J2000", "NONE", "EARTH", moonstate, &lt);/* Get the Earth-Moon position vector*/
      
      /* Assign position vectors and rotation matrices to time-series arrays */
      for(j=0; j<3; j++) {
	(*resun)[i][j] = sunstate[j];
	(*remoon)[i][j] = moonstate[j];
	for(k=0; k<3; k++) {
	  (*RITJ)[i][j][k] = mat[j][k];
	}
      }

    }
#endif
    
  }

  /* Define special value so we can easily find end of the array */
  (*doubletime)[*dnp] = -1;
  (*doubletime)[*dnp+1] = -1;
  
}



void propagator_core(char *atm_dir, const char *data_dir, struct config_struct config, data_read get_data, 
		     double initSV[6], double deltaT, double finalSV[6], double psig, double vsig, 
		     double **state_cart, double **state_elem, double state0[6],
		     double **initstate, double *t, double deltat, double *doubletime,  
		     double (*RITJ)[3][3], double (*resun)[3], double (*remoon)[3], double (*CMAT)[MAXGDO], 
		     double (*SMAT)[MAXGDO], double **msismod, int numpoints, int *natm, int isat, int n, 
		     int dnp, int nobs, struct densities_struct *input, struct initCd_struct *initCd, 
		     struct sat_struct *psatout, struct msis_struct *pmsis, struct RSM_struct *RSMdata, 
		     gsl_rng *rr) {

  int i;                /* General index       */
  int atmcounter;       /* Atmospheric counter */
  int it;               /* Timestep counter    */

  double var[6];        /* Variances for orbital state vector */
  double pvsig[6];      /* Standard deviations for orbital state vector */
  double OrbitT[6];     /* Perturbation from mean orbital state vector   [Units: {km, km, km, km/s, km/s, km/s}] */
  double kepstate[6];   /* Keplerian state vector                        [Units: {km, unitlesss, rad, rad, rad, rad}]*/
  double nstate[6];     /* Updated state vector through RK integration   [Units: {km, km, km, km/s, km/s, km/s}] */
  double rvec[3];       /* Position vector                               [Units: {km, km, km}] */
  double vvec[3];       /* Velocity vector                               [Units: {km/s, km/s, km/s}] */
  double et;            /* Ephemeris Time (seconds since J2000)          [Units: seconds] */

  struct sat_struct *satout; /* Structure for satellite propagation output variables */

  /* Initialize state vectors */
  for(i=0; i<6; i++) {
    kepstate[i] = 0.0;
    nstate[i] = 0.0;
  }

  /* Define position and velocity standard deviations */
  pvsig[0] = pvsig[1] = pvsig[2] = psig;
  pvsig[3] = pvsig[4] = pvsig[5] = vsig;

  /* Compute perturbation on mean state */
  for(i=0; i<6; i++) {
    var[i] = ranfgs(rr, 1.0);
    OrbitT[i] = var[i]*pvsig[i];
  }
  
  /* Apply uncertainty in state */
  for(i=0; i<6; i++) {
#if DIORAMA
    state_cart[0][i] = initSV[i];
#else /* DIORAMA */
    if(StateEnsembles) {
      state_cart[0][i] = initstate[isat][i];
    } else {
      state_cart[0][i] = state0[i] + OrbitT[i];
    }
#endif /* DIORAMA */
  }
  
  /* Define position vector from state vector */
  for(i=0; i<3; i++) {
    rvec[i] = state_cart[0][i];
  }
  
  /* Define velocity vector from state vector */
  for(i=3; i<6; i++) {
    vvec[i-3] = state_cart[0][i];
  }

  /* Convert cartesian state to Keplerian state */
  if(Keplerian) {    
    cart2kep(rvec, vvec, kepstate);
    for(i=0; i<6; i++) {
      state_elem[0][i] = kepstate[i];
    }
  }
  
  atmcounter = 0;
  
  /********************* Propagate the cartesian state ***********************/
  for(it=0; it<numpoints; it++) {
    
    et = t[it];   /* Current Ephemeris Time */

    /* Import atmospheric properties from data files */
    if (UseDrag && !DynamicMSIS) {

      if (get_data != NULL) {
	if(it==0) {
	  get_data(atm_dir, data_dir, config, numpoints, input, t, natm, isat, et);
	  atmcounter++;
	} else if (et >= input->etlist[atmcounter]) {
	  get_data(atm_dir, data_dir, config, numpoints, input, t, natm, isat, et);
	  atmcounter++;
	}
      }

    }
    
    /* Define satellite variables in satellite output structure */
    satout = psatout + it;   
    for(i=0; i<3; i++) {
      satout->pos[i] = state_cart[it][i];
      satout->vel[i] = state_cart[it][i+3];
      satout->m = initCd->sat_mass;
      satout->Cr = initCd->SRP_Cr;
    }
    
    if (RungeKutta == 4) {
      rk4(rhseom, deltat, state_cart[it], t[it], n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT,
	  input, dnp, psatout, it, pmsis, initCd, msismod, isat, atmcounter, RSMdata, nobs);
    } 
    
    if (RungeKutta == 8) {
      rk8(rhseom, deltat, state_cart[it], t[it], n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT,
	  input, dnp, psatout, it, pmsis, initCd, msismod, isat, atmcounter, RSMdata, nobs);
    } 

    /* Update state vector */
    for(i=0; i<6; i++) {
      state_cart[it+1][i] = nstate[i];
    }
    
    /* Define position vector from state vector */
    for(i=0; i<3; i++) {
      rvec[i] = state_cart[it+1][i];
    }
    
    /* Define velocity vector from state vector */
    for(i=3; i<6; i++) {
      vvec[i-3] = state_cart[it+1][i];
    }
    
    /* Compute Keplerian elements from Cartesian State Vector */
    if(Keplerian) {    
      cart2kep(rvec, vvec, kepstate);  
      for(i=0; i<6; i++) {
	state_elem[it+1][i] = kepstate[i];
	satout->kep[i] = kepstate[i];
      }    
    } else {  
      for(i=0; i<6; i++) {
	satout->kep[i] = 0.0;
      }   
    }
    
  } /******************************* END TIME LOOP ******************************/

#if DIORAMA
  for(i=0; i<6; i++) {
    finalSV[i] = state_cart[numpoints-1][i];
  }
#endif /* DIORAMA */
  
}


void propagator_write(char *atm_dir, const char *data_dir, data_store record_data, int numpoints, 
			struct constants_struct *psatconstant,double *t, struct sat_struct *psatout,
			int scn, int isat, double **state_cart) {

  int it;
  char outfilename[100];

  if (record_data != NULL) {
    record_data(atm_dir, psatconstant, numpoints, t, psatout);
  }
  else {
    /* Open file for data output - Temporary until HDF5 is incorporated */
    sprintf(outfilename, "%s/propdata_%d_%d.dat", data_dir, scn, isat);
    FILE *fout = fopen(outfilename, "w");
    if (fout == NULL) {
      printf("%s: %s\n", strerror(errno), outfilename);
      exit(1);
    }
    
    /* Write X, Y, Z, U, V, and W to a file as a function of time */
    printf("Writing Cartesian Position and Velocity for Satellite #%d and Ensemble #%d to File %s\n", scn, isat, outfilename);
    fprintf(fout, "Time [s]     X [km]     Y [km]     Z [km]     U [km/s]     V [km/s]     W [km/s]\n");
    for(it=0; it<numpoints; it++) {
      fprintf(fout, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e \n", t[it]-t[0],
	      state_cart[it][0], state_cart[it][1], state_cart[it][2],
	      state_cart[it][3], state_cart[it][4], state_cart[it][5]);
    }
    
    fclose(fout);
  }
  
}


/* Free all of the allocated memory and unload cspice kernels*/
void propagator_close(struct initCd_struct *initCd, struct RSM_struct *RSMdata, double *t, double *doubletime,
		      double (*CMAT)[MAXGDO], double (*SMAT)[MAXGDO], struct densities_struct *input, 
		      double (*RITJ)[3][3], double (*resun)[3], double (*remoon)[3], 
		      struct constants_struct *psatconstant, struct sat_struct *psatout,
		      double **state_cart, double **state_elem, double **initstate, int *natm, int dnp, int n, 
		      double **msismod, struct msis_struct *pmsis, gsl_rng *rr, char kernelpath[NKER][PATH_MAX]) {

  int i;

  free(initCd);
  free(RSMdata);

  free(t);
  free(doubletime);

  free(CMAT);
  free(SMAT);

  if (!DynamicMSIS) {
    free(input);
  }

  free(RITJ);

  free(resun);
  free(remoon);

  free(psatconstant);
  free(psatout);

  for(i=0; i<dnp; i++) {
    free(state_cart[i]);
    free(state_elem[i]);
  }

  free(state_cart);
  free(state_elem);
  
  if(StateEnsembles) {
    for(i=0; i<n; i++) {
      free(initstate[i]);
    }
    
    free(initstate);

    free(natm);
  }

  for(i=0; i<n; i++) {
    free(msismod[i]);
  }

  free(msismod);

  free(pmsis);

  gsl_rng_free(rr);

#if !DIORAMA
  if( UseSphHarmon || UseSun || UseMoon || UseDrag) {

    /* Unload kernels */
    for(i=0; i<NKER; i++) {
      unload_c(kernelpath[i]);
    }

  }
#endif
  
}




