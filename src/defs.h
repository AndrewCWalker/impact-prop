#ifndef DEFS_H
#define DEFS_H
/******************************************************************************/
/*************************PROPAGATOR CONTROLS *********************************/
/******************************************************************************/
#ifndef HAVE_CONFIG
extern int Use2Body;
extern int UseSphHarmon;
extern int GDO;
extern int EarthGrav;
extern int UseSRP;
extern int UseDrag;
extern int DensityModel;
extern int DynamicMSIS;
extern int ReadSW;
extern int UseBCscale;
extern int UseHFCd;
extern int RSM;
extern int UseSTimeAvg;
extern int UseSun;
extern int UseMoon;
extern int DensityTimeLookup;
extern int RungeKutta;
extern int StateEnsembles;
extern int Keplerian;

extern int NVERT;
extern int NLAT;
extern int NLON;
extern int RKN;
#endif


#ifndef PARALLEL
#define PARALLEL       0                /* Turns on MPI parallelization for multiple satellites */
#endif

#ifndef DIORAMA
#define DIORAMA        0                /* Special compiling for DIORAMA */
#endif

#define DIORAMA_NTS     3000             /* Number of DIORAMA timesteps */

/************************************************************************/
/****************************CONVERSION FACTORS**************************/
/************************************************************************/
#define DEG2ARCSEC     3600.0
#define ARCSEC2DEG     (1.0/DEG2ARCSEC)
#define RAD2DEG        (180.0/M_PI)
#define DEG2RAD        (1.0/RAD2DEG)
#define ARCSEC2RAD     (ARCSEC2DEG*DEG2RAD)
#define RAD2ARCSEC     (1.0/ARCSEC2RAD)
#define CSOLV          299792.458      /* Speed of light in a vacuum [km/s] */


/********************************************************************************/
/***********************EARTH & GRAVITY PARAMETERS*******************************/
/********************************************************************************/
#define JGM3Re         6378.137         /* Equatorial Radius of the Earth [km], JGM-3 */
#define JGM3Rp         6356.7516        /* Polar Radius of the Earth [km], Vallado */
#define JGM3GMe        3.986004418e5    /* GM of Earth [km^3/s^2], JGM-3 */
#define JGM3ome        7.292115e-5      /* Omega (angular velocity) of the Earth [rad/s] - (Is this really JGM-3?) */
#define EGM96GMe       3.986004418e5    /* GM of Earth [km^3,s^2], EGM96 */ 
#define EGM96Re        6378.137         /* Equatorial Radius of the Earth [km], EGM96 */
#define MoonGM         4.902801e3       /* GM of Moon [km^3/s^2], Montenbruck */
#define SunGM          1.32712440018e11 /* GM of Sun [km^3/s^2], Montenbruck */
#define EECC           0.081819221456   /* Earth Eccentricity, Vallado*/


/********************************************************************************/
/**************************ATMOSPHERIC PARAMETERS**********************************/
/********************************************************************************/

#define MAXVERT        429
#define MAXLAT         37
#define MAXLON         73
#define MAXSPECIES     6
#define MAXGDO         27


/********************************************************************************/
/*****************************OTHER CONSTANTS************************************/
/********************************************************************************/
#define NSPECIES       6                /* Number of species */
#define NKER           4                /* Number of CSPICE kernels */
#define MAXFILES       500               /* Maximum number of files that can be read in for propagation */
#define length(arr)   (sizeof(arr)/sizeof(arr[0]))


#endif
