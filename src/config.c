#define HAVE_CONFIG 1
#include "defs.h"

/* Configuration values */
int Use2Body;     
int UseSphHarmon; 
int GDO;          
int EarthGrav;    
int UseSRP;       
int UseDrag;      
int DensityModel; 
int DynamicMSIS;  
int ReadSW;       
int UseHFCd;      
int RSM;        
int UseSun;       
int UseMoon;      
int DensityTimeLookup; 
int RungeKutta;   
int StateEnsembles; 
int Keplerian;

int NVERT;          /* Number of altitudinal cells */
int NLAT;           /* Number of latitudinal cells */
int NLON;           /* Number of longitudinal cells */
int RKN;            /* Number of Runge-Kutta iterations */


int rkcounter    = 0;
