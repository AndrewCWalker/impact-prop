#include <stdio.h>
#include <limits.h>
#include <string.h>

#include "defs.h"
#include "structs.h"
#include "propagator.h"
#include "io.h"

#if PARALLEL
#include "mpi.h"
#endif /* PARALLEL */

#ifndef DATA_DIR
#define DATA_DIR "../data/"
#endif


void read_density(char *atm_dir, const char *data_dir, struct config_struct config,
		  int numpoints, struct densities_struct* densities, double *t, int *natm, 
		  int isat, double et) {

  if (config.DensityModel==1 || config.DensityModel==2 || config.DensityModel==4) { /* GITM or MSIS */
    
    double tsearch[2];
    
      /* Define start and end points of simulation */
    tsearch[0] = t[0];
    tsearch[1] = t[numpoints-1] + 3600.0*3.0; /* Add an extra 3 hours buffer */
    
    read_files(tsearch, atm_dir, densities->etlist, densities->files, &densities->Nfiles);
    
#if 0 //FIXME
    if (config.DensityModel==1) /* MSIS */
    	importmsis(atm_dir, densities->etlist, densities->files, densities->AltKm, densities->LatDeg, densities->LonDeg, densities->Rho, densities->ndenvec, densities->Temperature, densities->Nfiles);
#endif    

    if (config.DensityModel==2) /* GITM */
      importgitm(atm_dir, densities->etlist, et, densities->files, densities->AltKm, densities->LatDeg, densities->LonDeg, densities->Rho, densities->ndenvec, densities->Temperature, densities->Nfiles, natm, isat);
    
    if (config.DensityModel==4) /* HASDM */
      importhasdm(atm_dir, densities->etlist, et, densities->files, densities->AltKm, densities->LatDeg, densities->LonDeg, densities->Rho, densities->ndenvec, densities->Temperature, densities->Nfiles);
    
  }
  
  if (config.DensityModel==3) /* US STANDARD ATMOSPHERE */
    importussa76(data_dir, densities->AltKm, densities->Rho[0][0][0]);
  
}


int main(int argc, char *argv[]) {

#if DIORAMA
  int i, j;
#endif /* DIORAMA */

  struct config_struct config;

  /* Define dummy variables only used if DIORAMA flag is on */
  double initSV[6];
  double finalSV[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double DioramaDeltaT = 180.0;

  /* Initialization for DIORAMA test */
  initSV[0] = -2.661791928199999893e+04;
  initSV[1] =  1.798946516100000053e+03;
  initSV[2] = -1.994912198999999831e+01;
  initSV[3] = -1.631212321299999879e-01;
  initSV[4] = -2.212360033700000006e+00;
  initSV[5] =  3.154939023599999892e+00;
    
  read_config(DATA_DIR, &config);

#if PARALLEL
  MPI_Init(&argc, &argv);  /* Initialize MPI */
#endif /* PARALLEL */

  /* Set data paths for atmospheric data files */
  char atm_dir[PATH_MAX] = {0};
  if (config.DensityModel==4) /* HASDM */
    strcpy(atm_dir, "/data0/IMPACT/data/HASDM/");
  if (config.DensityModel==2) /* GITM */
    strcpy(atm_dir, "/data0/IMPACT/data/analysis/");
  if (config.DensityModel==1) /* MSIS (non-dynamic) */
    strcpy(atm_dir, "/data0/tmp/NEWdata/SimpleMSIS/");

  propagate(atm_dir, DATA_DIR, config, &read_density, NULL, initSV, DioramaDeltaT, finalSV);

#if PARALLEL
  MPI_Finalize();
#endif /* PARALLEL */
 
  return 0;
}
