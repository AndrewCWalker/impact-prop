#ifndef IO_H
#define IO_H

#include "structs.h"
#include "defs.h"

void read_files(double tsearch[2], char *dirname, double etlist[MAXFILES], char files[MAXFILES][100], int *Nfiles);

void importgitm(char *dirname, double etlist[MAXFILES], double et, char files[MAXFILES][100], double AltKm[MAXVERT], double LatDeg[MAXLAT], double LonDeg[MAXLON], double Rho[2][MAXLON][MAXLAT][MAXVERT], double ndenvec[MAXSPECIES][2][MAXLON][MAXLAT][MAXVERT], double Temperature[2][MAXLON][MAXLAT][MAXVERT], int Nfiles, int *natm, int isat);

void importmsis(char *dirname, double etlist[MAXFILES], char files[MAXFILES][100], double AltKm[MAXVERT], double LatDeg[MAXLAT], double LonDeg[MAXLON], double Rho[2][MAXLON][MAXLAT][MAXVERT], double ndenvec[MAXSPECIES][2][MAXLON][MAXLAT][MAXVERT], double Temperature[2][MAXLON][MAXLAT][MAXVERT], int Nfiles);

void importhasdm(char *dirname, double etlist[MAXFILES], double et, char files[MAXFILES][100], double AltKm[MAXVERT], double LatDeg[MAXLAT], double LonDeg[MAXLON], double Rho[2][MAXLON][MAXLAT][MAXVERT], double ndenvec[MAXSPECIES][2][MAXLON][MAXLAT][MAXVERT], double Temperature[2][MAXLON][MAXLAT][MAXVERT], int Nfiles);

void read_earth_grav(const char *data_dir, double (*CMAT)[MAXGDO], double (*SMAT)[MAXGDO]);

int read_num_lines(const char *data_dir, char *filename);

void read_input_file(const char *data_dir, double *a0, double *e0, double *i0, double *Om0, double *om0, double *f0, char **utc0, int *maxt, double *deltat, int *scn, int *n, double *sigma, double *mean, double *psig, double *vsig, double *cf107, double *cf107A, double *cAp, struct initCd_struct *initCd);

void read_config(const char *data_dir, struct config_struct *config);

void importussa76(const char *data_dir, double *altitude, double *density);

void create_constants_structure(struct constants_struct *psatconstant, struct initCd_struct *initCd);

void read_state_ensemble(const char *data_dir, int nsat, int scn, double **initstate, int *natm, double **msismod);

void read_RSM(const char *data_dir, struct RSM_struct *RSMdata);

#endif
