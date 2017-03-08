#ifndef TRANSFORMSTION_H
#define TRANSFORMSTION_H

void cart2sph(double x, double y, double z, double *R, double *theta, double *phi);

void sph2cart(double theta, double phi, double R, double *x, double *y, double *z);

void cart2kep(double rvec[3], double vvec[3], double kep[6]);

void kep2cart(double kep[6], double state[6]);


#endif
