#ifndef MISC_H
#define MISC_H


double hellipsoid(double rvec[]);

int findnearest(double et, double dt, double et_doubletime[], int dnp);

double dot(double V[], double W[]);

void cross(double V[], double W[], double VEC[]);

double norm(double V[]);

void matrix_multiply(double A[3][3], double V[3], double VEC[3]);

double ranf0(gsl_rng *rr);

double ranfgs(gsl_rng *rr, double sigma);


#endif
