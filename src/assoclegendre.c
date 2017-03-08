#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_gamma.h>
#include <time.h>
#include "defs.h"
#include "assoclegendre.h"


void assoclegendre(double phigc, double P[][MAXGDO]) {

/*     Calculate the associated legendre functions, where phigc is the          */
/*     satellite's geocentric latitude in radians, and N is the degree and      */
/*     order. The return matrix contains the resulting function values, where   */
/*     the python indices correspond with the degree and order numbers          */
/*     (e.g. P(n=4,m=0) would be accessed in python as P[4,0]). Note also that  */
/*     this matrix indexing is simple, but the return matrix is somewhat sparse */
/*     (has zeros in the upper triangular part).                                */
    
/*     n = degree (letter "l" in Vallado, letter "n" in Montebruck)             */
/*     m = order                                                                */
    
/*     note: important to initialize with 0, because the P_n_minus_2_m term is  */
/*     required below, and Vallado (pp. 557) says this term should be 0. (so    */
/*     don't initiailze as, e.g. NaN).                                          */
  
  int n, m;
  int N = GDO;

  for(n=0; n<MAXGDO; n++) {
    for(m=0; m<MAXGDO; m++) {
      P[n][m] = 0.0;
    }
  }

  /* DEFINE LOWEST DEGREE AND ORDER POLYNOMIALS */
  P[0][0] = 1.0;
  P[1][0] = sin(phigc);
  P[1][1] = cos(phigc);

  /* USE RECURSION RELATIONSHIPS TO CALCULATE HIGHER POLYNOMIALS */
  for(n=2; n<N+1; n++) { 
    P[n][0] = ((2.0*n-1.0)*P[1][0]*P[n-1][0] - (n-1.0)*P[n-2][0])/n;    
    P[n][n] = (2.0*n - 1.0)*P[1][1]*P[n-1][n-1];
  }

  for(m=1; m<N+1; m++) {
    for(n=m+1; n<N+1; n++) {
      P[n][m] = P[n-2][m] + (2.0*n-1.0)*P[1][1]*P[n-1][m-1];
    }
  }

}

  

double normcoef(int L, int M) {

/*    Calculate the normalization parameter (Pi) for deg (L) and order (M) for  */
/*    a grav coeficient or legendre polynomial (see Vallado pp. 519).           */

  double k, Pi;

  if(M==0) {
    k=1.0;
  } else {
    k=2.0;
  }

  Pi = sqrt(gsl_sf_fact(L+M)/(gsl_sf_fact(L-M)*k*(2.0*L+1.0)));

  return Pi;

}
 

  
    

