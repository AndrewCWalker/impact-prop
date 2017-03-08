#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include "defs.h"
#include "structs.h"
#include "rhseom.h"
#include "rk.h"


void rk4(rhseom_t rhseom, 
	 double dt, double x[6], double t, int n, 
	 double nstate[6], double doubletime[], 
	 double (*resun)[3], double (*remoon)[3], double (*RITJ)[3][3], 
	 double (*CMAT)[MAXGDO], double (*SMAT)[MAXGDO], struct densities_struct *input, 
	 int dnp, struct sat_struct *psatout, int it, struct msis_struct *pmsis,
	 struct initCd_struct *initCd, double **msismod, int isat, int atmcounter,
	 struct RSM_struct *RSMdata, int nobs) {

  /* A homemade 4th order Runge Kutta routine */

  int i;
  double k1[6];
  double k2[6];
  double k3[6];
  double k4[6];
  double ns1[6];
  double ns2[6];
  double ns3[6];

  rhseom(x, t, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    ns1[i] = x[i] + nstate[i]*dt/2.0;
    k1[i] = nstate[i];
  }

  rhseom(ns1, t+dt/2.0, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    ns2[i] = x[i] + nstate[i]*dt/2.0;
    k2[i] = nstate[i];
  }

  rhseom(ns2, t+dt/2.0, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    ns3[i] = x[i] + nstate[i]*dt;
    k3[i] = nstate[i];
  }

  rhseom(ns3, t+dt, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    k4[i] = nstate[i];
  }

  for(i=0; i<6; i++) {
    nstate[i] = x[i] + dt*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
  }

}



void rk8(rhseom_t rhseom, 
	 double dt, double x[6], double t, int n, 
	 double nstate[6], double doubletime[], 
	 double (*resun)[3], double (*remoon)[3], double (*RITJ)[3][3], 
	 double (*CMAT)[MAXGDO], double (*SMAT)[MAXGDO], struct densities_struct *input,
	 int dnp, struct sat_struct *psatout, int it, struct msis_struct *pmsis,
	 struct initCd_struct *initCd, double **msismod, int isat, int atmcounter,
	 struct RSM_struct *RSMdata, int nobs) {

  int i;
  double k1[6], k2[6], k3[6], k4[6], k5[6], k6[6], k7[6], k8[6], k9[6], k10[6], k11[6], k12[6], k13[6];
  double ns1[6], ns2[6], ns3[6], ns4[6], ns5[6], ns6[6], ns7[6], ns8[6], ns9[6], ns10[6], ns11[6], ns12[6], ns13[6];

  double h = dt;
   
  double c1 = 13.0 / 288.0;
  double c6 = 32.0 / 125.0;
  double c7 = 31213.0 / 144000.0;
  double c8 = 2401.0 / 12375.0;
  double c9 = 1701.0 / 14080.0;
  double c10 = 2401.0 / 19200.0;
  double c11 = 19.0 / 450.0;

  double a2 = 1.0 / 4.0;
  double a3 = 1.0 / 12.0;
  double a4 = 1.0 / 8.0;
  double a5 = 2.0 / 5.0;
  double a6 = 1.0 / 2.0;
  double a7 = 6.0 / 7.0;
  double a8 = 1.0 / 7.0;
  double a9 = 2.0 / 3.0;
  double a10 = 2.0 / 7.0;
  double a12 = 1.0 / 3.0;

  double b31 = 5.0 / 72.0;
  double b32 = 1.0 / 72.0;
  double b41 = 1.0 / 32.0;
  double b43 = 3.0 / 32.0;
  double b51 = 106.0 / 125.0;
  double b53 = -408.0 / 125.0;
  double b54 = 352.0 / 125.0;
  double b61 = 1.0 / 48.0;
  double b64 = 8.0 / 33.0;
  double b65 = 125.0 / 528.0;
  double b71 = -13893.0 / 26411.0;
  double b74 =  39936.0 / 26411.0;
  double b75 = -64125.0 / 26411.0;
  double b76 =  60720.0 / 26411.0;
  double b81 = 37.0/392.0;
  double b85 = 1625.0/9408.0;
  double b86 = -2.0/15.0;
  double b87 = 61.0/6720.0;
  double b91 = 17176.0 / 25515.0;
  double b94 = -47104.0 / 25515.0;
  double b95 = 1325.0 / 504.0;
  double b96 = -41792.0 / 25515.0;
  double b97 = 20237.0 / 145800.0;
  double b98 = 4312.0 / 6075.0;
  double b10_1 = -23834.0 / 180075.0;
  double b10_4 = -77824.0 / 1980825.0;
  double b10_5 = -636635.0 / 633864.0;
  double b10_6 = 254048.0 / 300125.0;
  double b10_7 = -183.0 / 7000.0;
  double b10_8 = 8.0 / 11.0;
  double b10_9 = -324.0 / 3773.0;
  double b11_1 = 12733.0 / 7600.0;
  double b11_4 = -20032.0 / 5225.0;
  double b11_5 = 456485.0 / 80256.0;
  double b11_6 = -42599.0 / 7125.0;
  double b11_7 = 339227.0 / 912000.0;
  double b11_8 = -1029.0 / 4180.0;
  double b11_9 = 1701.0 / 1408.0;
  double b11_10 = 5145.0 / 2432.0;
  double b12_1 = -27061.0 / 204120.0;
  double b12_4 = 40448.0 / 280665.0;
  double b12_5 = -1353775.0 / 1197504.0;
  double b12_6 = 17662.0 / 25515.0;
  double b12_7 = -71687.0 / 1166400.0;
  double b12_8 = 98.0 / 225.0;
  double b12_9 = 1.0 / 16.0;
  double b12_10 = 3773.0 / 11664.0;
  double b13_1 = 11203.0 / 8680.0;
  double b13_4 = -38144.0 / 11935.0;
  double b13_5 = 2354425.0 / 458304.0;
  double b13_6 = -84046.0 / 16275.0;
  double b13_7 = 673309.0 / 1636800.0;
  double b13_8 = 4704.0 / 8525.0;
  double b13_9 = 9477.0 / 10912.0;
  double b13_10 = -1029.0 / 992.0;
  double b13_12 = 729.0 / 341.0;
   
  //double e1 = -6600.0 / 3168000.0;
  //double e6 = -135168.0 / 3168000.0;
  //double e7 = -14406.0 / 3168000.0;
  //double e8 = 57624.0 / 3168000.0;
  //double e9 = 54675.0 / 3168000.0;
  //double e10 = -396165.0 / 3168000.0;
  //double e11 = -133760.0 / 3168000.0;
  //double e12 = 437400.0 / 3168000.0;
  //double e13 = 136400.0 / 3168000.0;

  double h4 = a2 * h;
  //double h6 = a3 * h;

  /* Compute k1 */
  rhseom(x, t, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    ns1[i] = x[i] + nstate[i]*h4;
    k1[i] = nstate[i];
  }

  /* Compute k2 */
  rhseom(ns1, t+h4, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    k2[i] = nstate[i];
    ns2[i] = x[i] + h * ( b31*k1[i] + b32*k2[i]);
  }

  /* Compute k3 */
  rhseom(ns2, t+a3*h, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    k3[i] = nstate[i];
    ns3[i] = x[i] + h * ( b41*k1[i] + b43*k3[i]);
  }

  /* Compute k4 */
  rhseom(ns3, t+a4*h, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    k4[i] = nstate[i];
    ns4[i] = x[i] +  h * ( b51*k1[i] + b53*k3[i] + b54*k4[i]);
  }

  /* Compute k5 */
  rhseom(ns4, t+a5*h, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    k5[i] = nstate[i];
    ns5[i] = x[i] +  h * ( b61*k1[i] + b64*k4[i] + b65*k5[i]);
  }

  /* Compute k6 */
  rhseom(ns5, t+a6*h, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    k6[i] = nstate[i];
    ns6[i] = x[i] +  h * ( b71*k1[i] + b74*k4[i] + b75*k5[i] + b76*k6[i]);
  }

  /* Compute k7 */
  rhseom(ns6, t+a7*h, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    k7[i] = nstate[i];
    ns7[i] = x[i] +  h * ( b81*k1[i] + b85*k5[i] + b86*k6[i] + b87*k7[i]);
  }

  /* Compute k8 */
  rhseom(ns7, t+a8*h, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    k8[i] = nstate[i];
    ns8[i] = x[i] +  h * ( b91*k1[i] + b94*k4[i] + b95*k5[i] + b96*k6[i] + b97*k7[i] + b98*k8[i]);
  }

  /* Compute k9 */
  rhseom(ns8, t+a9*h, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    k9[i] = nstate[i];
    ns9[i] = x[i] +  h * ( b10_1*k1[i] + b10_4*k4[i] + b10_5*k5[i] + b10_6*k6[i] + b10_7*k7[i] +
			   b10_8*k8[i] + b10_9*k9[i]);
  }

  /* Compute k10 */
  rhseom(ns9, t+a10*h, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    k10[i] = nstate[i];
    ns10[i] = x[i] +  h * ( b11_1*k1[i] + b11_4*k4[i] +b11_5*k5[i] + b11_6*k6[i] + b11_7*k7[i] +
			    b11_8*k8[i] + b11_9*k9[i] + b11_10*k10[i]);
  }

  /* Compute k11 */
  rhseom(ns10, t+h, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    k11[i] = nstate[i];
    ns11[i] = x[i] +  h * (b12_1*k1[i] + b12_4*k4[i] + b12_5*k5[i] + b12_6*k6[i] + b12_7*k7[i] +
			   b12_8*k8[i] + b12_9*k9[i] + b12_10*k10[i]);
  }

  /* Compute k12 */
  rhseom(ns11, t+a12*h, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    k12[i] = nstate[i];
    ns12[i] = x[i] + h * (b12_1*k1[i] + b12_4*k4[i] + b12_5*k5[i] + b12_6*k6[i] + b12_7*k7[i] +
			  b12_8*k8[i] + b12_9*k9[i] + b12_10*k10[i]);
  }

  /* Compute k13 */
  rhseom(ns12, t+h, dt, n, nstate, doubletime, resun, remoon, RITJ, CMAT, SMAT, input, dnp, psatout, it, pmsis, 
	 initCd, msismod, isat, atmcounter, RSMdata, nobs);
  for(i=0; i<6; i++) {
    k13[i] = nstate[i];
    ns13[i] = x[i] + h * (b13_1*k1[i] + b13_4*k4[i] + b13_5*k5[i] + b13_6*k6[i] + b13_7*k7[i] +
			  b13_8*k8[i] + b13_9*k9[i] + b13_10*k10[i] + b13_12*k12[i]);
  }


  for(i=0; i<6; i++) {
    nstate[i] = x[i] + h * (c1*k1[i] + c6*k6[i] + c7*k7[i] + c8*k8[i] + c9*k9[i] + c10*k10[i] + c11*k11[i]);
  }

}






