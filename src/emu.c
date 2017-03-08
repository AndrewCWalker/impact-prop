/*
 *  emu.c
 *  
 *  Emulator for GRACE drag prediction.
 *
 *  
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

#include "structs.h"
#include "emu.h"

// Species order will always be H, He, N, N2, O, O2

// Number of simulations, number of inputs, number of chemical species.
static int m=1000, p=7, nspec=6;
// Kriging basis computed by emuInit, sizes should be nspec and m
static double KrigBasis[7][1000];


// Initialization function that computes the Kriging basis
void emuInit(struct RSM_struct *RSMdata) {
    int i,j,k,l;
    double cov;
    gsl_matrix *SigmaSim = gsl_matrix_alloc(m,m);
    gsl_vector *b = gsl_vector_alloc(m);
    
    // Do these one principal component at a time
    for(i=0; i<nspec; i++) {
		
        // Fill in the covariance matrix for the principals components
        // Also make a gsl_vector with the weights.
        for(j=0; j<m; j++) {
            // Diagonal
            gsl_matrix_set(SigmaSim, j, j, (1.0/RSMdata->lamz[i]) + (1.0/RSMdata->lamws[i]));
            // Off-diagonals
            for(k=0; k<j; k++) {
                // Compute the covariance
                cov = 0.0;
                for(l=0; l<p; l++) {
                    cov -= RSMdata->beta[i][l]*pow(RSMdata->x[i][j][l]-RSMdata->x[i][k][l], 2.0);
                }
                cov = exp(cov)/RSMdata->lamz[i];
                gsl_matrix_set(SigmaSim, j, k, cov);
                gsl_matrix_set(SigmaSim, k, j, cov);
            } // for(k=0; k<j; k++)
            gsl_vector_set(b, j, RSMdata->w[i][j]);
        } // for(j=0; j<m; j++)
        
        // Cholesky and solve
        gsl_linalg_cholesky_decomp(SigmaSim);
        gsl_linalg_cholesky_svx(SigmaSim, b);
        
        // Copy into the Kriging Basis
        for(j=0; j<m; j++) {
            KrigBasis[i][j] = gsl_vector_get(b, j);
        }
    } // for(i=0; i<nspec; i++)
    gsl_matrix_free(SigmaSim);
    gsl_vector_free(b);
    
}

// The actual emulation
// input parameters, species weights, output
void emu(struct RSM_struct *RSMdata, double *xstar, double *specw, double *ystar) {
    static int inited=0;
    int i, j, k;
    double xstarstd[nspec][p], wstar[nspec], Sigmastar[nspec][m], logc;

    double mass[nspec];
    double mavg = 0.0;
    /* Define atomic/molecular masses of each species */
    mass[0] = 1.674e-27;  /* hydrogen */
    mass[1] = 6.646e-27;  /* helium */
    mass[2] = 2.326e-26;  /* atomic nitrogen */
    mass[3] = 4.652e-26;  /* diatomic nitrogen */
    mass[4] = 2.657e-26;  /* atomic oxygen */
    mass[5] = 5.314e-26;  /* diatomic oxygen */

    /* Compute average mass */
    /* Compute the average mass */
    for(i=0; i<nspec; i++) {
      mavg += specw[i]*mass[i];
    }
    		
    // Iinitialize if necessary
    if(inited==0) {
        emuInit(RSMdata);
        inited=1;
    }
    
    // Standardize the inputs
    // The different species require slightly different standardizations.
    for(i=0; i<nspec; i++) {
        for(j=0; j<p; j++) {
            xstarstd[i][j] = (xstar[j] - RSMdata->xmin[i][j]) / RSMdata->xrange[i][j];
        }
    }
    
    // Compute the covariances between the new input and sims for all PCs
    for(i=0; i<nspec; i++) {
        for(j=0; j<m; j++) {
            logc = 0.0;
            for(k=0; k<p; k++) {
                logc -= RSMdata->beta[i][k]*pow(RSMdata->x[i][j][k]-xstarstd[i][k], 2.0);
            }
            Sigmastar[i][j] = exp(logc)/RSMdata->lamz[i];
        }
    }
    
    // Compute wstar, the predicted species results for the new input
    for(i=0; i<nspec; i++) {
        wstar[i]=0.0;
        for(j=0; j<m; j++) {
            wstar[i] += Sigmastar[i][j] * KrigBasis[i][j];
        }
		wstar[i] = wstar[i]*RSMdata->sd[i] + RSMdata->mean[i];
    }
    
    // Compute ystar, the new output
    *ystar = 0.0;
    for(j=0; j<nspec; j++) {
      *ystar += specw[j]*wstar[j]*mass[j]/mavg;
    }

}


