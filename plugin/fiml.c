/*
 *  Copyright (c) 2004 by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h"
#include "gretl_matrix.h"
#include "gretl_matrix_private.h"

typedef struct fiml_system_ fiml_system;

struct fiml_system_ {
    int n;                  /* number of observations per equation */
    int g;                  /* number of (stochastic) equations */
    int gn;                 /* convenience: g * n = number of obs in stacked vectors */
    gretl_equation_system *sys; /* pointer to "parent" equation system */
    gretl_vector *udot;     /* stacked gn-vector: LHS of artificial regression */
    gretl_matrix *xdot;     /* RHS: transformed stacked matrix of indep vars */
    gretl_matrix *uhat;     /* structural-form residuals, all equations */
    gretl_matrix *sigma;    /* cross-equation covariance matrix */
};

static int 
fill_fiml_udot (fiml_system *fsys, const double **Z, int t1)
{
    int i, t;

    /* stack the g dependent variables vertically */
    for (i=0; i<fsys->g; i++) {
	int k = system_get_depvar(fsys->sys, i);

	for (t=0; t<fsys->n; t++) {
	    int gt = i * fsys->n + t;

	    gretl_vector_set(fsys->udot, gt, Z[k][t + t1]);
	}
    }

    /* subtract ML fitted values to get ML residuals */
    

    return 0;
}

static void fiml_system_destroy (fiml_system *fsys)
{
    if (fsys->udot != NULL) {
	gretl_vector_free(fsys->udot);
    }

    if (fsys->xdot != NULL) {
	gretl_matrix_free(fsys->xdot);
    }    

    if (fsys->uhat != NULL) {
	gretl_matrix_free(fsys->uhat);
    } 

    if (fsys->sigma != NULL) {
	gretl_matrix_free(fsys->sigma);
    }     

    free(fsys);
}

static fiml_system *fiml_system_new (gretl_equation_system *sys)
{
    fiml_system *fsys;

    fsys = malloc(sizeof *fsys);
    if (fsys == NULL) return NULL;

    fsys->sys = sys;

    fsys->g = system_n_equations(sys);
    fsys->n = system_n_obs(sys);
    fsys->gn = fsys->g * fsys->n;

    fsys->udot = NULL;
    fsys->xdot = NULL;
    fsys->uhat = NULL;
    fsys->sigma = NULL;

    fsys->udot = gretl_column_vector_alloc(fsys->gn);
    fsys->xdot = gretl_matrix_alloc(fsys->gn, fsys->gn); /* FIXME dimensions */
    fsys->uhat = gretl_matrix_alloc(fsys->n, fsys->g);
    fsys->sigma = gretl_matrix_alloc(fsys->g, fsys->g);

    if (fsys->udot == NULL || fsys->sdot == NULL || 
	fsys->uhat == NULL || fsys->sigma == NULL) {
	fiml_system_destroy(fsys);
	fsys = NULL;
    }

    return fsys;
}

static int form_fiml_psi (fiml_system *fsys)
{
    int err = 0;

    err = gretl_matrix_multiply_mod(fsys->uhat, GRETL_MOD_TRANSPOSE,
				    fsys->uhat, GRETL_MOD_NONE,
				    fsys->sigma);

    if (!err) {
	err = gretl_invert_symmetric_matrix(fsys->sigma);
    }

    if (!err) {
	err = gretl_matrix_cholesky_decomp(fsys->sigma);
	gretl_square_matrix_transpose(fsys->sigma);
	gretl_matrix_zero_lower(fsys->sigma);
    }

    return err;
}

static int fiml_transcribe_results (fiml_system *fsys)
{
    int err = 0;

    /* to be written */

    return err;
}

/* Below: this is meant to grow into a driver function for FIML as
   described in Davidson and McKinnon, ETM, chap 12, sect 5.
*/

int fiml_driver (gretl_equation_system *sys, const double **Z, const DATAINFO *pdinfo)
{
    fiml_system *fsys;
    int err = 0;

    fsys = fiml_system_new(sys);
    if (fsys == NULL) {
	return E_ALLOC;
    }

    /* work in progress!! */

    err = fill_fiml_udot(fsys, Z, pdinfo->t1); 

    /* initialize uhat based on 3SLS coefficient estimates */

    /* form \hat{Sigma} (ETM, equation 12.81), invert, 
       and Cholesky-decompose to get \Psi */

    /* form "udot" and "xdot" */

    /* run artificial regression (ETM, equation 12.86) */
    gretl_matrix_ols(fsys->udot, fsys->xdot, fsys->b, NULL);

    /* need a loop here, with appropriate convergence/fail conditions */

    /* write the results into the parent system */
    fiml_transcribe_results(fsys);

    
    /* clean up */
    fim_system_destroy(fsys);

    return err;
}



    
	

    
    
    

    
    
