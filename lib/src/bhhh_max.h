/* 
 * Copyright (C) 2004 Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

enum {
    PRESERVE_OPG_MODEL = 1 << 0,
    FULL_VCV_MATRIX    = 1 << 1
} bhhh_opts;

typedef struct _model_info model_info;

struct _model_info {

    /* members that should be set by caller of bhhh_max: */

    int k;              /* number of parameters */
    int p, q, r;        /* for use with ARMA: AR, MA orders and number of
                           other regressors */
    int t1, t2;         /* starting and ending point of sample */
    int n_series;       /* number of additional series needed in the
                           likelihood and/or score calculations */
    double tol;         /* tolerance for convergence */
    unsigned char opts; /* options from among bhhh_opts */

    /* members set within bhhh_max: */

    int n;            /* length of series */
    double ll;        /* log-likelihood */
    double ll2;       /* temporary log-likelihood value */
    double s2;        /* error variance */
    int *list;        /* OPG regression list */
    double *theta;    /* vector of parameters */
    double **series;  /* additional series */

    /* full VCV matrix from OPG regression, if (opts & FULL_VCV_MATRIX) */
    gretl_matrix *VCV; 

    /* pointer to OPG model, if (opts & PRESERVE_OPG_MODEL) */
    MODEL *pmod;
};

void model_info_free (model_info *model);

int bhhh_max (int (*loglik) (double *, const double **, double **,
			     model_info *, int), 
	      const double **X, const double *init_coeff,
	      model_info *model, PRN *prn);
