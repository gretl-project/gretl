/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 2005 Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
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

#ifndef JOHANSEN_H_
#define JOHANSEN_H_

#include "gretl_matrix.h"

typedef enum {
    J_NO_CONST = 0,
    J_REST_CONST,
    J_UNREST_CONST,
    J_REST_TREND,
    J_UNREST_TREND
} JohansenCode;

typedef struct JVAR_ JVAR;

struct JVAR_ {
    JohansenCode code;    /* see above */
    int *list;            /* list of endogenous and exogenous vars */
    int order;            /* order of VAR (order of VECM will be 1 less) */
    int neqns;            /* number of equations = number of endogenous vars */
    int t1;               /* starting observation */
    int t2;               /* ending observation */
    int rank;             /* if specified, chosen cointegration rank, else 0 */
    gretl_matrix *u;      /* resids, VAR in differences */
    gretl_matrix *v;      /* resids, second regressions */
    gretl_matrix *Suu;    /* matrix of cross-products of residuals */
    gretl_matrix *Svv;    /* matrix of cross-products of residuals */
    gretl_matrix *Suv;    /* matrix of cross-products of residuals */
    gretl_matrix **Pi;    /* matrices of coefficients, VAR in differences */
    gretl_matrix **Theta; /* matrices of coefficients, second regressions */
    gretl_matrix *A;      /* matrix of eigenvectors */
    int err;              /* error code */
};

void johansen_VAR_free (JVAR *jv);

JVAR *johansen_test (int order, const int *list, double ***pZ, DATAINFO *pdinfo,
		     gretlopt opt, PRN *prn);

int johansen_test_simple (int order, const int *list, double ***pZ, DATAINFO *pdinfo,
			  gretlopt opt, PRN *prn);
    
#endif /* JOHANSEN_H_ */

