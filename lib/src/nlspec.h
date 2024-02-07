/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

/* Private header for sharing info between nls.c and gmm.c */

#include "libgretl.h" 

typedef struct parm_ parm;
typedef struct ocset_ ocset;

typedef enum {
    NL_ANALYTICAL  = 1 << 0,
    NL_AUTOREG     = 1 << 1,
    NL_AHESS       = 1 << 2,
    NL_NEWTON      = 1 << 3,
    NL_SMALLSTEP   = 1 << 4,
    NL_NAMES_ARRAY = 1 << 5
} nl_flags;

struct nlspec_ {
    int ci;             /* NLS, MLE or GMM */
    int generr;         /* error from genr */
    nl_flags flags;     /* numeric or analytic derivatives, etc. */
    gretlopt opt;       /* can include OPT_V for verbose output; if ci = MLE
			   can also include OPT_H (Hessian) or OPT_R (QML)
			   to control the estimator of the variance matrix;
			   can also include OPT_N to force use of 
			   numerical derivatives.
			*/
    int dv;             /* ID number of dependent variable (NLS) */
    int lhtype;         /* type of the LHS variable */
    char lhname[VNAMELEN]; /* name of LHS var in criterion function */
    char hname[VNAMELEN];  /* name of Hessian matrix, if present */
    char *parnames;     /* user-set names for parameters */
    int lhv;            /* ID number of LHS series in function being
			   minimized or maximized... */
    gretl_matrix *lvec; /* or LHS vector */
    char *nlfunc;       /* string representation of function,
			   expressed in terms of the residuals (NLS,
			   GMM) or the log-likelihood (MLE)
			*/
    int nparam;         /* number of parameters */
    int ncoeff;         /* number of coefficients (allows for vector params) */
    int nmat;           /* number of matrix parameters */
    int naux;           /* number of auxiliary commands */
    int ngenrs;         /* number of variable-generating formulae */
    int iters;          /* number of iterations performed */
    int fncount;        /* number of function evaluations (BFGS) */
    int grcount;        /* number of gradient evaluations (BFGS) */
    int t1;             /* starting observation */
    int t2;             /* ending observation */
    int real_t1;        /* real starting observation (if sub-sampled) */
    int real_t2;        /* real ending observation (if sub-sampled) */
    int nobs;           /* number of observations used */
    double crit;        /* criterion (minimand or maximand) */
    double tol;         /* tolerance for stopping iteration */
    parm *params;       /* array of information on function parameters
			   (see the parm_ struct above) */
    double *fvec;       /* function vector */
    double *coeff;      /* coefficient estimates */
    gretl_matrix *J;    /* Jacobian matrix */
    gretl_matrix *Hinv; /* negative inverse of Hessian */
    char **aux;         /* auxiliary commands */
    char *hesscall;     /* function call for Hessian */
    GENERATOR **genrs;  /* variable-generation pointers */
    GENERATOR *hgen;    /* generator for Hessian */
    DATASET *dset;      /* pointer to dataset */
    PRN *prn;           /* printing aparatus */
    ocset *oc;          /* orthogonality info (GMM) */
    char *missmask;     /* mask for missing observations */
};

void nlspec_destroy_arrays (nlspec *s);

void oc_set_destroy (ocset *oc);

int nl_calculate_fvec (nlspec *s);

int update_coeff_values (const double *x, nlspec *s);

int check_gmm_requirements (nlspec *spec);

int nlspec_add_orthcond (nlspec *s, const char *str,
			 const DATASET *dset);

int nlspec_add_ivreg_oc (nlspec *s, int lhv, const int *rlist,
			 const double **Z);

int nlspec_add_weights (nlspec *s, const char *str);

void nlspec_print_gmm_info (const nlspec *spec, PRN *prn);

void maybe_add_gmm_residual (MODEL *pmod, const nlspec *spec, 
			     const DATASET *dset);

int gmm_add_vcv (MODEL *pmod, nlspec *spec);

int gmm_calculate (nlspec *s, PRN *prn);

int gmm_missval_check_etc (nlspec *s);
