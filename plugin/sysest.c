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

#include "libgretl.h"
#include "gretl_matrix.h"
#include "system.h"
#include "sysml.h"

#define SDEBUG 0
#define SYS_USE_SVD 0

/* insert the elements of sub-matrix M, multuplied by scale, in the
   appropriate position within the big matrix X */

static void 
kronecker_place (gretl_matrix *X, const gretl_matrix *M,
		 int startrow, int startcol, double scale)
{
    int i, j;
    int row, col;
    double x;
    
    for (i=0; i<M->rows; i++) {
	row = startrow + i;
	for (j=0; j<M->cols; j++) {
	    col = startcol + j;
	    x = gretl_matrix_get(M, i, j);
	    gretl_matrix_set(X, row, col, x * scale);
	}
    }
}

/* retrieve the special k-class transformed data wanted for LIML
   estimation */

static int make_liml_X_block (gretl_matrix *X, const MODEL *pmod,
			      double **Z, int t1)
{
    int i, t;
    const double *Xi;

    X->cols = pmod->ncoeff;

    for (i=0; i<X->cols; i++) {
	Xi = tsls_get_Xi(pmod, Z, i);
	if (Xi == NULL) {
	    return 1;
	}
	for (t=0; t<X->rows; t++) {
	    gretl_matrix_set(X, t, i, Xi[t+t1]);
	}
    }

    return 0;
}

/* construct the X data block pertaining to a specific equation, using
   either the original data or fitted values from regression on a set
   of instruments */

static int 
make_sys_X_block (gretl_matrix *X, const MODEL *pmod,
		  double **Z, int t1, int method)
{
    int i, t;
    const double *Xi;

    X->cols = pmod->ncoeff;

    for (i=0; i<X->cols; i++) {
	if (method == SYS_METHOD_3SLS || 
	    method == SYS_METHOD_FIML || 
	    method == SYS_METHOD_TSLS) {
	    Xi = tsls_get_Xi(pmod, Z, i);
	} else {
	    Xi = Z[pmod->list[i+2]];
	}
	if (Xi == NULL) {
	    return E_DATA;
	}
	for (t=0; t<X->rows; t++) {
	    gretl_matrix_set(X, t, i, Xi[t+t1]);
	}
    }

    return 0;
}

/* populate the cross-equation covariance matrix based on the
   per-equation residuals */

static int
gls_sigma_from_uhat (equation_system *sys, gretl_matrix *sigma)
{
    const gretl_matrix *e = sys->uhat;
    int m = sys->neqns;
    int T = sys->T;
    int geomean = system_vcv_geomean(sys);
    int i, j, t;
    double xx;

    for (i=0; i<m; i++) {
	for (j=i; j<m; j++) {
	    xx = 0.0;
	    for (t=0; t<T; t++) {
		xx += gretl_matrix_get(e, t, i) * gretl_matrix_get(e, t, j);
	    }
	    if (geomean) {
		xx /= system_vcv_denom(sys, i, j);
	    } else {
		xx /= T;
	    }
	    gretl_matrix_set(sigma, i, j, xx);
	    if (j != i) {
		gretl_matrix_set(sigma, j, i, xx);
	    }
	}
    }

    if (sys->method == SYS_METHOD_OLS && sys->diag == 0.0) {
	double sii, sij, sjj;

	for (i=1; i<m; i++) {
	    sii = gretl_matrix_get(sigma, i, i);
	    for (j=0; j<i; j++) {
		sij = gretl_matrix_get(sigma, i, j);
		sjj = gretl_matrix_get(sigma, j, j);
		sys->diag += (sij * sij) / (sii * sjj);
	    }
	}
	sys->diag *= T;
    }

    return 0;
}

/* compute residuals, for all cases other than FIML */

static void 
sys_resids (equation_system *sys, int eq, const double **Z)
{
    MODEL *pmod = sys->models[eq];
    double yh;
    int i, t;

    pmod->ess = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	yh = 0.0;
	for (i=0; i<pmod->ncoeff; i++) {
	    yh += pmod->coeff[i] * Z[pmod->list[i+2]][t];
	}
	pmod->yhat[t] = yh;
	pmod->uhat[t] = Z[pmod->list[1]][t] - yh;
	/* for cross-equation vcv */
	gretl_matrix_set(sys->uhat, t - pmod->t1, pmod->ID, pmod->uhat[t]);
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }

    /* df correction? */
    if (system_want_df_corr(sys)) {
	pmod->sigma = sqrt(pmod->ess / pmod->dfd);
    } else {
	pmod->sigma = sqrt(pmod->ess / pmod->nobs);
    }
}

static void
liml_scale_vcv (equation_system *sys, gretl_matrix *vcv)
{
    double s2, vij;
    int vmin = 0;
    int vi, vj;
    int i, j, k;

    for (i=0; i<sys->neqns; i++) {
	s2 = sys->models[i]->sigma * sys->models[i]->sigma;
	for (j=0; j<sys->models[i]->ncoeff; j++) {
	    for (k=j; k<sys->models[i]->ncoeff; k++) {
		vi = j + vmin;
		vj = k + vmin;
		vij = gretl_matrix_get(vcv, vi, vj);
		vij *= s2;
		gretl_matrix_set(vcv, vi, vj, vij);
		gretl_matrix_set(vcv, vj, vi, vij);
	    }
	}
	vmin += sys->models[i]->ncoeff;
    }		
}

/* calculate the standard error of the residuals for the system as a
   whole and use this for scaling the covariance matrix */

static void 
single_eq_common_sigma_scale_vcv (equation_system *sys, gretl_matrix *vcv)
{
    double s, ess = 0.0;
    int nr = 0, dfc = 0;
    int den = 0;
    int i;

    if (system_want_df_corr(sys)) {
	nr = system_n_restrictions(sys);
	dfc = 1;
    }

    /* is this right? */

    for (i=0; i<sys->neqns; i++) {
	ess += sys->models[i]->ess;
	den += sys->models[i]->nobs;
	if (dfc) {
	    den -= sys->models[i]->ncoeff;
	}
    }

    den += nr;

    s = sqrt(ess / den);

    gretl_matrix_multiply_by_scalar(vcv, s * s);
}

/* use the per-equation residual standard deviations for scaling
   the covariance matrix, block by diagonal block
*/

static void 
single_eq_scale_vcv (equation_system *sys, gretl_matrix *vcv)
{
    MODEL *pmod;
    double s2, vij;
    int ioff = 0, joff = 0;
    int i, k, ii, jj;

    for (i=0; i<sys->neqns; i++) {
	pmod = sys->models[i];
	s2 = pmod->sigma * pmod->sigma;
	k = pmod->ncoeff;
	for (ii=0; ii<k; ii++) {
	    for (jj=0; jj<k; jj++) {
		vij = gretl_matrix_get(vcv, ioff + ii, joff + jj);
		gretl_matrix_set(vcv, ioff + ii, joff + jj, vij * s2);
	    }
	}
	ioff += k;
	joff += k;
    }    
}

/* compute SUR, 3SLS or LIML parameter estimates (or restricted OLS,
   TSLS, WLS) */

static int 
calculate_sys_coeffs (equation_system *sys,
		      const double **Z,
		      gretl_matrix *X, gretl_matrix *y, 
		      int mk, int nr, int do_iteration)
{
    int do_bdiff = ((sys->method == SYS_METHOD_3SLS) && do_iteration);
    double bij, oldb, bnum = 0.0, bden = 0.0;
    gretl_matrix *vcv = NULL;
    gretl_matrix *b = NULL;
    int i, j, k, j0;
    int err = 0;

    vcv = gretl_matrix_copy(X);
    if (vcv == NULL) {
	return 1;
    }

    err = gretl_LU_solve(X, y);
    if (err) {
	return err;
    }

#if SDEBUG > 1
    gretl_matrix_print(y, "in calc_coeffs, betahat");
#endif

#if SYS_USE_SVD
    /* memory-intensive for big matrices */
    err = gretl_SVD_invert_matrix(vcv);
#else    
    err = gretl_invert_general_matrix(vcv);
#endif

#if SDEBUG
    fprintf(stderr, "calculate_sys_coeffs: invert, err=%d\n", err);
#endif

    if (err) {
	return err;
    }

#if SDEBUG > 1
    gretl_matrix_print(vcv, "in calc_coeffs, vcv");
#endif

    j0 = 0;
    for (i=0; i<sys->neqns; i++) {
	for (j=0; j<sys->models[i]->ncoeff; j++) {
	    k = j0 + j;
	    bij = gretl_vector_get(y, k);
	    if (do_bdiff) {
		oldb = sys->models[i]->coeff[j];
		bnum += (bij - oldb) * (bij - oldb);
		bden += oldb * oldb;
	    }
	    sys->models[i]->coeff[j] = bij;
	}
	sys_resids(sys, i, Z);
	j0 += sys->models[i]->ncoeff;
    }

    if (do_bdiff) {
	sys->bdiff = sqrt(bnum / bden);
    }

    /* simple single-equation methods: need to multiply by an estimate
       of sigma.  Should this really be the system sigma, and not
       equation-specific?
    */

    if (sys->method == SYS_METHOD_OLS || 
	sys->method == SYS_METHOD_TSLS) {
	if (nr > 0) {
	    single_eq_common_sigma_scale_vcv(sys, vcv);
	} else {
	    single_eq_scale_vcv(sys, vcv);
	}
    } else if (sys->method == SYS_METHOD_LIML) {
	liml_scale_vcv(sys, vcv);
    }

    /* now set the model standard errors */
    j0 = 0;
    for (i=0; i<sys->neqns; i++) {
	for (j=0; j<sys->models[i]->ncoeff; j++) {
	    k = j0 + j;
	    sys->models[i]->sderr[j] = sqrt(gretl_matrix_get(vcv, k, k));
	}
	j0 += sys->models[i]->ncoeff;
    }

    /* save the coefficient vector and covariance matrix */
    b = gretl_matrix_copy(y);
    system_attach_coeffs(sys, b);
    system_attach_vcv(sys, vcv);

    return err;
}

static int in_list (const int *list, int k)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (k == list[i]) return 1;
    }

    return 0;
}

static int *
system_model_list (equation_system *sys, int i, int *freeit)
{
    int *list = NULL;

    *freeit = 0;

    if (sys->method == SYS_METHOD_SUR || 
	sys->method == SYS_METHOD_3SLS ||
	sys->method == SYS_METHOD_OLS || 
	sys->method == SYS_METHOD_TSLS ||
	sys->method == SYS_METHOD_WLS) {
	list = system_get_list(sys, i);
    }

    if (sys->method == SYS_METHOD_3SLS || 
	sys->method == SYS_METHOD_TSLS) {
	/* is list already in tsls form? */
	if (list != NULL && !in_list(list, LISTSEP)) {
	    list = NULL;
	}
    }

    if (sys->method == SYS_METHOD_FIML || 
	sys->method == SYS_METHOD_LIML ||
	((sys->method == SYS_METHOD_3SLS || 
	  sys->method == SYS_METHOD_TSLS) && list == NULL)) {
	list = compose_tsls_list(sys, i);
	*freeit = 1;
    }

    return list;
}

/* Hansen-Sargan overidentification test for the system as a whole, as
   in Davidson and MacKinnon, ETM: p. 511 and equation (12.25) for the
   case of SUR; p. 532 and equation (12.61) for the case of 3SLS.  See
   also D & M, Estimation and Inference in Econometrics, equation
   (18.60), for a more computation-friendly statement of the criterion
   function.
*/

static int hansen_sargan_test (equation_system *sys, 
			       const double **Z)
{
    const int *exlist = system_get_instr_vars(sys);
    int nx = exlist[0];
    int m = sys->neqns;
    int T = sys->T;
    int df = system_get_overid_df(sys);

    const double *Wi, *Wj;

    gretl_matrix *WTW = NULL;
    gretl_matrix *eW = NULL;
    gretl_matrix *tmp = NULL;

    double x, X2;
    int i, j, t;
    int err = 0;

    if (df <= 0) return 1;

    WTW = gretl_matrix_alloc(nx, nx);
    eW = gretl_matrix_alloc(m, nx);
    tmp = gretl_matrix_alloc(m, nx);

    if (WTW == NULL || eW == NULL || tmp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* construct W-transpose W */
    for (i=0; i<nx; i++) {
	Wi = Z[exlist[i+1]] + sys->t1;
	for (j=i; j<nx; j++) {
	    Wj = Z[exlist[j+1]] + sys->t1;
	    x = 0.0;
	    for (t=0; t<sys->T; t++) {
		x += Wi[t] * Wj[t];
	    }
	    gretl_matrix_set(WTW, i, j, x);
	    if (i != j) {
		gretl_matrix_set(WTW, j, i, x);
	    }
	}
    }

    err = gretl_invert_symmetric_matrix(WTW);
#if SDEBUG
    fprintf(stderr, "hansen_sargan: on invert, err=%d\n", err);
#endif
    if (err) {
	sys->X2 = NADBL;
	goto bailout;
    }

    /* set up vectors of SUR or 3SLS residuals, transposed, times W:
       these are stacked in an m * nx matrix
    */
    for (i=0; i<m; i++) {
	for (j=0; j<nx; j++) {
	    Wj = Z[exlist[j+1]] + sys->t1;
	    x = 0.0;
	    for (t=0; t<T; t++) {
		x += gretl_matrix_get(sys->uhat, t, i) * Wj[t];
	    }
	    gretl_matrix_set(eW, i, j, x);
	}
    }

    /* multiply these vectors into (WTW)^{-1} */
    for (i=0; i<m; i++) {
	for (j=0; j<nx; j++) {
	    x = 0.0;
	    for (t=0; t<nx; t++) {
		x += gretl_matrix_get(eW, i, t) * 
		    gretl_matrix_get(WTW, t, j);
	    }
	    gretl_matrix_set(tmp, i, j, x);
	}
    }   

    /* cumulate the Chi-square value */
    X2 = 0.0;
    for (i=0; i<m; i++) {
	for (j=0; j<m; j++) {
	    x = 0.0;
	    for (t=0; t<nx; t++) {
		x += gretl_matrix_get(tmp, i, t) * 
		    gretl_matrix_get(eW, j, t); /* transposed */
	    }
	    X2 += gretl_matrix_get(sys->sigma, i, j) * x;
	}
    }

#if SDEBUG
    fprintf(stderr, "Hansen-Sargan: Chi-square(%d) = %g (p-value %g)\n", 
	    df, X2, chisq_cdf_comp(X2, df));
#endif
    sys->X2 = X2;

 bailout:

    gretl_matrix_free(WTW);
    gretl_matrix_free(eW);
    gretl_matrix_free(tmp);

    return err;
}

static int basic_system_allocate (equation_system *sys,
				  int mk, int nr, 
				  gretl_matrix **X,
				  gretl_matrix **y)
{
    int m = sys->neqns;
    int T = sys->T;
    int ldx = mk + nr;

    /* allocate a model for each stochastic equation */
    sys->models = gretl_model_array_new(m);

    sys->uhat = gretl_matrix_alloc(T, m);
    if (sys->uhat == NULL) {
	return E_ALLOC;
    }

    sys->sigma = gretl_matrix_alloc(m, m);
    if (sys->sigma == NULL) {
	return E_ALLOC;
    }

    *X = gretl_matrix_alloc(ldx, ldx);
    if (*X == NULL) {
	return E_ALLOC;
    }

    *y = gretl_column_vector_alloc(ldx);
    if (*y == NULL) {
	return E_ALLOC;
    }

    return 0;
}

/* compute log-likelihood for iterated SUR estimator */

double sur_ll (equation_system *sys)
{
    int m = sys->neqns;
    int T = sys->T;
    gretl_matrix *sigtmp;
    double ldet;

    sigtmp = gretl_matrix_alloc(m, m);
    if (sigtmp == NULL) return NADBL;

    gls_sigma_from_uhat(sys, sigtmp);
    ldet = gretl_vcv_log_determinant(sigtmp);

    if (na(ldet)) {
	sys->ll = NADBL;
    } else {
	sys->ll = -(m * T / 2.0) * (LN_2_PI + 1.0);
	sys->ll -= (T / 2.0) * ldet;
    }

    gretl_matrix_free(sigtmp);

    return sys->ll;
}

/* if we're estimating with a specified set of linear restrictions,
   Rb = q, augment the X matrix with R and R-transpose */

static int 
augment_X_with_restrictions (gretl_matrix *X, int mk, 
			     equation_system *sys)
{
    double rij;
    int nr, nc;
    int i, j;

    if (sys->R == NULL) return 1;

    nr = sys->R->rows;
    nc = sys->R->cols;

    /* place the R matrix */
    kronecker_place(X, sys->R, mk, 0, 1.0);

    /* place R-transpose */
    for (i=0; i<nr; i++) {
	for (j=0; j<nc; j++) {
	    rij = gretl_matrix_get(sys->R, i, j);
	    gretl_matrix_set(X, j, i + mk, rij);
	}
    }

    /* zero the bottom right-hand block */
    for (i=mk; i<mk+nr; i++) {
	for (j=mk; j<mk+nr; j++) {
	    gretl_matrix_set(X, i, j, 0.0);
	}
    }

    return 0;
}

/* when estimating with a specified set of linear restrictions,
   Rb = q, augment the y matrix with q */

static int 
augment_y_with_restrictions (gretl_matrix *y, int mk, int nr,
			     equation_system *sys)
{
    int i;

    if (sys->q == NULL) return 1;

    for (i=0; i<nr; i++) {
	gretl_vector_set(y, mk + i, gretl_vector_get(sys->q, i));
    }

    return 0;
}

#define SYS_MAX_ITER 100
#define SYS_LL_TOL 1.0e-12
#define SYS_BDIFF_TOL 1.0e-9

/* check for convergence of iteration: we use the change in the
   log-likelihood when iterating SUR to the ML solution, or
   a measure of the change in the coefficients when iterating
   three-stage least squares
*/

static int converged (equation_system *sys, 
		      double *llbak, int *err,
		      PRN *prn)
{
    double crit, tol = 0.0;
    int met = 0;

    if (sys->method == SYS_METHOD_SUR || 
	sys->method == SYS_METHOD_WLS) {
	double ll = sur_ll(sys);

	tol = SYS_LL_TOL;
	crit = ll - *llbak;

#if SDEBUG
	printf("SUR iteration %d, ll = %.8g\n", sys->iters, ll);
#endif
	if (crit <= tol) {
	    met = 1;
	} else if (sys->iters < SYS_MAX_ITER) {
	    *llbak = ll;
	} 
    } else if (sys->method == SYS_METHOD_3SLS) {
	tol = SYS_BDIFF_TOL;
	crit = sys->bdiff;

#if SDEBUG
	printf("3SLS iteration %d, crit = %.8g\n", sys->iters, crit);
#endif
	if (crit <= tol) {
	    met = 1;
	} 
    }

    if (!met && sys->iters >= SYS_MAX_ITER) {
	pprintf(prn, "reached %d iterations without meeting "
		"tolerance of %g\n", sys->iters, tol);
	*err = E_NOCONV;
    }	

    return met;
}

static void clean_up_models (equation_system *sys)
{
    double ess = 0.0;
    int i;

    for (i=0; i<sys->neqns; i++) {
	ess += sys->models[i]->ess;
	if (sys->method == SYS_METHOD_3SLS || 
	    sys->method == SYS_METHOD_FIML || 
	    sys->method == SYS_METHOD_TSLS || 
	    sys->method == SYS_METHOD_LIML) {
	    tsls_free_data(sys->models[i]);
	}
	gretl_model_free(sys->models[i]);
    }

    free(sys->models);
    sys->models = NULL;

    sys->ess = ess;
}

static int drop_redundant_instruments (equation_system *sys,
				       const int *droplist, int i)
{
    int j, k, pos, err = 0;

    for (j=1; j<=droplist[0]; j++) {
	pos = gretl_list_position(droplist[j], sys->ilist);
	if (pos > 0) {
	    gretl_list_delete_at_pos(sys->ilist, pos);
	} else {
	    err = 1;
	}
    }

    pos = gretl_list_separator_position(sys->lists[i]);
    if (pos > 0) {
	for (j=1; j<=droplist[0]; j++) {
	    for (k=pos+1; k<=sys->lists[i][0]; k++) {
		if (sys->lists[i][k] == droplist[j]) {
		    gretl_list_delete_at_pos(sys->lists[i], k);
		    break;
		}
	    }
	}
    }

    return err;
}

/* general function that forms the basis for all specific system
   estimators */

int system_estimate (equation_system *sys, double ***pZ, DATAINFO *pdinfo, 
		     gretlopt opt, PRN *prn)
{
    int i, j, k, T, t, m = 0;
    int v, l, mk, krow, nr;
    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;

    gretl_matrix *X = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *Xi = NULL;
    gretl_matrix *Xj = NULL;
    gretl_matrix *M = NULL;
    MODEL **models = NULL;

    int method = sys->method;
    double llbak = -1.0e9;
    int single_equation = 0;
    int do_iteration = 0;
    int rtsls = 0;
    int err = 0;

    sys->iters = 0;

    if (sys->flags & SYSTEM_ITERATE) {
	do_iteration = 1;
    }
    
    nr = system_n_restrictions(sys);

    if (method == SYS_METHOD_OLS || method == SYS_METHOD_TSLS || 
	method == SYS_METHOD_LIML || method == SYS_METHOD_WLS) {
	single_equation = 1;
    }

    if (nr > 0 && method == SYS_METHOD_3SLS) {
	/* doing 3SLS with restrictions: we want to obtain
	   restricted TSLS estimates as a starting point */
	rtsls = 1;
    }

    /* get uniform sample starting and ending points */
    if (system_adjust_t1t2(sys, &pdinfo->t1, &pdinfo->t2, 
			   (const double **) *pZ)) {
	err = E_DATA;
	goto cleanup;
    } 

    /* number of equations */
    m = sys->neqns;

    /* max indep vars per equation */
    k = system_max_indep_vars(sys);

    /* total indep vars, all equations */
    mk = system_n_indep_vars(sys);

    /* number of observations per series */
    T = sys->T;

    /* allocate models etc */
    err = basic_system_allocate(sys, mk, nr, &X, &y);
    if (err) goto cleanup;

    /* convenience pointers */
    models = sys->models;
	    
    if ((method == SYS_METHOD_FIML || 
	 method == SYS_METHOD_LIML) && !(opt & OPT_Q)) {
	print_equation_system_info(sys, pdinfo, OPT_H, prn);
    }

    /* First estimate the equations separately (either by OLS or
       TSLS), and put the single-equation residuals into the uhat
       matrix.  Note that at this stage we are not in a position to
       impose any cross-equation restrictions, since we're doing
       straight equation-by-equation estimation.
    */

    for (i=0; i<m; i++) {
	int freeit = 0;
	int *list = system_model_list(sys, i, &freeit);
	const int *droplist = NULL;
	MODEL *pmod = models[i];

	if (list == NULL) {
	    err = 1;
	    break;
	}

	if (method == SYS_METHOD_SUR || 
	    method == SYS_METHOD_OLS || 
	    method == SYS_METHOD_WLS) {
	    *pmod = lsq(list, pZ, pdinfo, OLS, OPT_A);
	} else if (method == SYS_METHOD_3SLS || 
		   method == SYS_METHOD_FIML || 
		   method == SYS_METHOD_TSLS) {
	    *pmod = tsls_func(list, SYSTEM, pZ, pdinfo, OPT_NONE);
	} else if (method == SYS_METHOD_LIML) {
	    *pmod = tsls_func(list, SYSTEM, pZ, pdinfo, OPT_N);
	}

	if (freeit) {
	    free(list);
	}

	if ((err = pmod->errcode)) {
	    fprintf(stderr, "system_estimate: failed to estimate equation %d: "
		    "err = %d\n", i+1, err);
	    break;
	} 

	droplist = gretl_model_get_data(pmod, "inst_droplist");
	if (droplist != NULL) {
	    drop_redundant_instruments(sys, droplist, i);
	}

	pmod->ID = i;
	pmod->aux = AUX_SYS;
	gretl_model_set_int(pmod, "method", method);

	/* save the sigma-squared for an LR test for a diagonal
	   covariance matrix */
	if (method == SYS_METHOD_SUR && do_iteration && nr == 0) {
	    gretl_model_set_double(pmod, "ols_sigma_squared",
				   pmod->ess / pmod->nobs);
	}

	for (t=0; t<T; t++) {
	    gretl_matrix_set(sys->uhat, t, i, pmod->uhat[t + sys->t1]);
	}
    }

    if (err) {
	fprintf(stderr, "system_estimate: after single-equation "
		"estimation, err = %d\n", err);
	goto cleanup;
    }

    if (method == SYS_METHOD_LIML) {
	/* compute the minimum eigenvalues and generate the
	   suitably transformed data matrices */
	err = liml_driver(sys, pZ, pdinfo, prn);
	if (err) goto cleanup;
    }

    /* marker for iterated versions of SUR, WLS, or 3SLS; also for
       loopback in case of restricted 3SLS, where we want to compute
       restricted TSLS estimates first
    */
 iteration_start:

    gls_sigma_from_uhat(sys, sys->sigma);

#if SDEBUG > 1
    gretl_matrix_print(sys->sigma, "gls_sigma_from_uhat");
#endif

    if (method == SYS_METHOD_WLS) {
	gretl_matrix_zero(X);
	err = gretl_invert_diagonal_matrix(sys->sigma);
    } else if (single_equation || rtsls) {
	gretl_matrix_zero(X);
    } else {
	err = gretl_invert_symmetric_matrix(sys->sigma);    
    }

#if SDEBUG
    fprintf(stderr, "system_estimate: on invert, err=%d\n", err);
#endif

    if (err) goto cleanup;

    /* the initial tests against NULL here allow for the possibility
       that we're iterating */

    if (Xi == NULL) {
	Xi = gretl_matrix_alloc(T, k);
    }
    if (Xj == NULL) {
	Xj = gretl_matrix_alloc(T, k);
    } 
    if (M == NULL) {
	M = gretl_matrix_alloc(k, k);
    }

    if (Xi == NULL || Xj == NULL || M == NULL) {
	err = E_ALLOC;
	goto cleanup;
    }

    /* form the big stacked X matrix: Xi = data matrix for equation i, 
       specified in lists[i] 
    */

    krow = 0;
    for (i=0; i<m && !err; i++) {
	int kcol = 0;

	err = make_sys_X_block(Xi, models[i], *pZ, sys->t1, method);

	for (j=0; j<m && !err; j++) { 
	    const gretl_matrix *Xk;
	    double sij;

	    if (i != j) {
		if (single_equation || rtsls) {
		    kcol += models[j]->ncoeff;
		    continue;
		}
		err = make_sys_X_block(Xj, models[j], *pZ, sys->t1, method);
		Xk = Xj;
	    } else if (method == SYS_METHOD_LIML) {
		err = make_liml_X_block(Xj, models[i], *pZ, sys->t1);
		Xk = Xj;
	    } else {
		Xk = Xi;
	    }

	    M->rows = Xi->cols;
	    M->cols = Xk->cols;

	    err = gretl_matrix_multiply_mod(Xi, GRETL_MOD_TRANSPOSE,
					    Xk, GRETL_MOD_NONE, 
					    M, GRETL_MOD_NONE);

	    if (rtsls || (single_equation && method != SYS_METHOD_WLS)) {
		sij = 1.0;
	    } else {
		sij = gretl_matrix_get(sys->sigma, i, j);
	    }

	    kronecker_place(X, M, krow, kcol, sij);
	    kcol += models[j]->ncoeff;
	}

	krow += models[i]->ncoeff;
    }

    if (err) {
	fprintf(stderr, "after trying to make X matrix: err = %d\n", err);
	goto cleanup;
    }

    if (nr > 0) {
	/* there are restrictions to be imposed */
	augment_X_with_restrictions(X, mk, sys);
    }

    if (!do_iteration && !rtsls) {
	/* we're not coming back this way, so free some storage */
	gretl_matrix_free(Xj);
	Xj = NULL;
	gretl_matrix_free(M);
	M = NULL;
    }

    /* form stacked Y column vector (m x k) */

    v = 0;
    for (i=0; i<m; i++) { 

	/* loop over the m vertically-arranged blocks */

	make_sys_X_block(Xi, models[i], *pZ, sys->t1, method);

	for (j=0; j<models[i]->ncoeff; j++) { 
	    /* loop over the rows within each of the m blocks */
	    double yv = 0.0;
	    int lmin = 0, lmax = m;

	    if (single_equation || rtsls) {
		/* no cross terms wanted */
		lmin = i;
		lmax = i + 1;
	    }

	    for (l=lmin; l<lmax; l++) { 
		/* loop over the components that must be
		   added to form each element */
		const double *yl = NULL;
		double sil, xx = 0.0;

		if (method == SYS_METHOD_LIML) {
		    yl = gretl_model_get_data(models[l], "liml_y");
		} else {
		    yl = (*pZ)[system_get_depvar(sys, l)];
		}

		/* multiply X'[l] into y */
		for (t=0; t<T; t++) {
		    xx += gretl_matrix_get(Xi, t, j) * yl[t + sys->t1];
		}

		if (rtsls || (single_equation && method != SYS_METHOD_WLS)) {
		    sil = 1.0;
		} else {
		    sil = gretl_matrix_get(sys->sigma, i, l);
		}

		yv += xx * sil;
	    }
	    gretl_vector_set(y, v++, yv);
	}
    }

    if (nr > 0) {
	/* there are restrictions */
	augment_y_with_restrictions(y, mk, nr, sys);
    }

#if SDEBUG > 1
    gretl_matrix_print(X, "sys X");
    gretl_matrix_print(y, "sys y");
#endif    

    /* The estimates calculated below will be SUR, 3SLS or LIML,
       depending on how the data matrices above were constructed --
       unless, that is, we're just doing restricted OLS, WLS or TSLS
       estimates.
    */
    calculate_sys_coeffs(sys, (const double **) *pZ, X, 
			 y, mk, nr, do_iteration);

    if (rtsls) {
	rtsls = 0;
	goto iteration_start;
    }

    if (do_iteration) {
	if (!converged(sys, &llbak, &err, prn)) {
	    if (err) {
		goto cleanup;
	    } else {
		sys->iters += 1;
		goto iteration_start;
	    }
	}
    }

    if (nr == 0 && (method == SYS_METHOD_3SLS || method == SYS_METHOD_SUR)) {
	/* compute this test while we have sigma-inverse available */
	hansen_sargan_test(sys, (const double **) *pZ);
    }

    /* refresh sigma (non-inverted) */
    gls_sigma_from_uhat(sys, sys->sigma);

    if (method == SYS_METHOD_FIML) {
	/* compute FIML estimates */
	err = fiml_driver(sys, pZ, pdinfo, opt, prn);
    }

    if (!err) {
	err = system_save_and_print_results(sys, pZ, pdinfo, 
					    opt, prn);
    }

 cleanup:

    gretl_matrix_free(Xi);
    gretl_matrix_free(Xj);
    gretl_matrix_free(M);
    gretl_matrix_free(X);
    gretl_matrix_free(y);

    if (models != NULL) {
	clean_up_models(sys);
    }

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    return err;
}
