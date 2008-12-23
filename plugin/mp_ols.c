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

/* mp_ols.c - gretl least squares with multiple precision (GMP) */

#include "libgretl.h"
#include "libset.h"

#include <float.h>
#include <gmp.h>

#ifdef HAVE_MPFR
#include <mpfr.h>
#endif

#define MP_DEBUG 0

static mpf_t MPF_ONE;
static mpf_t MPF_ZERO;
static mpf_t MPF_MINUS_ONE;
static mpf_t MPF_TINY;

typedef struct {
    int ID;                      /* ID number for model */
    int t1, t2, nobs;            /* starting observation, ending
                                    observation, and number of obs */
    int ncoeff, dfn, dfd;        /* number of coefficents; degrees of
                                    freedom in numerator and denominator */
    int *list;                   /* list of variables by (translated) ID number */
    int *varlist;                /* counterpart of list, using "real" ID #s */
    const int *polylist;         /* list of polynomial powers */
    int ifc;                     /* = 1 if the equation includes a constant,
				    else = 0 */
    mpf_t *coeff;                /* array of coefficient estimates */
    mpf_t *sderr;                /* array of estimated std. errors */
    mpf_t *xpx;                  /* X'X */
    mpf_t ess, tss;              /* Error and Total Sums of Squares */
    mpf_t sigma;                 /* Standard error of regression */
    mpf_t rsq, adjrsq;           /* Unadjusted and adjusted R^2 */     
    mpf_t fstt;                  /* F-statistic */
    int errcode;                 /* Error code in case of failure */
    int polyvar;                 /* number of the variable to be raised to
				    specified powers, if any */
} MPMODEL;

typedef struct {
    mpf_t *xpx;
    mpf_t *xpy;
    int ivalue;
    int nv;
    int errcode;
} MPXPXXPY;	

typedef struct {
    MPXPXXPY xpxxpy;
    mpf_t *coeff;
    mpf_t rss;
    int errcode;
} MPCHOLBETA;

static void set_gretl_mp_bits (void);
#ifdef HAVE_MPFR
static void set_gretl_mpfr_bits (void);
#endif

static void mpf_constants_init (void)
{
    mpf_init_set_d(MPF_ONE, 1.0);
    mpf_init_set_d(MPF_ZERO, 0.0);
    mpf_init_set_d(MPF_MINUS_ONE, -1.0);
    mpf_init_set_d(MPF_TINY, 1.0e-25);
}

static void mpf_constants_clear (void)
{
    mpf_clear(MPF_ONE);
    mpf_clear(MPF_ZERO);
    mpf_clear(MPF_MINUS_ONE);
    mpf_clear(MPF_TINY);
}

static void free_mpZ (mpf_t **mpZ, int v, int n)
{
    int i, t;

    for (i=0; i<v; i++) {
	if (mpZ[i] != NULL) {
	    for (t=0; t<n; t++) {
		mpf_clear(mpZ[i][t]);
	    }
	    free(mpZ[i]);
	}
    }
    free(mpZ);
}

/* reject the incoming data (a) if any values are missing, (b) if any
   vars are all zero */

static int data_problems (const int *list, const double **Z, 
			  const DATAINFO *pdinfo, char *errbuf)
{
    int i, t, allzero;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) {
	    continue; 
	}
	allzero = 1;
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) { 
	    if (!floateq(Z[list[i]][t], 0.0)) {
		allzero = 0;
	    }
	}
	if (allzero) {
	    sprintf(errbuf, _("Variable '%s' is all zeros"), 
		    pdinfo->varname[list[i]]);
	    return 1;
	}
    }

    return 0;
}

static void make_poly_series (MPMODEL *pmod, mpf_t **mpZ,
			      int pli, int ppos, int mpi)
{
    unsigned long pwr = pmod->polylist[pli];
    int t, s = 0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
#if MP_DEBUG
	printf("generating mpZ[%d][%d] from mpZ[%d][%d],\n"
	       "using power %lu taken from polylist[%d]\n", 
	       mpi, s, ppos, s, pwr, pli);
#endif
	mpf_init(mpZ[mpi][s]);
	mpf_pow_ui(mpZ[mpi][s],            /* target */
		   mpZ[ppos][s],           /* source */ 
		   pwr);                   /* power */
	s++;
    }
}

static void fill_mp_series (MPMODEL *pmod, const double **Z, mpf_t **mpZ,
			    const int *zdigits, int i, int mpi)
{
    char numstr[64];
    int t, s = 0; 

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (zdigits != NULL && zdigits[i] > 0 && zdigits[i] < DBL_DIG) {
	    /* do trick with strings */
	    sprintf(numstr, "%.*g", zdigits[i], Z[i][t]);
	    mpf_init_set_str(mpZ[mpi][s], numstr, 10);
	} else { 
	    /* do straight conversion */
	    mpf_init_set_d(mpZ[mpi][s], Z[i][t]);
	}
	s++;
    }
}

/* Given a data set containing doubles, build a data set using
   GMP's mpf_t floating point type.  

   For "ordinary" data we either (a) simply initialize the mpf_t
   straight from the corresponding double, or (b) attempt a "clever
   trick", namely print the double to a string using a precision that
   was recorded at the time the original data were read, then set the
   mpf_t from that string.  This is designed to avoid the transmission
   to the mpf_t of garbage lying beyond DBL_DIG into the double.
   The trick is applicable only for data read from some original
   source, e.g. the NIST data files; also, it works only if the
   original data had a precision of not more than DBL_DIG.  

   Besides converting ordinary data, this function is also used to
   generate powers of x in the case of a polynomial regression of y on
   x, as found in several of the NIST examples.  Accuracy of the
   regression results may suffer if the generation of successive
   powers of x is done in regular double precision, so we do it here.
*/

static mpf_t **make_mpZ (MPMODEL *mpmod, const int *zdigits,
			 const double **Z, const DATAINFO *pdinfo, 
			 char **xnames)
{
    int i, s, t;
    int n = mpmod->t2 - mpmod->t1 + 1;
    int l0 = mpmod->list[0];
    int npoly, mp_poly_pos = 0;
    int listpt, nvars = 0, v = 0;
    mpf_t **mpZ = NULL;
    int err = 0;

    if (n <= 0) {
	return NULL;
    }

    /* "varlist" holds the regression specification, using
       the numbering of variables from the original dataset 
    */

    mpmod->varlist = gretl_list_new(l0);
    if (mpmod->varlist == NULL) {
	return NULL;
    }

    mpZ = malloc(l0 * sizeof *mpZ);
    if (mpZ == NULL) {
	return NULL;
    }

    for (i=0; i<l0; i++) {
	mpZ[i] = NULL;
    }

    if (mpmod->ifc) { 
	mpZ[0] = malloc(n * sizeof **mpZ);
	s = 0;
	for (t=mpmod->t1; t<=mpmod->t2; t++) {
	    mpf_init_set_d(mpZ[0][s++], 1.0);
	}
	if (xnames != NULL) { 	 
	    strcpy(xnames[v++], pdinfo->varname[0]); 	 
	}
	nvars++;
    } 

    /* number of polynomial terms to be generated */
    if (mpmod->polylist != NULL) {
	npoly = mpmod->polylist[0];
    } else {
	npoly = 0;
    }

    /* process the ordinary data */
    for (i=1; i<=l0-npoly; i++) {
	if (mpmod->list[i] == 0) {
	    /* the constant is already handled */	    
	    mpmod->varlist[i] = 0;
	    continue;
	}
	mpZ[nvars] = malloc(n * sizeof **mpZ);
	if (mpZ[nvars] == NULL) {
	    err = 1;
	    break;
	}
	if (mpmod->list[i] == mpmod->polyvar) {
	    /* record position in mpZ of the var to be raised 
	       to various powers, if applicable */
#if MP_DEBUG
	    fprintf(stderr, "var to be raised to powers: it's "
		   "at position %d in the regression list,\n"
		   "and at slot %d in mpZ\n", i, nvars);
#endif
	    mp_poly_pos = nvars;
	}
	    
	fill_mp_series(mpmod, Z, mpZ, zdigits, mpmod->list[i], nvars); 
	mpmod->varlist[i] = mpmod->list[i];
        if (xnames != NULL && i > 1) { 	 
	    strcpy(xnames[v++], pdinfo->varname[mpmod->list[i]]); 	 
	}
	mpmod->list[i] = nvars;
	nvars++;
    } /* end processing ordinary data */

    listpt = i;

    /* generate polynomial data (if applicable) */

    for (i=0; i<npoly && !err; i++) {  
	mpZ[nvars] = malloc(n * sizeof **mpZ);
	if (mpZ[nvars] == NULL) {
	    err = 1;
	    break;
	}

	make_poly_series(mpmod, mpZ, i+1, mp_poly_pos, nvars);
	mpmod->varlist[i+listpt] = mpmod->polyvar;

        if (xnames != NULL) { 	 
	    sprintf(xnames[v++], "%s^%d", 	 
		    pdinfo->varname[mpmod->polyvar], 	 
		    mpmod->polylist[i+1]); 	 
	}

	mpmod->list[i+listpt] = nvars;
	nvars++;
    }

    if (err) {
	free_mpZ(mpZ, nvars, n);
	mpZ = NULL;
    }

    return mpZ;
}

/* budget version of "make_mpZ" when we're loading data from
   given y and X matrices */

static mpf_t **mpZ_from_matrices (const gretl_matrix *y, 
				  const gretl_matrix *X,
				  int *err)
{
    mpf_t **mpZ = NULL;
    int T = y->rows;
    int k = X->cols + 1;
    double x;
    int i, t;

    mpZ = malloc(k * sizeof *mpZ);
    if (mpZ == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<k; i++) {
	mpZ[i] = NULL;
    }

    for (i=0; i<k && !*err; i++) {
	mpZ[i] = malloc(T * sizeof **mpZ);
	if (mpZ[i] == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*err) {
	free_mpZ(mpZ, k, T);
	mpZ = NULL;
    } else {
	for (t=0; t<T; t++) {
	    x = gretl_vector_get(y, t);
	    mpf_init_set_d(mpZ[0][t], x);
	}
	for (i=1; i<k; i++) {
	    for (t=0; t<T; t++) {
		x = gretl_matrix_get(X, t, i-1);
		mpf_init_set_d(mpZ[i][t], x);
	    }
	}
    }

    return mpZ;
}

static void mp_model_free (MPMODEL *mpmod)
{
    int i, nx = 0, nt = 0;

    if (mpmod->list != NULL) {
	nx = mpmod->list[0] - 1;
	nt = nx * (nx + 1) / 2;
    }

    free(mpmod->list);
    free(mpmod->varlist);

    if (mpmod->coeff != NULL) {
	for (i=0; i<mpmod->ncoeff; i++) { 
	    mpf_clear(mpmod->coeff[i]);
	}
	free(mpmod->coeff);
    }

    if (mpmod->sderr != NULL) {    
	for (i=0; i<mpmod->ncoeff; i++) {
	    mpf_clear(mpmod->sderr[i]);
	}
	free(mpmod->sderr);
    }

    if (mpmod->xpx != NULL) {
	for (i=0; i<nt; i++) {
	    mpf_clear(mpmod->xpx[i]);
	}
	free(mpmod->xpx);
    }

    mpf_clear(mpmod->ess);
    mpf_clear(mpmod->tss);
    mpf_clear(mpmod->sigma);
    mpf_clear(mpmod->rsq);
    mpf_clear(mpmod->adjrsq);
    mpf_clear(mpmod->fstt); 
}

static void mp_model_init (MPMODEL *mpmod)
{
    mpmod->ID = 0;
    mpmod->t1 = 0;
    mpmod->t2 = 0;

    mpmod->list = NULL;
    mpmod->varlist = NULL;
    mpmod->polylist = NULL; /* don't free, the caller does that */
    mpmod->ifc = 1;
    mpmod->coeff = NULL;
    mpmod->sderr = NULL;
    mpmod->xpx = NULL;

    mpf_init(mpmod->ess);
    mpf_init(mpmod->tss);
    mpf_init(mpmod->sigma);
    mpf_init(mpmod->rsq);
    mpf_init(mpmod->adjrsq);
    mpf_init(mpmod->fstt);

    mpmod->errcode = 0;
    mpmod->polyvar = 0;
}

/**
 * mp_vector_raise_to_power:
 * @srcvec: source vector (doubles)
 * @targvec: vector to be filled in with results
 * @n: length of vector
 * @power: integer power to which elements of @srcvec should
 * be raised, using multiple precision arithmetic.
 *
 * Returns: 0 on success, error code on failure.
 */

int mp_vector_raise_to_power (const double *srcvec, double *targvec,
			      int n, unsigned pwr)
{
    int t;
    mpf_t src, targ;

    set_gretl_mp_bits();

    mpf_init(src);
    mpf_init(targ);

    for (t=0; t<n; t++) {
	if (na(srcvec[t])) {
	    targvec[t] = NADBL;
	    continue;
	}
	mpf_set_d(src, srcvec[t]);
	mpf_pow_ui(targ, src, (unsigned long) pwr);
	targvec[t] = mpf_get_d(targ);
    }

    mpf_clear(src);
    mpf_clear(targ);

    return 0;
}

#ifdef HAVE_MPFR

/**
 * mp_vector_ln:
 * @srcvec: source vector (doubles)
 * @targvec: vector to be filled in with results
 * @n: length of vector
 *
 * Returns: 0 on success, error code on failure.
 */

int mp_vector_ln (const double *srcvec, double *targvec, int n)
{
    int t;
    mpfr_t src, targ;

    set_gretl_mpfr_bits();

    mpfr_init(src);
    mpfr_init(targ);

    for (t=0; t<n; t++) {
	if (na(srcvec[t])) {
	    targvec[t] = NADBL;
	    continue;
	}
	mpfr_set_d(src, srcvec[t], GMP_RNDN);
	mpfr_log(targ, src, GMP_RNDN);
	targvec[t] = mpfr_get_d(targ, GMP_RNDN);
    }

    mpfr_clear(src);
    mpfr_clear(targ);

    return 0;
}

#endif

static int poly_check (MPMODEL *mpmod, const int *list)
{
    int i;

    /* check that all powers are > 1 */

    for (i=1; i<=mpmod->polylist[0]; i++) {
	if (mpmod->polylist[i] < 2) {
	    return 1;
	}
    }

    /* take the rightmost var in the regression list (other than 
       the constant) as the one to be raised to various powers */

    for (i=list[0]; i>1; i--) {
	if (list[i] != 0) {
	    mpmod->polyvar = list[i];
	    break;
	}
    }

    if (mpmod->polyvar == 0) {
	return 1;
    }

    return 0;
}

static int *poly_copy_list (const int *list, const int *poly)
{
    int *targ;
    int i;

    targ = gretl_list_new(list[0] + poly[0]);
    if (targ == NULL) {
	return NULL;
    }

    for (i=1; i<=list[0]; i++) {
	targ[i] = list[i];
    }

    for (i=1; i<=poly[0]; i++) {
	targ[list[0] + i] = list[0] + i - 1; 
    }  

    return targ;
}

static void set_gretl_mp_bits (void)
{
    char *user_bits = getenv("GRETL_MP_BITS");
    unsigned long bits = get_mp_bits();
    
    if (user_bits != NULL) {
	bits = strtoul(user_bits, NULL, 10);
    }

    fprintf(stderr, "GMP: using %d bits\n", (int) bits);

    mpf_set_default_prec(bits);
}

#ifdef HAVE_MPFR

static void set_gretl_mpfr_bits (void)
{
    char *user_bits = getenv("GRETL_MP_BITS");
    unsigned long bits = get_mp_bits();
    
    if (user_bits != NULL) {
	bits = strtoul(user_bits, NULL, 10);
    }

    mpfr_set_default_prec(bits);
} 

/* compute log-likelihood etc. in multiple precision, using the
   MPFR library
*/

static void mp_ll_stats (const MPMODEL *mpmod, MODEL *pmod)
{
    double n = mpmod->nobs;
    double k = mpmod->ncoeff;
    mpfr_t mll, mln, ll2;
    mpfr_t mx1, mx2;
    mpfr_t ln2pi1, crit;

    set_gretl_mpfr_bits();

    mpfr_init(mll);
    mpfr_init(mln);
    mpfr_init(ll2);
    mpfr_init(mx1);
    mpfr_init(mx2);
    mpfr_init(ln2pi1);
    mpfr_init(crit);

    mpfr_const_pi(ln2pi1, GMP_RNDN);          /* pi */
    mpfr_set_d(mx1, 2.0, GMP_RNDN);           /* 2 */
    mpfr_mul(ln2pi1, ln2pi1, mx1, GMP_RNDN);  /* 2 * pi */
    mpfr_log(ln2pi1, ln2pi1, GMP_RNDN);       /* log(2*pi) */
    mpfr_set_d(mx1, 1.0, GMP_RNDN);           /* 1 */
    mpfr_add(ln2pi1, ln2pi1, mx1, GMP_RNDN);  /* log(2*pi) + 1 */

    mpfr_set_f(mll, mpmod->ess, GMP_RNDN);
    mpfr_set_d(mx1, -.5, GMP_RNDN);
    mpfr_set_d(mx2, n, GMP_RNDN);

    mpfr_log(mll, mll, GMP_RNDN);          /* log(ess) */
    mpfr_log(mln, mx2, GMP_RNDN);          /* log(n) */

    mpfr_mul(mll, mx2, mll, GMP_RNDN);     /* n * log(ess) */
    mpfr_mul(mll, mx1, mll, GMP_RNDN);     /* -.5 * n * log(ess) */

    mpfr_mul(mx1, mx1, mx2, GMP_RNDN);     /* -.5 * n */
    mpfr_sub(mx2, ln2pi1, mln, GMP_RNDN);  /* log(2*pi) - log(n) */
    mpfr_mul(mx2, mx1, mx2, GMP_RNDN);     /* -.5 * n * (log(2*pi) - log(n)) */

    mpfr_add(mll, mll, mx2, GMP_RNDN);     /* now actual log-likelihood */  

    pmod->lnL = mpfr_get_d(mll, GMP_RNDN);

    if (xna(pmod->lnL)) {
	pmod->lnL = NADBL;
	mle_criteria(pmod, 0);
    } else {
	mpfr_set_d(mx1, -2.0, GMP_RNDN);
	mpfr_mul(ll2, mx1, mll, GMP_RNDN);    /* -2.0 * ll */

	mpfr_set_d(mx1, 2.0, GMP_RNDN);
	mpfr_set_d(mx2, k, GMP_RNDN);
	mpfr_mul(mx1, mx1, mx2, GMP_RNDN);    /* 2 * k */

	mpfr_add(crit, ll2, mx1, GMP_RNDN);   /* -2.0 * ll + 2 * k */
	pmod->criterion[C_AIC] = mpfr_get_d(crit, GMP_RNDN);

	mpfr_mul(mx1, mx2, mln, GMP_RNDN);    /* k * log(n) */
	mpfr_add(crit, ll2, mx1, GMP_RNDN);   /* -2.0 * ll + k * log(n) */
	pmod->criterion[C_BIC] = mpfr_get_d(crit, GMP_RNDN);

	mpfr_set_d(mx1, 2.0, GMP_RNDN);
	mpfr_set_d(mx2, (double) k, GMP_RNDN);
	mpfr_mul(mx1, mx1, mx2, GMP_RNDN);    /* 2 * k, again */
	mpfr_log(mx2, mln, GMP_RNDN);         /* log(log(n)) */
	mpfr_mul(mx1, mx1, mx2, GMP_RNDN);    /* 2 * k * log(log(n) */
	mpfr_add(crit, ll2, mx1, GMP_RNDN);   /* -2.0 * ll + 2 * k * log(log(n)) */
	pmod->criterion[C_HQC] = mpfr_get_d(crit, GMP_RNDN);
    }

    mpfr_clear(mll);
    mpfr_clear(mln);
    mpfr_clear(ll2);
    mpfr_clear(mx1);
    mpfr_clear(mx2);
    mpfr_clear(ln2pi1);
    mpfr_clear(crit);

    mpfr_free_cache();
}

#else

#include <errno.h>

/* compute log-likelihood etc., based on the ESS from the
   multiple-precision model but using ordinary double-precision
   arithmetic: a fallback if the MPFR library is not available
*/

static void mp_ll_stats (const MPMODEL *mpmod, MODEL *pmod)
{
    int k = mpmod->ncoeff;
    int n = mpmod->nobs;

    pmod->ess = mpf_get_d(mpmod->ess);

    if (pmod->ess < 0 || xna(pmod->ess)) {
	pmod->ess = pmod->lnL = NADBL;
    } else {
	const double ln2pi1 = 2.837877066409345;

	errno = 0;

	pmod->lnL = -.5 * n * log(pmod->ess);

	if (errno == EDOM || errno == ERANGE) {
	    pmod->lnL = NADBL;
	    errno = 0;
	} else {
	    pmod->lnL += -.5 * n * (ln2pi1 - log((double) n));
	}
    }

    pmod->ncoeff = k;
    pmod->nobs = n;

    mle_criteria(pmod, 0);
}

#endif

static void mp_dwstat (const MPMODEL *mpmod, MODEL *pmod,
		       mpf_t *uhat)
{
    mpf_t num, x;
    mpf_t ut1, u11;
    int t;

    mpf_init(num);
    mpf_init(x);
    mpf_init(ut1);
    mpf_init(u11);

    for (t=1; t<mpmod->nobs; t++)  {
	mpf_sub(x, uhat[t], uhat[t-1]);
	mpf_pow_ui(x, x, 2);
	mpf_add(num, num, x);
	mpf_mul(x, uhat[t], uhat[t-1]);
	mpf_add(ut1, ut1, x);
	mpf_mul(x, uhat[t-1], uhat[t-1]);
	mpf_add(u11, u11, x);
    }

    mpf_div(x, num, mpmod->ess);
    pmod->dw = mpf_get_d(x);

    if (isnan(pmod->dw) || isinf(pmod->dw)) {
	pmod->dw = NADBL;
    }    

    mpf_div(x, ut1, u11);
    pmod->rho = mpf_get_d(x);

    if (isnan(pmod->rho) || isinf(pmod->rho)) {
	pmod->dw = NADBL;
	pmod->rho = NADBL;
    }

    mpf_clear(num);
    mpf_clear(x);
    mpf_clear(ut1);
    mpf_clear(u11);
}

/* compute coefficient covariance matrix in multiple precision */

static int mp_makevcv (const MPMODEL *mpmod, MODEL *pmod,
		       gretl_matrix *V, double *ps2)
{
    mpf_t *vcv;
    int dec, mst, kk, i, j, kj, icnt, m, k, l = 0;
    const int nv = mpmod->ncoeff;
    const int nxpx = (nv * nv + nv) / 2; 
    mpf_t d, x, s2;

    if (pmod == NULL && V == NULL) {
	return 0;
    }

    if (mpmod->xpx == NULL) {
	return 1;
    }

    mpf_init(d);
    mpf_init(x);
    mpf_init(s2);

    mst = nxpx;
    kk = nxpx - 1;

    vcv = malloc(nxpx * sizeof *vcv);
    if (vcv == NULL) {
	return E_ALLOC;
    }

    if (pmod != NULL) {
	pmod->vcv = malloc(nxpx * sizeof *pmod->vcv);
	if (pmod->vcv == NULL) {
	    free(vcv);
	    return E_ALLOC;
	}
    }

    for (i=0; i<nxpx; i++) {
	mpf_init(vcv[i]);
    }

    for (i=0; i<nv; i++) {
	mst -= i;
	/* find diagonal element */
	mpf_set(d, mpmod->xpx[kk]);
	if (i > 0) {
	    for (j=kk+1; j<=kk+i; j++) {
		mpf_mul(x, mpmod->xpx[j], vcv[j]);
		mpf_sub(d, d, x);
	    }
	}
	mpf_mul(vcv[kk], d, mpmod->xpx[kk]);
	/* find off-diagonal elements indexed by kj */
	kj = kk;
	kk = kk - i - 2;
	if (i > nv - 2) {
	    continue;
	}
	for (j=i+1; j<nv; j++) {
	    icnt = i+1;
	    kj -= j;
	    mpf_set(d, MPF_ZERO);
	    m = mst + 1;
	    for (k=0; k<=j-1; k++) {
		if (icnt > 0) {
		    dec = 1;
		    icnt--;
		} else {
		    dec = k;
		}
		m -= dec;
		l = kj + i - k;
		mpf_mul(x, vcv[m-1], mpmod->xpx[l]);
		mpf_add(d, d, x);
	    }
	    mpf_mul(x, d, mpmod->xpx[l-1]);
	    mpf_neg(vcv[kj], x);
	}
    }

    if (pmod != NULL || ps2 != NULL) {
	/* get residual variance */
	mpf_mul(s2, mpmod->sigma, mpmod->sigma);
    }

    if (pmod != NULL) {
	for (i=0; i<nxpx; i++) {
	    mpf_mul(x, vcv[i], s2);
	    pmod->vcv[i] = mpf_get_d(x);
	    mpf_clear(vcv[i]);
	}
    } else {
	double dval;

	for (i=0; i<nv; i++) {
	    for (j=0; j<=i; j++) {
		k = ijton(i, j, nv);
		if (ps2 != NULL) {
		    mpf_mul(x, vcv[k], s2);
		    dval = mpf_get_d(x);
		} else {
		    dval = mpf_get_d(vcv[k]);
		}
		gretl_matrix_set(V, i, j, dval);
		gretl_matrix_set(V, j, i, dval);
		mpf_clear(vcv[k]);
	    }
	}
    }

    mpf_clear(d);
    mpf_clear(x);
    mpf_clear(s2);
    free(vcv);

    return 0;
}

/* compute mean and s.d. of dependent variable in multiple precision */

static void mp_depvarstats (const MPMODEL *mpmod, MODEL *pmod,
			    mpf_t **mpZ)
{
    mpf_t xbar, ssx, diff, mn;
    double xn = mpmod->nobs;
    int yno = mpmod->list[1];
    int t;
    
    mpf_init(xbar);
    mpf_init(ssx);
    mpf_init(diff);
    mpf_init(mn);

    mpf_set(xbar, MPF_ZERO);
    mpf_set(ssx, MPF_ZERO);
    mpf_set_d(mn, xn);

    for (t=0; t<mpmod->nobs; t++) {
	mpf_add(xbar, xbar, mpZ[yno][t]);
    }

    mpf_div(xbar, xbar, mn);
    pmod->ybar = mpf_get_d(xbar);

    for (t=0; t<mpmod->nobs; t++) {
	mpf_sub(diff, mpZ[yno][t], xbar);
	mpf_pow_ui(diff, diff, 2);
	mpf_add(ssx, ssx, diff);
    }

    mpf_sub(mn, mn, MPF_ONE);
    mpf_div(ssx, ssx, mn);
    mpf_sqrt(ssx, ssx);     
    pmod->sdy = mpf_get_d(ssx);

    mpf_clear(xbar);
    mpf_clear(ssx);
    mpf_clear(diff);
    mpf_clear(mn);
}

/* compute residuals and fitted values in multiple precision */

static void mp_hatvars (const MPMODEL *mpmod, MODEL *pmod,
			gretl_vector *uvec, mpf_t **mpZ, 
			int tseries)
{
    mpf_t *uhat = NULL;
    mpf_t yht, uht, xbi;
    int yno = mpmod->list[1];
    int i, t;

    if (tseries) {
	uhat = malloc(mpmod->nobs * sizeof *uhat);
	if (uhat != NULL) {
	    for (t=0; t<mpmod->nobs; t++) {
		mpf_init(uhat[t]);
	    }
	}
    }
    
    mpf_init(yht);
    mpf_init(uht);
    mpf_init(xbi);

    for (t=0; t<mpmod->nobs; t++) {
	mpf_set_d(yht, 0.0);
	for (i=0; i<mpmod->ncoeff; i++) {
	    mpf_mul(xbi, mpmod->coeff[i], mpZ[mpmod->list[i+2]][t]);
	    mpf_add(yht, yht, xbi);
	}
	mpf_sub(uht, mpZ[yno][t], yht);
	if (pmod != NULL) {
	    pmod->yhat[t + mpmod->t1] = mpf_get_d(yht);
	    pmod->uhat[t + mpmod->t1] = mpf_get_d(uht);
	} else if (uvec != NULL) {
	    uvec->val[t] = mpf_get_d(uht);
	}
	if (uhat != NULL) {
	    mpf_set(uhat[t], uht);
	}
    }

    mpf_clear(yht);
    mpf_clear(uht);
    mpf_clear(xbi);

    if (uhat != NULL) {
	mp_dwstat(mpmod, pmod, uhat);
	for (t=0; t<mpmod->nobs; t++) {
	    mpf_clear(uhat[t]);
	}
	free(uhat);
    } else if (pmod != NULL) {
	pmod->rho = pmod->dw = NADBL;
    }
}

static int copy_mp_results (const MPMODEL *mpmod, MODEL *pmod,
			    const DATAINFO *pdinfo, mpf_t **mpZ,
			    char **xnames, gretlopt opt)
{
    int tseries = dataset_is_time_series(pdinfo);
    int i, err = 0;

    pmod->ncoeff = mpmod->ncoeff;
    pmod->full_n = pdinfo->n;
    pmod->ci = MPOLS;

    err = gretl_model_allocate_storage(pmod);
    if (err) {
	if (xnames != NULL) { 	 
	    free_strings_array(xnames, pmod->ncoeff); 	 
	} 	 
	return err;
    }

    if (xnames != NULL) { 	 
	gretl_model_add_allocated_varnames(pmod, xnames); 	 
    } 	 

    for (i=0; i<mpmod->ncoeff; i++) {
	pmod->coeff[i] = mpf_get_d(mpmod->coeff[i]);
	pmod->sderr[i] = mpf_get_d(mpmod->sderr[i]);
    }

    pmod->sigma = mpf_get_d(mpmod->sigma);
    pmod->ess = mpf_get_d(mpmod->ess);
    pmod->rsq = mpf_get_d(mpmod->rsq);
    pmod->fstt = mpf_get_d(mpmod->fstt);
    pmod->chisq = NADBL;

    if (opt & OPT_S) {
	/* saving additional results */
	pmod->t1 = mpmod->t1;
	pmod->t2 = mpmod->t2;
	pmod->nobs = mpmod->nobs;
	pmod->ifc = mpmod->ifc;
	pmod->dfn = mpmod->dfn;
	pmod->dfd = mpmod->dfd;
	pmod->adjrsq = mpf_get_d(mpmod->adjrsq);
	pmod->list = gretl_list_copy(mpmod->varlist);
	if (pmod->list == NULL) {
	    err = E_ALLOC;
	} else {
	    mp_depvarstats(mpmod, pmod, mpZ);
	    mp_hatvars(mpmod, pmod, NULL, mpZ, tseries);
	    mp_ll_stats(mpmod, pmod);
	    mp_makevcv(mpmod, pmod, NULL, NULL);
	}
    }

    return err;
}

static int matrix_copy_mp_results (const MPMODEL *mpmod, mpf_t **mpZ,
				   gretl_vector *b, gretl_matrix *vcv, 
				   gretl_vector *uhat, double *s2)
{
    int i, err = 0;

    for (i=0; i<mpmod->ncoeff; i++) {
	b->val[i] = mpf_get_d(mpmod->coeff[i]);
    }

    if (vcv != NULL) {
	err = mp_makevcv(mpmod, NULL, vcv, s2);
    } else if (s2 != NULL) {
	mpf_t ms2;

	mpf_init(ms2);
	mpf_mul(ms2, mpmod->sigma, mpmod->sigma);
	*s2 = mpf_get_d(ms2);
	mpf_clear(ms2);
    }

    if (uhat != NULL) {
	mp_hatvars(mpmod, NULL, uhat, mpZ, 0);
    }

    return err;
}

static char **allocate_xnames (const int *list) 	 
{ 	 
    int n = list[0] - 1; 	 
    char **s = strings_array_new_with_length(n, VNAMELEN + 6);

    return s;
}

static void mp_xpxxpy_init (MPXPXXPY *m)
{
    m->xpy = NULL;
    m->xpx = NULL;
    m->errcode = 0;
    m->nv = 0;
    m->ivalue = 0;
}

static MPXPXXPY mp_xpxxpy_func (const int *list, int n, mpf_t **mpZ)
{
    int i, j, li, lj, m, t;
    const int l0 = list[0];
    const int yno = list[1];
    mpf_t xx, yy, z1, z2, tmp;
    MPXPXXPY xpxxpy;

    mp_xpxxpy_init(&xpxxpy);

    i = l0 - 1;
    m = i * (i + 1) / 2;

    if ((xpxxpy.xpy = malloc((l0 + 1) * sizeof *xpxxpy.xpy)) == NULL ||
	(xpxxpy.xpx = malloc(m * sizeof *xpxxpy.xpx)) == NULL) {
        xpxxpy.errcode = E_ALLOC;
        return xpxxpy;
    }

    for (i=0; i<=l0; i++) {
	mpf_init(xpxxpy.xpy[i]);
    }

    for (i=0; i<m; i++) {
	mpf_init(xpxxpy.xpx[i]);
    }

    mpf_init(xx);
    mpf_init(yy);
    mpf_init(z1);
    mpf_init(z2);
    mpf_init(tmp);

    xpxxpy.nv = l0 - 1;

    for (t=0; t<n; t++) {
        mpf_set(xx, mpZ[yno][t]);
	mpf_add(xpxxpy.xpy[0], xpxxpy.xpy[0], xx);
	mpf_mul(yy, xx, xx);
	mpf_add(xpxxpy.xpy[l0], xpxxpy.xpy[l0], yy);
    }

    if (mpf_sgn(xpxxpy.xpy[l0]) == 0) {
         xpxxpy.ivalue = yno; 
         return xpxxpy; 
    }    

    m = 0;

    for (i=2; i<=l0; i++) {
        li = list[i];
        for (j=i; j<=l0; j++) {
            lj = list[j];
            mpf_set(xx, MPF_ZERO);
            for (t=0; t<n; t++) {
		mpf_mul(tmp, mpZ[li][t], mpZ[lj][t]);
		mpf_add(xx, xx, tmp);
	    }
            if ((mpf_sgn(xx) == 0) && li == lj)  {
                xpxxpy.ivalue = li;
                return xpxxpy;  
            }
            mpf_set(xpxxpy.xpx[m++], xx);
        }
        mpf_set(xx, MPF_ZERO);
        for (t=0; t<n; t++) {
	    mpf_mul(tmp, mpZ[yno][t], mpZ[li][t]);
	    mpf_add(xx, xx, tmp);
	}
        mpf_set(xpxxpy.xpy[i-1], xx);
    }

    xpxxpy.ivalue = 0;

    mpf_clear(xx);
    mpf_clear(yy);
    mpf_clear(z1);
    mpf_clear(z2);
    mpf_clear(tmp);

    return xpxxpy; 
}

static MPCHOLBETA mp_cholbeta (MPXPXXPY xpxxpy)
{
    int i, j, k, kk, l, jm1, nv;
    mpf_t e, d, d1, test, rtest, xx, tmp;
    MPCHOLBETA cb;

    nv = xpxxpy.nv; 
    cb.errcode = 0; 
    mpf_init (cb.rss);

    if ((cb.coeff = malloc(nv * sizeof *cb.coeff)) == NULL) {
        cb.errcode = E_ALLOC;
        return cb;
    }

    for (j=0; j<nv; j++) {
	mpf_init (cb.coeff[j]);
    }

    mpf_init(e);
    mpf_init(d);
    mpf_init(d1);
    mpf_init(test);
    mpf_init(rtest);
    mpf_init(xx);
    mpf_init(tmp);

    cb.xpxxpy = xpxxpy;

    mpf_sqrt(tmp, xpxxpy.xpx[0]);
    mpf_div(e, MPF_ONE, tmp);
    mpf_set(xpxxpy.xpx[0], e);
    mpf_mul(xpxxpy.xpy[1], xpxxpy.xpy[1], e);
    for (i=1; i<nv; i++) { 
	mpf_mul(xpxxpy.xpx[i], xpxxpy.xpx[i], e);
    }

    kk = nv;

    for (j=2; j<=nv; j++) {
	/* diagonal elements */
	mpf_set(d, MPF_ZERO);
	mpf_set(d1, MPF_ZERO);
        k = jm1 = j - 1;
        for (l=1; l<=jm1; l++) {
	    mpf_set(xx, xpxxpy.xpx[k]);
	    mpf_mul(tmp, xx, xpxxpy.xpy[l]);
	    mpf_add(d1, d1, tmp);
	    mpf_mul(tmp, xx, xx);
	    mpf_add(d, d, tmp);
            k += nv-l;
        }
	mpf_sub(test, xpxxpy.xpx[kk], d);
	mpf_div(rtest, test, xpxxpy.xpx[kk]);
        if (mpf_sgn(test) != 1 || mpf_cmp(rtest, MPF_TINY) < 0) {
	    fprintf(stderr, "mp_cholbeta: rtest = %g\n", mpf_get_d(rtest));
	    mpf_set(cb.rss, MPF_MINUS_ONE); 
	    goto mp_cholbeta_abort;
        }   
	mpf_sqrt(tmp, test);
	mpf_div(e, MPF_ONE, tmp);
	mpf_set(xpxxpy.xpx[kk], e);
	mpf_sub(tmp, xpxxpy.xpy[j], d1);
	mpf_mul(xpxxpy.xpy[j], tmp, e);

        /* off-diagonal elements */
        for (i=j+1; i<=nv; i++) {
            kk++;
            mpf_set(d, MPF_ZERO);
            k = j - 1;
            for (l=1; l<=jm1; l++) {
		mpf_mul(tmp, xpxxpy.xpx[k], xpxxpy.xpx[k-j+i]);
		mpf_add(d, d, tmp);
                k += nv - l;
            }
	    mpf_sub(tmp, xpxxpy.xpx[kk], d);
	    mpf_mul(xpxxpy.xpx[kk], tmp, e);
        }
        kk++;
    }

    kk--;
    mpf_set(d, MPF_ZERO);

    for (j=1; j<=nv; j++) {
	mpf_mul(tmp, xpxxpy.xpy[j], xpxxpy.xpy[j]);
	mpf_add(d, d, tmp);
    }

    mpf_set(cb.rss, d);
    mpf_mul(cb.coeff[nv-1], xpxxpy.xpy[nv], xpxxpy.xpx[kk]);

    for (j=nv-1; j>=1; j--) {
	mpf_set(d, xpxxpy.xpy[j]);
        for (i=nv-1; i>=j; i--) {
            kk--;
	    mpf_mul(tmp, cb.coeff[i], xpxxpy.xpx[kk]);
	    mpf_sub(d, d, tmp);
        }
        kk--;
	mpf_mul(cb.coeff[j-1], d, xpxxpy.xpx[kk]);
    }  

 mp_cholbeta_abort:

    mpf_clear(e);
    mpf_clear(d);
    mpf_clear(d1);
    mpf_clear(test);
    mpf_clear(rtest);
    mpf_clear(xx);
    mpf_clear(tmp);
 
    return cb; 
}

static void mp_diaginv (MPXPXXPY xpxxpy, mpf_t *diag)
{
    int kk, l, m, k, i, j;
    const int nv = xpxxpy.nv;
    const int nxpx = nv * (nv + 1) / 2;
    mpf_t d, e, tmp;

    mpf_init(d);
    mpf_init(e);
    mpf_init(tmp);

    kk = 0;

    for (l=1; l<=nv-1; l++) {
	mpf_set(d, xpxxpy.xpx[kk]);
	mpf_set(xpxxpy.xpy[l], d);
	mpf_mul(e, d, d);
        m = 0;
        if (l > 1) {
	    for (j=1; j<=l-1; j++) {
		m += nv - j;
	    }
	}
        for (i=l+1; i<=nv; i++) {
	    mpf_set(d, MPF_ZERO); 
            k = i + m - 1;
            for (j=l; j<=i-1; j++) {
		mpf_mul(tmp, xpxxpy.xpy[j], xpxxpy.xpx[k]);
		mpf_add(d, d, tmp);
                k += nv - j;
            }
	    mpf_mul(d, d, xpxxpy.xpx[k]);
	    mpf_mul(d, d, MPF_MINUS_ONE);
	    mpf_set(xpxxpy.xpy[i], d);
	    mpf_mul(tmp, d, d);
	    mpf_add(e, e, tmp);
        }
        kk += nv + 1 - l;
        mpf_set(diag[l-1], e);
    }

    mpf_mul(diag[nv-1], xpxxpy.xpx[nxpx-1], xpxxpy.xpx[nxpx-1]);

    mpf_clear(d);
    mpf_clear(e);
    mpf_clear(tmp);
}

static void mp_regress (MPMODEL *pmod, MPXPXXPY xpxxpy, char *errbuf)
{
    int i, v, nobs, nv, yno;
    mpf_t *diag, ysum, ypy, zz, rss, tss;
    mpf_t den, sgmasq, tmp;
    double ess;
    MPCHOLBETA cb;

    nv = xpxxpy.nv;
    yno = pmod->list[1];

    pmod->sderr = malloc(nv * sizeof *pmod->sderr);
    if (pmod->sderr == NULL) {
        pmod->errcode = E_ALLOC;
        return;
    }

    for (i=0; i<nv; i++) {
	mpf_init (pmod->sderr[i]);
    }

    mpf_init(den);
    mpf_init(sgmasq);
    mpf_init(ysum);
    mpf_init(ypy);
    mpf_init(zz);
    mpf_init(rss);
    mpf_init(tss);
    mpf_init(tmp);

    nobs = pmod->nobs;
    pmod->ncoeff = nv;
    pmod->dfd = nobs - nv;

    if (pmod->dfd < 0) { 
       pmod->errcode = E_DF; 
       return; 
    }

    pmod->dfn = nv - pmod->ifc;
    mpf_set(ysum, xpxxpy.xpy[0]);
    mpf_set(ypy, xpxxpy.xpy[nv + 1]);
    if (mpf_sgn(ypy) == 0) { 
        pmod->errcode = E_ZERO;
        return; 
    }

    mpf_mul(zz, ysum, ysum);
    mpf_set_d(tmp, (double) nobs);
    mpf_div(zz, zz, tmp);
    mpf_sub(tss, ypy, zz);
    if (mpf_sgn(tss) < 0) { 
        pmod->errcode = E_TSS; 
        return; 
    }

    /* Choleski-decompose X'X and find the coefficients */
    cb = mp_cholbeta(xpxxpy);
    pmod->coeff = cb.coeff;
    pmod->xpx = cb.xpxxpy.xpx;

    if (cb.errcode) {
        pmod->errcode = E_ALLOC;
        return;
    } 

    mpf_set(rss, cb.rss);
    mpf_clear(cb.rss);
    if (mpf_cmp(rss, MPF_MINUS_ONE) == 0) { 
        pmod->errcode = E_SINGULAR;
        return; 
    }

    mpf_sub(pmod->ess, ypy, rss);
    ess = mpf_get_d(pmod->ess);
    if (fabs(ess) < DBL_EPSILON) {
	mpf_set(pmod->ess, MPF_ZERO);
    }

    if (mpf_sgn(pmod->ess) < 0) { 
	if (errbuf != NULL) {
	    sprintf(errbuf, _("Error sum of squares is not >= 0"));
	}
	pmod->errcode = E_DATA;
        return; 
    }

    if (pmod->dfd == 0) {
	mpf_set(pmod->sigma, MPF_ZERO);
	mpf_set_d(pmod->adjrsq, NADBL);
    } else {
	mpf_set_d(tmp, (double) pmod->dfd);
	mpf_div(sgmasq, pmod->ess, tmp);
	mpf_sqrt(pmod->sigma, sgmasq);
	mpf_mul(den, tss, tmp);
    }

    if (mpf_sgn(tss) <= 0) {
	mpf_set_d(pmod->rsq, NADBL);
	mpf_set_d(pmod->adjrsq, NADBL);
	pmod->errcode = E_TSS;
	return;
    }       

    if (pmod->errcode) {
	fprintf(stderr, "mp_ols: pmod->errcode = %d\n", pmod->errcode);
	return;
    }

    mpf_div(tmp, pmod->ess, tss);
    mpf_sub(pmod->rsq, MPF_ONE, tmp);

    if (pmod->dfd > 0) {
	mpf_set_d(tmp, (double) (nobs - 1));
	mpf_div(tmp, tmp, den);
	mpf_mul(tmp, tmp, pmod->ess);
	mpf_sub(pmod->adjrsq, MPF_ONE, tmp);
	if (!pmod->ifc) { 
	    mpf_t df;

	    mpf_div(tmp, pmod->ess, ypy);
	    mpf_sub(pmod->rsq, MPF_ONE, tmp);
	    mpf_sub(tmp, MPF_ONE, pmod->rsq);
	    mpf_init_set_d(df, (double) (nobs - 1));
	    mpf_mul(tmp, tmp, df);
	    mpf_set_d(df, (double) pmod->dfd);
	    mpf_div(tmp, tmp, df);
	    mpf_sub(pmod->adjrsq, MPF_ONE, tmp);
	    mpf_clear(df);
	}
    }

    if (pmod->ifc && nv == 1) {
        mpf_set(zz, MPF_ZERO);
        pmod->dfn = 1;
    }

    if (mpf_sgn(sgmasq) != 1 || pmod->dfd == 0) {
	mpf_set_d(pmod->fstt, NADBL);
    } else { 
	mpf_set_d(tmp, (double) pmod->ifc);
	mpf_mul(tmp, zz, tmp);
	mpf_sub(pmod->fstt, rss, tmp);
	mpf_div(pmod->fstt, pmod->fstt, sgmasq);
	mpf_set_d(tmp, (double) pmod->dfn);
	mpf_div(pmod->fstt, pmod->fstt, tmp);
    }

    diag = malloc(nv * sizeof *diag); 
    if (diag == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }

    for (i=0; i<nv; i++) {
	mpf_init(diag[i]);
    }

    mp_diaginv(xpxxpy, diag);

    for (v=0; v<nv; v++) { 
	mpf_sqrt(zz, diag[v]);
	mpf_mul(pmod->sderr[v], pmod->sigma, zz);
    }

    for (i=0; i<nv; i++) {
	mpf_clear(diag[i]);
    }

    free(diag); 

    mpf_clear(den);
    mpf_clear(sgmasq);
    mpf_clear(ysum);
    mpf_clear(ypy);
    mpf_clear(zz);
    mpf_clear(rss);
    mpf_clear(tss);
    mpf_clear(tmp);
    
    return;  
}

/* checks a list for a constant term (ID # 0), and if present, 
   move it to the first indep var position.  Return 1 if constant found,
   else 0.
*/

static int mp_rearrange (int *list)
{
    int i, v;

    for (v=list[0]; v>2; v--) {
        if (list[v] == 0)  {
	    for (i=v; i>2; i--) {
		list[i] = list[i-1];
	    }
	    list[2] = 0;
	    return 1;
        }
    }

    return (list[2] == 0);
}

/* public interface */

/**
 * mplsq:
 * @list: dependent variable plus list of regressors.
 * @polylist: list of polynomial terms (or %NULL).
 * @zdigits: list of digits of input precision for data series
 * (or %NULL).
 * @Z: data array.
 * @pdinfo: information on the data set.
 * @errbuf: where to print any error message.
 * @pmod: MODEL pointer to hold results.
 * @opt: if contains %OPT_S, save additional model
 * information (including the names of parameters in 
 * @pmod, if required).
 *
 * Computes multiple-precision OLS estimates of the model 
 * specified by @list, and stores them in @pmod.
 * 
 * Returns: 0 on success, error code on failure.
 */

int mplsq (const int *list, const int *polylist, const int *zdigits,
	   const double **Z, DATAINFO *pdinfo, 
	   char *errbuf, MODEL *pmod, gretlopt opt) 
{
    int l0, i;
    mpf_t **mpZ = NULL;
    char **xnames = NULL;
    MPXPXXPY xpxxpy;
    MPMODEL mpmod;
    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;
    int err = 0;

    *errbuf = 0;

    if (list == NULL || Z == NULL || pdinfo == NULL ||
	list[0] < 2 || pdinfo->v < 2) {
	return E_DATA;
    }

    set_gretl_mp_bits();

    mp_model_init(&mpmod);

    mpmod.t1 = pdinfo->t1;
    mpmod.t2 = pdinfo->t2;

    if (polylist == NULL) {
	mpmod.list = gretl_list_copy(list);
    } else {
	mpmod.list = poly_copy_list(list, polylist);
    }

    if (mpmod.list == NULL) {
	return E_ALLOC;
    }

    mpmod.polylist = polylist; /* attached for convenience */

    if (polylist != NULL && poly_check(&mpmod, list)) {
	err = E_DATA;
	goto bailout;
    }

    /* check for missing obs in sample */
    err = list_adjust_t1t2(list, Z, pdinfo);
    if (err) {
	goto bailout;
    }
    
    /* in case of any changes */
    mpmod.t1 = pdinfo->t1;
    mpmod.t2 = pdinfo->t2;

    /* check for other data issues */
    if (data_problems(list, Z, pdinfo, errbuf)) {
	err = E_DATA;
	goto bailout;
    }

    /* enable names for polynomial terms? */
    if (polylist != NULL && (opt & OPT_S)) { 	 
	xnames = allocate_xnames(mpmod.list); 	 
	if (xnames == NULL) { 
	    err = E_ALLOC;
	    goto bailout;
	} 	 
    }

    /* see if the regressor list contains a constant */
    mpmod.ifc = mp_rearrange(mpmod.list);

    /* construct multiple-precision data matrix */
    mpZ = make_mpZ(&mpmod, zdigits, Z, pdinfo, xnames);

    if (mpZ == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    mpf_constants_init();

    l0 = mpmod.list[0];  
    mpmod.ncoeff = l0 - 1; 
    mpmod.nobs = mpmod.t2 - mpmod.t1 + 1;

    /* check degrees of freedom */
    if (mpmod.nobs < mpmod.ncoeff) { 
        sprintf(errbuf, _("No. of obs (%d) is less than no. "
			  "of parameters (%d)"), mpmod.nobs, mpmod.ncoeff);
	free_mpZ(mpZ, l0, mpmod.nobs);
	mpf_constants_clear();
	err = E_DF;
	goto bailout;
    }

    /* calculate regression results */

    xpxxpy = mp_xpxxpy_func(mpmod.list, mpmod.nobs, mpZ);
    mpf_set(mpmod.tss, xpxxpy.xpy[l0]);

    mp_regress(&mpmod, xpxxpy, errbuf);

    for (i=0; i<=l0; i++) {
	mpf_clear(xpxxpy.xpy[i]);
    }
    free(xpxxpy.xpy);

    err = mpmod.errcode;
    if (!err) {
	err = copy_mp_results(&mpmod, pmod, pdinfo, mpZ, xnames, opt);
    } 

    /* free all the mpf stuff */
    free_mpZ(mpZ, l0, mpmod.nobs);
    mpf_constants_clear();

 bailout:

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    mp_model_free(&mpmod);

    return err;
}

int matrix_mp_ols (const gretl_vector *y, const gretl_matrix *X,
		   gretl_vector *b, gretl_matrix *vcv, 
		   gretl_vector *uhat, double *s2)
{
    mpf_t **mpZ;
    MPXPXXPY xpxxpy;
    MPMODEL mpmod;
    int *list;
    int i, k, T, l0;
    int err = 0;

    k = X->cols;
    T = y->rows;

    if (X->rows != T) {
	return E_NONCONF;
    }

    if (T < k) { 
	return E_DF;
    }

    list = gretl_consecutive_list_new(0, k);
    if (list == NULL) {
	return E_ALLOC;
    }

    set_gretl_mp_bits();
    mp_model_init(&mpmod);
    mpmod.t2 = T - 1;
    mpmod.list = list;

    /* construct multiple-precision data matrix */
    mpZ = mpZ_from_matrices(y, X, &err);
    if (err) {
	goto bailout;
    }

    mpf_constants_init();
    l0 = mpmod.list[0];
    mpmod.ncoeff = k; 
    mpmod.nobs = T;

    /* calculate regression results */

    xpxxpy = mp_xpxxpy_func(mpmod.list, T, mpZ);
    mpf_set(mpmod.tss, xpxxpy.xpy[l0]);

    mp_regress(&mpmod, xpxxpy, NULL);

    for (i=0; i<=l0; i++) {
	mpf_clear(xpxxpy.xpy[i]);
    }
    free(xpxxpy.xpy);

    err = mpmod.errcode;
    if (!err) {
	err = matrix_copy_mp_results(&mpmod, mpZ, b, vcv, uhat, s2);
    } 

    /* free all the mpf stuff */
    free_mpZ(mpZ, l0, mpmod.nobs);
    mpf_constants_clear();

 bailout:

    mp_model_free(&mpmod);

    return err;
}
