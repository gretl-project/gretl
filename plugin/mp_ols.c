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
#include "version.h"
#include "libset.h"

#include <float.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>

#define MP_DEBUG 0

static mpf_t MPF_ONE;
static mpf_t MPF_ZERO;
static mpf_t MPF_MINUS_ONE;
static mpf_t MPF_TINY;

typedef struct {
    int ID;                      /* ID number for model */
    int t1, t2, nobs;            /* starting observation, ending
                                    observation, and number of obs */
    char *mask;                  /* missing obs mask */
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

static void real_set_mp_bits (void)
{
    unsigned long bits = get_mp_bits();

#if 0
    fprintf(stderr, "GMP: using %d bits\n", (int) bits);
#endif
    mpf_set_default_prec(bits);
}

static void real_set_mpfr_bits (void)
{
    unsigned long bits = get_mp_bits();

#if 0
    fprintf(stderr, "MPFR: using %d bits\n", (int) bits);
#endif
    mpfr_set_default_prec(bits);
}

static void mpf_constants_init (void)
{
    mpf_init_set_d(MPF_ONE, 1.0);
    mpf_init_set_d(MPF_ZERO, 0.0);
    mpf_init_set_d(MPF_MINUS_ONE, -1.0);
    mpf_init_set_d(MPF_TINY, 1.0e-30);
}

static void mpf_constants_clear (void)
{
    mpf_clear(MPF_ONE);
    mpf_clear(MPF_ZERO);
    mpf_clear(MPF_MINUS_ONE);
    mpf_clear(MPF_TINY);
}

static void mpf_2d_array_free (mpf_t **X, int v, int n)
{
    int i, t;

    if (X == NULL) {
	return;
    }

    for (i=0; i<v; i++) {
	if (X[i] != NULL) {
	    for (t=0; t<n; t++) {
		mpf_clear(X[i][t]);
	    }
	    free(X[i]);
	}
    }
    free(X);
}

static void mpfr_2d_array_free (mpfr_t **X, int v, int n)
{
    int i, t;

    if (X == NULL) {
	return;
    }

    for (i=0; i<v; i++) {
	if (X[i] != NULL) {
	    for (t=0; t<n; t++) {
		mpfr_clear(X[i][t]);
	    }
	    free(X[i]);
	}
    }
    free(X);
}

static mpf_t **mpf_2d_array_alloc (int v, int n)
{
    mpf_t **X = malloc(v * sizeof *X);
    int i, j;

    if (X != NULL) {
	for (i=0; i<v; i++) {
	    X[i] = malloc(n * sizeof **X);
	    if (X[i] == NULL) {
		for (j=0; j<i; j++) {
		    free(X[j]);
		}
		free(X);
		return NULL;
	    }
	}
    }

    return X;
}

/* somewhat arbitrary */

#define eqzero(x) (fabs(x) < 1.0e-300)

/* reject the incoming data if any vars are all-zero */

static int data_problems (const int *list, const DATASET *dset)
{
    int i, t, allzero;
    double xit;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) {
	    continue;
	}
	allzero = 1;
	for (t=dset->t1; t<=dset->t2; t++) {
	    xit = dset->Z[list[i]][t];
	    if (!na(xit) && !eqzero(xit)) {
		allzero = 0;
		break;
	    }
	}
	if (allzero) {
	    gretl_errmsg_sprintf(_("Variable '%s' is all zeros"),
				 dset->varname[list[i]]);
	    return E_DATA;
	}
    }

    return 0;
}

#define mpmod_missing(m,t) (m->mask != NULL && m->mask[t-m->t1] != 0)

static void make_poly_series (MPMODEL *mpmod, mpf_t **mpZ,
			      int i, int ppos, int mpi)
{
    unsigned long pwr = mpmod->polylist[i];
    int t, s = 0;

    for (t=mpmod->t1; t<=mpmod->t2; t++) {
	if (mpmod_missing(mpmod, t)) {
	    continue;
	}
#if MP_DEBUG
	printf("generating mpZ[%d][%d] from mpZ[%d][%d],\n"
	       "using power %lu taken from polylist[%d]\n",
	       mpi, s, ppos, s, pwr, i);
#endif
	mpf_init(mpZ[mpi][s]);
	mpf_pow_ui(mpZ[mpi][s],   /* target */
		   mpZ[ppos][s],  /* source */
		   pwr);          /* power */
	s++;
    }
}

static void fill_mp_series (MPMODEL *mpmod, const DATASET *dset,
			    mpf_t **mpZ, const int *zdigits,
			    int i, int mpi)
{
    char numstr[64];
    int t, s = 0;

    for (t=mpmod->t1; t<=mpmod->t2; t++) {
	if (mpmod_missing(mpmod, t)) {
	    continue;
	}
	if (zdigits != NULL && zdigits[i] > 0 && zdigits[i] < DBL_DIG) {
	    /* do trick with strings */
	    sprintf(numstr, "%.*g", zdigits[i], dset->Z[i][t]);
	    mpf_init_set_str(mpZ[mpi][s], numstr, 10);
	} else {
	    /* do straight conversion */
	    mpf_init_set_d(mpZ[mpi][s], dset->Z[i][t]);
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
			 const DATASET *dset, char **xnames)
{
    int i, j, k, t, xn;
    int l0 = mpmod->list[0];
    int n = mpmod->nobs;
    int nv, npoly, poly_pos = 0;
    mpf_t **mpZ = NULL;

    /* "varlist" holds the regression specification, using
       the numbering of variables from the original dataset
    */

    mpmod->varlist = gretl_list_copy(mpmod->list);
    if (mpmod->varlist == NULL) {
	return NULL;
    }

    mpZ = mpf_2d_array_alloc(l0, n);
    if (mpZ == NULL) {
	return NULL;
    }

#if 0
    fprintf(stderr, "mpZ: allocated %d vars, %d obs\n", l0, n);
#endif

    /* number of polynomial terms to be generated */
    npoly = mpmod->polylist != NULL ? mpmod->polylist[0] : 0;

    /* number of ordinary vars in list */
    nv = l0 - npoly;

    /* start Z-column and x-name counters */
    j = xn = 0;

    if (mpmod->ifc || mpmod->list[1] == 0) {
	/* handle the constant first */
	for (t=0; t<n; t++) {
	    mpf_init_set_d(mpZ[j][t], 1.0);
	}
	if (xnames != NULL) {
	    strcpy(xnames[xn++], dset->varname[0]);
	}
	j++;
    }

    /* start polylist read-position counter */
    k = 1;

    for (i=1; i<=l0; i++) {
	int vi = mpmod->list[i];

	if (vi == 0) {
	    /* the constant */
	    mpmod->list[i] = 0;
	} else if (i <= nv) {
	    /* a regular variable */
	    if (vi == mpmod->polyvar) {
		/* record column in mpZ of the var to be raised
		   to various powers, if applicable */
		poly_pos = j;
	    }
	    fill_mp_series(mpmod, dset, mpZ, zdigits, vi, j);
	    if (xnames != NULL && i > 1) {
		strcpy(xnames[xn++], dset->varname[vi]);
	    }
	    mpmod->list[i] = j++;
	} else {
	    /* a polynomial term */
	    make_poly_series(mpmod, mpZ, k, poly_pos, j);
	    mpmod->varlist[i] = mpmod->polyvar;
	    if (xnames != NULL) {
		sprintf(xnames[xn++], "%s^%d",
			dset->varname[mpmod->polyvar],
			mpmod->polylist[k]);
	    }
	    mpmod->list[i] = j++;
	    k++;
	}
    }

    return mpZ;
}

/* budget version of "make_mpZ" when we're loading data from
   given y and X matrices */

static mpf_t **mpZ_from_matrices (const gretl_matrix *y,
				  const gretl_matrix *X,
				  int *err)
{
    mpf_t **mpZ;
    int T = y->rows;
    int k = X->cols + 1;
    double x;
    int i, t;

    mpZ = mpf_2d_array_alloc(k, T);

    if (mpZ == NULL) {
	*err = E_ALLOC;
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
    free(mpmod->mask);

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
    mpmod->mask = NULL;

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

    real_set_mp_bits();

    mpf_init(src);
    mpf_init(targ);

    for (t=0; t<n; t++) {
	if (na(srcvec[t])) {
	    targvec[t] = NADBL;
	} else {
	    mpf_set_d(src, srcvec[t]);
	    mpf_pow_ui(targ, src, (unsigned long) pwr);
	    targvec[t] = mpf_get_d(targ);
	}
    }

    mpf_clear(src);
    mpf_clear(targ);

    return 0;
}

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

    real_set_mpfr_bits();

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

static int poly_check (MPMODEL *mpmod, const int *polylist, const int *list)
{
    int i;

    /* check that all powers of x are > 1 */

    for (i=1; i<=polylist[0]; i++) {
	if (polylist[i] < 2) {
	    return E_INVARG;
	}
    }

    /* take the rightmost var in the regression list (other than
       the constant) as the one to be raised to various powers
    */

    for (i=list[0]; i>1; i--) {
	if (list[i] != 0) {
	    mpmod->polyvar = list[i];
	    break;
	}
    }

    return mpmod->polyvar == 0 ? E_DATA : 0;
}

static int *poly_copy_list (const int *list, const int *poly)
{
    int *targ;
    int i, n = list[0] + poly[0];

    targ = gretl_list_new(n);

    if (targ == NULL) {
	return NULL;
    }

    for (i=1; i<=n; i++) {
	if (i <= list[0]) {
	    targ[i] = list[i];
	} else {
	    targ[i] = i - 1;
	}
    }

    return targ;
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

    real_set_mpfr_bits();

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

    mpfr_set(mll, mpmod->ess, GMP_RNDN);
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

    if (na(pmod->lnL)) {
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

static void mp_dwstat (const MPMODEL *mpmod, MODEL *pmod,
		       mpf_t *uhat)
{
    mpf_t num, x;
    mpf_t ut1, u11;
    int t;

    if (mpmod->mask != NULL || mpf_sgn(mpmod->ess) <= 0) {
	pmod->dw = pmod->rho = NADBL;
	return;
    }

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

    if (mpf_sgn(u11) <= 0) {
	pmod->rho = NADBL;
    } else {
	mpf_div(x, ut1, u11);
	pmod->rho = mpf_get_d(x);
	if (isnan(pmod->rho) || isinf(pmod->rho)) {
	    pmod->dw = NADBL;
	    pmod->rho = NADBL;
	}
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
    int i, vi, t, s;

    if (tseries && mpmod->mask == NULL) {
	uhat = malloc(mpmod->nobs * sizeof *uhat);
	if (uhat != NULL) {
	    for (s=0; s<mpmod->nobs; s++) {
		mpf_init(uhat[s]);
	    }
	}
    }

    mpf_init(yht);
    mpf_init(uht);
    mpf_init(xbi);

    s = 0;

    for (t=mpmod->t1; t<=mpmod->t2; t++) {
	if (mpmod_missing(mpmod, t)) {
	    continue;
	}
	mpf_set_d(yht, 0.0);
	for (i=0; i<mpmod->ncoeff; i++) {
	    vi = mpmod->list[i+2];
	    mpf_mul(xbi, mpmod->coeff[i], mpZ[vi][s]);
	    mpf_add(yht, yht, xbi);
	}
	mpf_sub(uht, mpZ[yno][s], yht);
	if (pmod != NULL) {
	    pmod->yhat[t] = mpf_get_d(yht);
	    pmod->uhat[t] = mpf_get_d(uht);
	} else if (uvec != NULL) {
	    uvec->val[s] = mpf_get_d(uht);
	}
	if (uhat != NULL) {
	    mpf_set(uhat[s], uht);
	}
	s++;
    }

    mpf_clear(yht);
    mpf_clear(uht);
    mpf_clear(xbi);

    if (uhat != NULL) {
	mp_dwstat(mpmod, pmod, uhat);
	for (s=0; s<mpmod->nobs; s++) {
	    mpf_clear(uhat[s]);
	}
	free(uhat);
    } else if (pmod != NULL) {
	pmod->rho = pmod->dw = NADBL;
    }
}

static void mp_uncentered_r_squared (MODEL *pmod,
				     MPMODEL *mpmod,
				     const DATASET *dset)
{
    double y0 = dset->Z[pmod->list[1]][pmod->t1];

    /* special computation for the case of TSS = 0, i.e.
       the dependent variable is a constant */

    if (y0 != 0) {
	mpf_t my0, tmp, tmp2;

	/* TSS = n * y_0^2 */
	mpf_init_set_d(my0, y0);
	mpf_init(tmp);
	mpf_mul(tmp, my0, my0);
	mpf_init_set_d(tmp2, (double) pmod->nobs);
	mpf_mul(mpmod->tss, tmp, tmp2);

	/* R^2 = 1 - ESS/TSS */
	mpf_div(tmp, mpmod->ess, mpmod->tss);
	mpf_set_d(tmp2, 1.0);
	mpf_sub(mpmod->rsq, tmp2, tmp);

	pmod->rsq = mpf_get_d(mpmod->rsq);
	gretl_model_set_int(pmod, "uncentered", 1);

	mpf_clear(my0);
	mpf_clear(tmp);
	mpf_clear(tmp2);
    }
}

static int copy_mp_results (MPMODEL *mpmod, MODEL *pmod,
			    const DATASET *dset, mpf_t **mpZ,
			    char **xnames, int lhconst,
			    gretlopt opt)
{
    int tseries = dataset_is_time_series(dset);
    int i, err = 0;

    pmod->ncoeff = mpmod->ncoeff;
    pmod->full_n = dset->n;
    pmod->ci = MPOLS;

    err = gretl_model_allocate_storage(pmod);
    if (err) {
	if (xnames != NULL) {
	    strings_array_free(xnames, pmod->ncoeff);
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

    if (opt & OPT_X) {
	/* saving extra results */
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
	if (pmod->ifc == 0) {
	    /* do what we do with regular ols */
	    gretl_model_set_int(pmod, "uncentered", 1);
	    gretl_model_set_double(pmod, "centered-R2",
				   mpf_get_d(mpmod->adjrsq));
	    pmod->adjrsq = NADBL;
	}
    }

    if (lhconst) {
	mp_uncentered_r_squared(pmod, mpmod, dset);
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
	    gmp_fprintf(stderr, "mp_cholbeta: rtest = %Fg\n", rtest);
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

static void mp_regress (MPMODEL *pmod, MPXPXXPY xpxxpy,
			int *plhconst)
{
    int i, v, nobs, nv;
    mpf_t *diag = NULL;
    mpf_t ysum, ypy, zz, rss, tss;
    mpf_t den, sgmasq, tmp;
    double ess;
    MPCHOLBETA cb;
    int lhconst = 0;

    nv = xpxxpy.nv;

    pmod->sderr = malloc(nv * sizeof *pmod->sderr);
    if (pmod->sderr == NULL) {
        pmod->errcode = E_ALLOC;
        return;
    }

    for (i=0; i<nv; i++) {
	mpf_init(pmod->sderr[i]);
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

    if (mpf_sgn(tss) == 0) {
	lhconst = 1;
    } else if (mpf_sgn(tss) < 0) {
	fprintf(stderr, "mpols: TSS = %g\n", mpf_get_d(tss));
        pmod->errcode = E_TSS;
	goto cleanup;
    }

    /* Choleski-decompose X'X and find the coefficients */
    cb = mp_cholbeta(xpxxpy);
    pmod->coeff = cb.coeff;
    pmod->xpx = cb.xpxxpy.xpx;

    if (cb.errcode) {
        pmod->errcode = E_ALLOC;
        goto cleanup;
    }

    mpf_set(rss, cb.rss);
    mpf_clear(cb.rss);
    if (mpf_cmp(rss, MPF_MINUS_ONE) == 0) {
        pmod->errcode = E_SINGULAR;
        goto cleanup;
    }

    mpf_sub(pmod->ess, ypy, rss);
    ess = mpf_get_d(pmod->ess);
    if (fabs(ess) < DBL_EPSILON) {
	mpf_set(pmod->ess, MPF_ZERO);
    }

    if (mpf_sgn(pmod->ess) < 0) {
	gretl_errmsg_set(_("Error sum of squares is not >= 0"));
	pmod->errcode = E_DATA;
        goto cleanup;
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

    if (lhconst) {
	mpf_set_d(pmod->rsq, NADBL);
	mpf_set_d(pmod->adjrsq, NADBL);
    }

    if (pmod->errcode) {
	fprintf(stderr, "mp_ols: pmod->errcode = %d\n", pmod->errcode);
	goto cleanup;
    }

    if (mpf_sgn(tss) > 0) {
	/* compute SSR/TSS */
	mpf_div(tmp, pmod->ess, tss);
	/* and subtract from 1 to get R-squared */
	mpf_sub(pmod->rsq, MPF_ONE, tmp);
    }

    if (pmod->dfd > 0 && !lhconst) {
	if (pmod->ifc) {
	    /* dfd =  n - 1 */
	    mpf_set_d(tmp, (double) (nobs - 1));
	    mpf_div(tmp, tmp, den);
	    mpf_mul(tmp, tmp, pmod->ess);
	    mpf_sub(pmod->adjrsq, MPF_ONE, tmp);
	} else {
	    /* record the centered R-squared as "adjusted" */
	    mpf_set(pmod->adjrsq, pmod->rsq);
	    /* compute SSR / sum(y^2) */
	    mpf_div(tmp, pmod->ess, ypy);
	    /* reset rsq as uncentered R-squared */
	    mpf_sub(pmod->rsq, MPF_ONE, tmp);
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
	goto cleanup;
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

    if (plhconst != NULL) {
	*plhconst = lhconst;
    }

 cleanup:

    mpf_clear(den);
    mpf_clear(sgmasq);
    mpf_clear(ysum);
    mpf_clear(ypy);
    mpf_clear(zz);
    mpf_clear(rss);
    mpf_clear(tss);
    mpf_clear(tmp);
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

static int add_missvals_mask (MPMODEL *mpmod, const int *list,
			      const DATASET *dset)
{
    int s, t, i;

    mpmod->mask = calloc(dset->t2 - dset->t1 + 1, 1);
    if (mpmod->mask == NULL) {
	return E_ALLOC;
    }

    s = 0;
    for (t=dset->t1; t<=dset->t2; t++) {
	for (i=1; i<=list[0]; i++) {
	    if (na(dset->Z[list[i]][t])) {
		mpmod->mask[s] = 1;
		break;
	    }
	}
	s++;
    }

    return 0;
}

/* public interface */

/**
 * mplsq:
 * @list: dependent variable plus list of regressors.
 * @polylist: list of polynomial terms (or %NULL).
 * @zdigits: list of digits of input precision for data series
 * (or %NULL).
 * @dset: dataset struct.
 * @pmod: MODEL pointer to hold results.
 * @opt: if contains %OPT_X, save extra model
 * information (including the names of parameters in
 * @pmod, if required).
 *
 * Computes multiple-precision OLS estimates of the model
 * specified by @list, and stores them in @pmod.
 *
 * Returns: 0 on success, error code on failure.
 */

int mplsq (const int *list, const int *polylist, const int *zdigits,
	   DATASET *dset, MODEL *pmod, gretlopt opt)
{
    int l0, i;
    mpf_t **mpZ = NULL;
    char **xnames = NULL;
    MPXPXXPY xpxxpy;
    MPMODEL mpmod;
    int orig_t1 = dset->t1;
    int orig_t2 = dset->t2;
    int missvals = 0;
    int lhconst = 0;
    int err = 0;

    gretl_error_clear();

    if (list == NULL || dset == NULL || dset->Z == NULL ||
	list[0] < 2 || dset->v < 2) {
	return E_DATA;
    }

    real_set_mp_bits();

    mp_model_init(&mpmod);

    if (polylist != NULL) {
	err = poly_check(&mpmod, polylist, list);
	if (err) {
	    return err;
	}
    }

    if (polylist == NULL) {
	mpmod.list = gretl_list_copy(list);
    } else {
	mpmod.list = poly_copy_list(list, polylist);
    }

    if (mpmod.list == NULL) {
	return E_ALLOC;
    }

    l0 = mpmod.list[0];
    mpmod.ncoeff = l0 - 1;

    mpmod.polylist = polylist; /* attached for convenience */

    /* check for missing obs in sample */
    list_adjust_sample(list, &dset->t1, &dset->t2, dset, &missvals);
    mpmod.nobs = dset->t2 - dset->t1 + 1 - missvals;

    if (mpmod.nobs < mpmod.ncoeff) {
	err = E_DF;
	goto bailout;
    }

    if (missvals > 0) {
	err = add_missvals_mask(&mpmod, list, dset);
	if (err) {
	    goto bailout;
	}
    }

    mpmod.t1 = dset->t1;
    mpmod.t2 = dset->t2;

    /* check for other data issues */
    if (data_problems(list, dset)) {
	err = E_DATA;
	goto bailout;
    }

    /* enable names for polynomial terms? */
    if (polylist != NULL && (opt & OPT_X)) {
	xnames = allocate_xnames(mpmod.list);
	if (xnames == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    /* see if the regressor list contains a constant */
    mpmod.ifc = mp_rearrange(mpmod.list);

    /* construct multiple-precision data matrix */
    mpZ = make_mpZ(&mpmod, zdigits, dset, xnames);

    if (mpZ == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    mpf_constants_init();

    /* calculate regression results */

    xpxxpy = mp_xpxxpy_func(mpmod.list, mpmod.nobs, mpZ);
    mpf_set(mpmod.tss, xpxxpy.xpy[l0]);

    mp_regress(&mpmod, xpxxpy, &lhconst);

    for (i=0; i<=l0; i++) {
	mpf_clear(xpxxpy.xpy[i]);
    }
    free(xpxxpy.xpy);

    err = mpmod.errcode;
    if (!err) {
	err = copy_mp_results(&mpmod, pmod, dset, mpZ, xnames,
			      lhconst, opt);
    }

    /* free all the mpf stuff */
    mpf_2d_array_free(mpZ, l0, mpmod.nobs);
    mpf_constants_clear();

 bailout:

    dset->t1 = orig_t1;
    dset->t2 = orig_t2;

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

    real_set_mp_bits();
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
    mpf_2d_array_free(mpZ, l0, mpmod.nobs);
    mpf_constants_clear();

 bailout:

    mp_model_free(&mpmod);

    return err;
}

static mpf_t *doubles_array_to_mp (const double *dy, int n)
{
    mpf_t *y = NULL;
    int i;

    y = malloc(n * sizeof *y);
    if (y == NULL) {
	return NULL;
    }

    for (i=0; i<n; i++) {
	mpf_init_set_d(y[i], dy[i]);
    }

    return y;
}

static mpf_t *mpf_array_new (int n)
{
    mpf_t *y = NULL;
    int i;

    y = malloc(n * sizeof *y);
    if (y == NULL) {
	return NULL;
    }

    for (i=0; i<n; i++) {
	mpf_init(y[i]);
    }

    return y;
}

static void mpf_array_free (mpf_t *y, int n)
{
    if (y != NULL) {
	int i;

	for (i=0; i<n; i++) {
	    mpf_clear(y[i]);
	}
	free(y);
    }
}

static mpfr_t *doubles_array_to_mpfr (const double *dy, int n)
{
    mpfr_t *y = NULL;
    int i;

    y = malloc(n * sizeof *y);
    if (y == NULL) {
	return NULL;
    }

    for (i=0; i<n; i++) {
	mpfr_init_set_d(y[i], dy[i], GMP_RNDN);
    }

    return y;
}

static mpfr_t *mpfr_array_new (int n)
{
    mpfr_t *y = NULL;
    int i;

    y = malloc(n * sizeof *y);
    if (y == NULL) {
	return NULL;
    }

    for (i=0; i<n; i++) {
	mpfr_init_set_d(y[i], 0.0, GMP_RNDN);
    }

    return y;
}

static void mpfr_array_free (mpfr_t *y, int n)
{
    if (y != NULL) {
	int i;

	for (i=0; i<n; i++) {
	    mpfr_clear(y[i]);
	}
	free(y);
    }
}

static mpf_t **mpf_2d_array_new (int n, int T)
{
    mpf_t **x = NULL;
    int i, t;

    x = malloc(n * sizeof *x);
    if (x == NULL) {
	return NULL;
    }

    for (i=0; i<n; i++) {
	x[i] = NULL;
    }

    for (i=0; i<n; i++) {
	x[i] = malloc(T * sizeof **x);
	if (x[i] == NULL) {
	    mpf_2d_array_free(x, n, T);
	    return NULL;
	} else {
	    for (t=0; t<T; t++) {
		mpf_init(x[i][t]);
	    }
	}
    }

    return x;
}

static mpfr_t **mpfr_2d_array_new (int n, int T)
{
    mpfr_t **x = NULL;
    int i, t;

    x = malloc(n * sizeof *x);
    if (x == NULL) {
	return NULL;
    }

    for (i=0; i<n; i++) {
	x[i] = NULL;
    }

    for (i=0; i<n; i++) {
	x[i] = malloc(T * sizeof **x);
	if (x[i] == NULL) {
	    mpfr_2d_array_free(x, n, T);
	    return NULL;
	} else {
	    for (t=0; t<T; t++) {
		mpfr_init(x[i][t]);
	    }
	}
    }

    return x;
}

#define Min(x,y) (((x) > (y))? (y) : (x))

static int mp_symm_toeplitz (mpf_t *g, mpf_t *y, int T, int q)
{
    mpf_t **mu = NULL;
    mpf_t tmp;
    int t, j, k, jmax;
    int err = 0;

    mpf_init(tmp);

    mu = mpf_2d_array_new(q+1, T);
    if (mu == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* factorize */
    for (t=0; t<T; t++) {
	for (k=Min(q, t); k>=0; k--) {
	    mpf_set(mu[k][t], g[k]);
	    jmax = q - k;
	    for (j=1; j<=jmax && t-k-j >= 0; j++) {
		mpf_mul(tmp, mu[j][t-k], mu[j+k][t]);
		mpf_mul(tmp, tmp, mu[0][t-k-j]);
		mpf_sub(mu[k][t], mu[k][t], tmp);
	    }
	    if (k > 0) {
		mpf_div(mu[k][t], mu[k][t], mu[0][t-k]);
	    }
	}
    }

    /* forward solve */
    for (t=0; t<T; t++) {
	jmax = Min(t, q);
	for (j=1; j<=jmax; j++) {
	    mpf_mul(tmp, mu[j][t], y[t-j]);
	    mpf_sub(y[t], y[t], tmp);
	}
    }

    /* divide by the diagonal */
    for (t=0; t<T; t++) {
	mpf_div(y[t], y[t], mu[0][t]);
    }

    /* backsolve */
    for (t=T-1; t>=0; t--) {
	jmax = Min(q, T - 1 - t);
	for (j=1; j<=jmax; j++) {
	    mpf_mul(tmp, mu[j][t+j], y[t+j]);
	    mpf_sub(y[t], y[t], tmp);
	}
    }

 bailout:

    mpf_clear(tmp);
    mpf_2d_array_free(mu, q+1, T);

    return err;
}

static int mp_gamma_y (mpf_t *g, mpf_t *y, mpf_t *tmp, int T, int n)
{
    mpf_t lx, rx, lrsum, addbit;
    int t, j;
    int err = 0;

    mpf_init(lx);
    mpf_init(rx);
    mpf_init(lrsum);
    mpf_init(addbit);

    for (t=0; t<T; t++) {
	for (j=n; j>0; j--) {
	    mpf_set(tmp[j], tmp[j-1]);
	}
	mpf_mul(tmp[0], g[0], y[t]);
	for (j=1; j<=n; j++) {
	    mpf_set(lx, (t - j < 0)? MPF_ZERO : y[t-j]);
	    mpf_set(rx, (t + j >= T)? MPF_ZERO : y[t+j]);
	    mpf_add(lrsum, lx, rx);
	    mpf_mul(addbit, g[j], lrsum);
	    mpf_add(tmp[0], tmp[0], addbit);
	}
	if (t >= n) {
	    mpf_set(y[t-n], tmp[n]);
	}
    }

    for (j=0; j<n; j++) {
	mpf_set(y[T-j-1], tmp[j]);
    }

    mpf_clear(lx);
    mpf_clear(rx);
    mpf_clear(lrsum);
    mpf_clear(addbit);

    return err;
}

static void mp_form_gamma (mpf_t *g, mpf_t *mu, int q)
{
    mpf_t tmp;
    int j, k;

    mpf_init(tmp);

    for (j=0; j<=q; j++) {
	mpf_set_ui(g[j], 0);
	for (k=0; k<=q-j; k++) {
	    mpf_mul(tmp, mu[k], mu[k+j]);
	    mpf_add(g[j], g[j], tmp);
	}
    }

    mpf_clear(tmp);
}

static void mp_sum_or_diff (mpf_t *theta, int n, int sign)
{
    int j, q;

    mpf_set_ui(theta[0], 1);

    for (q=1; q<=n; q++) {
	mpf_set_ui(theta[q], 0);
	for (j=q; j>0; j--) {
	    if (sign > 0) {
		mpf_add(theta[j], theta[j], theta[j-1]);
	    } else {
		mpf_sub(theta[j], theta[j], theta[j-1]);
	    }
	}
    }
}

static void mp_form_mvec (mpf_t *g, mpf_t *mu, int n)
{
    mp_sum_or_diff(mu, n, 1);
    mp_form_gamma(g, mu, n);
}

static void mp_form_svec (mpf_t *g, mpf_t *mu, int n)
{
    mp_sum_or_diff(mu, n, -1);
    mp_form_gamma(g, mu, n);
}

/* g is the target, mu and tmp are used as workspace */

static void mp_form_wvec (mpf_t *g, mpf_t *mu, mpf_t *tmp,
			  int n, mpf_t *lambda)
{
    mpf_t x;
    int i;

    mpf_init(x);

    /* svec = Q'SQ where Q is the 2nd-diff operator */

    mp_form_svec(tmp, mu, n);
    for (i=0; i<=n; i++) {
	mpf_mul(g[i], lambda[0], tmp[i]);
    }

    mp_form_mvec(tmp, mu, n);
    for (i=0; i<=n; i++) {
	mpf_mul(x, lambda[1], tmp[i]);
	mpf_add(g[i], g[i], x);
    }

    mpf_clear(x);
}

static void mp_qprime_y (mpf_t *y, int T)
{
    mpf_t x;
    int t;

    mpf_init(x);

    for (t=0; t<T-2; t++) {
	mpf_add(y[t], y[t], y[t+2]);
	mpf_mul_ui(x, y[t+1], 2);
	mpf_sub(y[t], y[t], x);
    }

    mpf_clear(x);
}

static void mp_form_Qy (mpf_t *y, int T)
{
    mpf_t tmp, lag1, lag2, x;
    int t;

    mpf_init(tmp);
    mpf_init(lag1);
    mpf_init(lag2);
    mpf_init(x);

    for (t=0; t<T-2; t++) {
	mpf_set(tmp, y[t]);
	mpf_mul_ui(x, lag1, 2);
	mpf_sub(x, lag2, x);
	mpf_add(y[t], y[t], x); /* y[t] += lag2 - 2*lag1 */
	mpf_set(lag2, lag1);
	mpf_set(lag1, tmp);
    }

    mpf_mul_ui(x, lag1, 2);
    mpf_sub(y[T-2], lag2, x);
    mpf_set(y[T-1], lag1);

    mpf_clear(tmp);
    mpf_clear(lag1);
    mpf_clear(lag2);
    mpf_clear(x);
}

static void calc_lambda (int n, double cutoff, mpf_t *lambda)
{
    double dlam;

    dlam = 1 / tan(cutoff / 2);
    mpf_set_d(lambda[0], dlam);
    mpf_pow_ui(lambda[0], lambda[0], n * 2);
    mpf_set_ui(lambda[1], 1);

    /* there's really only one "lambda": one out of
       lam1, lam2 = 1 and has no effect on the
       calculation */

    if (mpf_cmp_ui(lambda[0], 1000000) > 0) {
	mpf_sqrt(lambda[0], lambda[0]);
	mpf_ui_div(lambda[1], 1, lambda[0]);
    }

    mpf_out_str(stderr, 10, 16, lambda[0]);
    fputc('\n', stderr);
    mpf_out_str(stderr, 10, 16, lambda[1]);
    fputc('\n', stderr);
}

/* note: for this function the cutoff is specified in radians */

int mp_bw_filter (const double *x, double *bw, int T, int n,
		  double cutoff)
{
    mpf_t *g, *ds, *tmp, *y;
    mpf_t mx, lambda[2];
    int t, m;
    int err = 0;

    real_set_mp_bits();
    mpf_constants_init();
    mpf_init(lambda[0]);
    mpf_init(lambda[1]);
    mpf_init(mx);

    m = 3 * (n+1);

    /* the workspace we need for everything except the
       Toeplitz solver */
    g = mpf_array_new(m);
    if (g == NULL) {
	return E_ALLOC;
    }

    ds = g + n + 1;
    tmp = ds + n + 1;

    calc_lambda(n, cutoff, lambda);

    /* copy the data into y */
    y = doubles_array_to_mp(x, T);
    if (y == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    mp_form_wvec(g, ds, tmp, n, lambda); /* W = M + lambda * Q'SQ */
    mp_qprime_y(y, T);

    /* solve (M + lambda*Q'SQ)x = d for x */
    err = mp_symm_toeplitz(g, y, T-2, n);

    if (!err) {
	mp_form_Qy(y, T);
	mp_form_svec(g, ds, n-2);
	mp_gamma_y(g, y, tmp, T, n-2);   /* Form SQg */
	/* write the trend into y (low-pass) */
	for (t=0; t<T; t++) {
	    mpf_mul(mx, lambda[0], y[t]);
	    bw[t] = x[t] - mpf_get_d(mx);
	}
    }

 bailout:

    mpf_clear(lambda[0]);
    mpf_clear(lambda[1]);
    mpf_clear(mx);

    mpf_constants_clear();

    mpf_array_free(g, m);
    mpf_array_free(y, T);

    return err;
}

/* Multiple precision versions of midas_weights() and
   midas_gradient(). We come here only if the standard
   double precision variants hit a range error -- see
   genfuncs.c for context.
*/

#define WDEBUG 0

int mp_midas_weights (const double *theta, int k,
		      gretl_matrix *w, int method)
{
    int p = gretl_vector_get_length(w);
    double eps = pow(2.0, -52);
    mpfr_t *mw, *mt;
    mpfr_t wsum, tmp;
    int i, j, err = 0;

    real_set_mpfr_bits();

    mw = mpfr_array_new(p);
    mt = doubles_array_to_mpfr(theta, k);

    if (mw == NULL || mt == NULL) {
	return E_ALLOC;
    }

    mpfr_init_set_d(wsum, 0.0, GMP_RNDN);
    mpfr_init(tmp);

    if (method == 1) {
	/* nealmon */
	mpfr_t incr;

	mpfr_init(incr);

	for (i=0; i<p; i++) {
	    mpfr_mul_ui(mw[i], mt[0], i+1, GMP_RNDN);
	    for (j=1; j<k; j++) {
		mpfr_ui_pow_ui(tmp, i+1, j+1, GMP_RNDN);
		mpfr_mul(incr, tmp, mt[j], GMP_RNDN);
		mpfr_add(mw[i], mw[i], incr, GMP_RNDN);
	    }
	    mpfr_set(tmp, mw[i], GMP_RNDN);
	    mpfr_exp(mw[i], tmp, GMP_RNDN);
	    mpfr_add(wsum, wsum, mw[i], GMP_RNDN);
	}

	mpfr_clear(incr);
    } else {
	/* beta */
	mpfr_t si, ai, bi;
	double dsi;

	mpfr_init(si);
	mpfr_init(ai);
	mpfr_init(bi);

	for (i=0; i<p; i++) {
	    dsi = i / (p - 1.0);
	    if (i == 0) {
		dsi += eps;
	    } else if (i == p-1) {
		dsi -= eps;
	    }
	    mpfr_set_d(si, dsi, GMP_RNDN);
	    mpfr_set_d(tmp, theta[0] - 1.0, GMP_RNDN);
	    mpfr_pow(ai, si, tmp, GMP_RNDN);
	    mpfr_set_d(si, 1.0 - dsi, GMP_RNDN);
	    mpfr_set_d(tmp, theta[1] - 1.0, GMP_RNDN);
	    mpfr_pow(bi, si, tmp, GMP_RNDN);
	    mpfr_mul(mw[i], ai, bi, GMP_RNDN);
	    mpfr_add(wsum, wsum, mw[i], GMP_RNDN);
	}

	mpfr_clear(si);
	mpfr_clear(ai);
	mpfr_clear(bi);
    }

    if (!err) {
	for (i=0; i<p; i++) {
	    /* normalize */
	    mpfr_div(tmp, mw[i], wsum, GMP_RNDN);
	}
	if (method == 3) {
	    /* beta with non-zero last lag: add theta[2]
	       and renormalize */
	    mpfr_set_d(wsum, 1.0 + p * theta[2], GMP_RNDN);
	    for (i=0; i<p; i++) {
		mpfr_add_d(mw[i], mw[i], theta[2], GMP_RNDN);
		mpfr_div(mw[i], mw[i], wsum, GMP_RNDN);
	    }
	}
	for (i=0; i<p; i++) {
	    w->val[i] = mpfr_get_d(tmp, GMP_RNDN);
	}
    }

#if WDEBUG
    gretl_matrix_print(w, "w, in mp module");
#endif

    mpfr_array_free(mw, p);
    mpfr_array_free(mt, k);
    mpfr_clear(wsum);
    mpfr_clear(tmp);

    mpfr_free_cache();

    return err;
}

int mp_midas_gradient (const double *theta,
		       gretl_matrix *G,
		       int method)
{
    int p = G->rows;
    int k = G->cols;
    double eps = pow(2.0, -52);
    mpfr_t *mw, *mt, **mg = NULL;
    mpfr_t wsum, tmp, gij;
    int i, j, err = 0;

    real_set_mpfr_bits();

    mw = mpfr_array_new(p);
    mt = doubles_array_to_mpfr(theta, k);

    if (mw == NULL || mt == NULL) {
	return E_ALLOC;
    }

    mpfr_init_set_d(wsum, 0.0, GMP_RNDN);
    mpfr_init(tmp);
    mpfr_init(gij);

    if (method > 1) {
	mg = mpfr_2d_array_new(k, p);
	if (mg == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    if (method == 1) {
	/* nealmon */
	mpfr_t *dsum = mpfr_array_new(k);
	double dgij;

	if (dsum == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}

	for (i=0; i<p; i++) {
	    mpfr_mul_ui(mw[i], mt[0], i+1, GMP_RNDN);
	    for (j=1; j<k; j++) {
		mpfr_ui_pow_ui(tmp, i+1, j+1, GMP_RNDN);
		mpfr_mul(tmp, tmp, mt[j], GMP_RNDN);
		mpfr_add(mw[i], mw[i], tmp, GMP_RNDN);
	    }
	    mpfr_set(tmp, mw[i], GMP_RNDN);
	    mpfr_exp(mw[i], tmp, GMP_RNDN);
	    mpfr_add(wsum, wsum, mw[i], GMP_RNDN);
	}
	for (i=0; i<p; i++) {
	    for (j=0; j<k; j++) {
		mpfr_ui_pow_ui(tmp, i+1, j+1, GMP_RNDN);
		mpfr_mul(tmp, tmp, mw[i], GMP_RNDN);
		mpfr_add(dsum[j], dsum[j], tmp, GMP_RNDN);
	    }
	}
	for (j=0; j<k; j++) {
	    mpfr_div(dsum[j], dsum[j], wsum, GMP_RNDN);
	}
	for (i=0; i<p; i++) {
	    mpfr_div(mw[i], mw[i], wsum, GMP_RNDN);
	    for (j=0; j<k; j++) {
		mpfr_ui_pow_ui(tmp, i+1, j+1, GMP_RNDN);
		mpfr_sub(tmp, tmp, dsum[j], GMP_RNDN);
		mpfr_mul(gij, mw[i], tmp, GMP_RNDN);
		dgij = mpfr_get_d(gij, GMP_RNDN);
		gretl_matrix_set(G, i, j, dgij);
	    }
	}

	mpfr_array_free(dsum, k);
    } else  {
	/* beta */
	mpfr_t si, ai, bi, ws2;
	mpfr_t g1sum, g2sum;
	double dsi, dgij;

	mpfr_init(si);
	mpfr_init(ai);
	mpfr_init(bi);
	mpfr_init(ws2);
	mpfr_init_set_d(g1sum, 0.0, GMP_RNDN);
	mpfr_init_set_d(g2sum, 0.0, GMP_RNDN);

	for (i=0; i<p; i++) {
	    dsi = i / (p - 1.0);
	    if (i == 0) {
		dsi += eps;
	    } else if (i == p - 1) {
		dsi -= eps;
	    }
	    mpfr_sub_ui(tmp, mt[0], 1, GMP_RNDN);
	    mpfr_set_d(si, dsi, GMP_RNDN);
	    mpfr_pow(ai, si, tmp, GMP_RNDN);
	    mpfr_sub_ui(tmp, mt[1], 1, GMP_RNDN);
	    mpfr_ui_sub(si, 1, si, GMP_RNDN);
	    mpfr_pow(bi, si, tmp, GMP_RNDN);
	    mpfr_mul(mw[i], ai, bi, GMP_RNDN);
	    mpfr_add(wsum, wsum, mw[i], GMP_RNDN);
	}
	mpfr_mul(ws2, wsum, wsum, GMP_RNDN);
	for (i=0; i<p; i++) {
	    dsi = i / (p - 1.0);
	    if (i == 0) {
		dsi += eps;
	    } else if (i == p - 1) {
		dsi -= eps;
	    }
	    mpfr_set_d(si, dsi, GMP_RNDN);
	    mpfr_log(tmp, si, GMP_RNDN);
	    mpfr_mul(ai, mw[i], tmp, GMP_RNDN);
	    mpfr_add(g1sum, g1sum, ai, GMP_RNDN);
	    mpfr_div(ai, ai, wsum, GMP_RNDN);
	    mpfr_set(mg[0][i], ai, GMP_RNDN);
	    mpfr_set_d(si, dsi, GMP_RNDN);
	    mpfr_ui_sub(si, 1, si, GMP_RNDN);
	    mpfr_log(tmp, si, GMP_RNDN);
	    mpfr_mul(bi, mw[i], tmp, GMP_RNDN);
	    mpfr_add(g2sum, g2sum, bi, GMP_RNDN);
	    mpfr_div(bi, bi, wsum, GMP_RNDN);
	    mpfr_set(mg[1][i], bi, GMP_RNDN);
	}
	for (i=0; i<p; i++) {
	    mpfr_set(ai, mg[0][i], GMP_RNDN);
	    mpfr_div(tmp, g1sum, ws2, GMP_RNDN);
	    mpfr_mul(tmp, tmp, mw[i], GMP_RNDN);
	    mpfr_sub(mg[0][i], ai, tmp, GMP_RNDN);
	    mpfr_set(bi, mg[1][i], GMP_RNDN);
	    mpfr_div(tmp, g2sum, ws2, GMP_RNDN);
	    mpfr_mul(tmp, tmp, mw[i], GMP_RNDN);
	    mpfr_sub(mg[1][i], bi, tmp, GMP_RNDN);
	}
	if (k == 3) {
	    /* not zero-terminated */
	    mpfr_t mm3;
	    double c3 = theta[2];
	    double m3 = 1 / (1 + p * c3);

	    mpfr_init_set_d(mm3, m3, GMP_RNDN);

	    /* scale the first two columns */
	    for (j=0; j<2; j++) {
		for (i=0; i<p; i++) {
		    mpfr_mul(mg[j][i], mg[j][i], mm3, GMP_RNDN);
		}
	    }
	    /* compute the third-col derivatives */
	    for (i=0; i<p; i++) {
		mpfr_div(tmp, mw[i], wsum, GMP_RNDN);
		mpfr_mul_ui(tmp, tmp, p, GMP_RNDN);
		mpfr_ui_sub(tmp, 1, tmp, GMP_RNDN);
		mpfr_mul(tmp, tmp, mm3, GMP_RNDN);
		mpfr_mul(mg[2][i], tmp, mm3, GMP_RNDN);
	    }

	    mpfr_clear(mm3);
	}

	/* transcribe results */
	for (j=0; j<k; j++) {
	    for (i=0; i<p; i++) {
		dgij = mpfr_get_d(mg[j][i], GMP_RNDN);
		gretl_matrix_set(G, i, j, dgij);
	    }
	}

	mpfr_clear(si);
	mpfr_clear(ai);
	mpfr_clear(bi);
	mpfr_clear(ws2);
	mpfr_clear(g1sum);
	mpfr_clear(g2sum);
    }

#if WDEBUG
    gretl_matrix_print(G, "G, in mp module");
#endif

 bailout:

    mpfr_array_free(mw, p);
    mpfr_array_free(mt, k);
    mpfr_2d_array_free(mg, k, p);
    mpfr_clear(wsum);
    mpfr_clear(tmp);
    mpfr_clear(gij);

    mpfr_free_cache();

    return err;
}
