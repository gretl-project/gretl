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

/* Native gretl code for ARMA estimation.  Much of the code here
   was contributed by Riccardo "Jack" Lucchetti, the rest is due
   to Allin Cottrell; thanks also to Stephen Moshier for cephes.
*/

#include "libgretl.h"
#include "bhhh_max.h"
#include "libset.h"
#include "kalman.h"
#include "matrix_extra.h"
#include "../cephes/libprob.h"

#define ARMA_DEBUG 0
#define AINIT_DEBUG 0

/* ln(sqrt(2*pi)) + 0.5 */
#define LN_SQRT_2_PI_P5 1.41893853320467274178

#define KALMAN_ALL 999

#include "arma_common.c"

static PRN *vprn;

struct bchecker {
    int qmax;
    double *temp;
    double *tmp2;
    cmplx *roots;
};

static void bchecker_free (struct bchecker *b)
{
    if (b != NULL) {
	free(b->temp);
	free(b->tmp2);
	free(b->roots);
	free(b);
    }
}

static struct bchecker *bchecker_allocate (struct arma_info *ainfo)
{
    struct bchecker *b;

    b = malloc(sizeof *b);
    if (b == NULL) {
	return NULL;
    }

    b->temp = NULL;
    b->tmp2 = NULL;
    b->roots = NULL;

    b->qmax = ainfo->q + ainfo->Q * ainfo->pd;

    b->temp  = malloc((b->qmax + 1) * sizeof *b->temp);
    b->tmp2  = malloc((b->qmax + 1) * sizeof *b->tmp2);
    b->roots = malloc(b->qmax * sizeof *b->roots);

    if (b->temp == NULL || b->tmp2 == NULL || b->roots == NULL) {
	bchecker_free(b);
	b = NULL;
    } 

    return b;
}

/* check whether the MA estimates have gone out of bounds in the
   course of iteration */

static int 
ma_out_of_bounds (struct arma_info *ainfo, const double *theta,
		  const double *Theta)
{
    static struct bchecker *b = NULL;
    double re, im, rt;
    int i, j, k, m, si, qtot;
    int tzero = 1, Tzero = 1;
    int err = 0, cerr = 0;

    if (ainfo == NULL) {
	/* signal for cleanup */
	bchecker_free(b);
	b = NULL;
	return 0;
    }

    k = 0;
    for (i=0; i<ainfo->q && tzero; i++) {
	if (MA_included(ainfo, i)) {
	    if (theta[k++] != 0.0) {
		tzero = 0;
	    }
	}
    }  

    for (i=0; i<ainfo->Q && Tzero; i++) {
	if (Theta[i] != 0.0) {
	    Tzero = 0;
	}    
    }  
    
    if (tzero && Tzero) {
	/* nothing to be done */
	return 0;
    }

    if (b == NULL) {
	b = bchecker_allocate(ainfo);
	if (b == NULL) {
	    return 1;
	}
    }

    b->temp[0] = 1.0;

    /* initialize to non-seasonal MA or zero */
    k = 0;
    for (i=0; i<b->qmax; i++) {
        if (i < ainfo->q && MA_included(ainfo, i)) {
	    b->temp[i+1] = theta[k++];
        } else {
            b->temp[i+1] = 0.0;
        }
    }

    /* add seasonal MA and interaction */
    if (Tzero) {
	qtot = ainfo->q;
    } else {
	qtot = b->qmax;
	for (j=0; j<ainfo->Q; j++) {
	    si = (j + 1) * ainfo->pd;
	    b->temp[si] += Theta[j];
	    k = 0;
	    for (i=0; i<ainfo->q; i++) {
		if (MA_included(ainfo, i)) {
		    m = si + (i + 1);
		    b->temp[m] += Theta[j] * theta[k++];
		} 
	    }
	}
    }

#if ARMA_DEBUG
    if (b->temp[qtot] == 0.0) {
	fprintf(stderr, "b->temp[%d] = 0; polrt won't work\n", qtot);
	fprintf(stderr, "q = %d, Q = %d, b->qmax = %d\n", 
		ainfo->q, ainfo->Q, b->qmax);
	for (i=0; i<=qtot; i++) {
	    fprintf(stderr, "b->temp[%d] = %g\n", i, b->temp[i]);
	}
    }
#endif

    cerr = polrt(b->temp, b->tmp2, qtot, b->roots);
    if (cerr) {
	fprintf(stderr, "ma_out_of_bounds: polrt returned %d\n", cerr);
	return 0; /* ?? */
    }

    for (i=0; i<qtot; i++) {
	re = b->roots[i].r;
	im = b->roots[i].i;
	rt = re * re + im * im;
	if (rt > DBL_EPSILON && rt <= 1.0) {
	    pprintf(vprn, "MA root %d = %g\n", i, rt);
	    err = 1;
	    break;
	}
    }

    return err;
}

static void bounds_checker_cleanup (void)
{
    ma_out_of_bounds(NULL, NULL, NULL);
}

static void do_MA_partials (double *drv,
			    struct arma_info *ainfo,
			    const double *theta,
			    const double *Theta,
			    int t)
{
    int i, j, k, s, p;

    k = 0;
    for (i=0; i<ainfo->q; i++) {
	if (MA_included(ainfo, i)) {
	    s = t - (i + 1);
	    if (s >= 0) {
		drv[t] -= theta[k] * drv[s];
	    }
	    k++;
	}
    }

    for (j=0; j<ainfo->Q; j++) {
	s = t - (j + 1) * ainfo->pd;
	if (s >= 0) {
	    drv[t] -= Theta[j] * drv[s];
	    k = 0;
	    for (i=0; i<ainfo->q; i++) {
		if (MA_included(ainfo, i)) {
		    p = s - (i + 1);
		    if (p >= 0) {
			drv[t] -= Theta[j] * theta[k] * drv[p];
		    }
		    k++;
		}
	    }
	}
    }
}

/* Calculate ARMA log-likelihood.  This function is passed to the
   bhhh_max() routine as a callback. */

static int arma_ll (double *coeff, 
		    const double **bhX, double **Z, 
		    model_info *arma,
		    int do_score)
{
    int i, j, k, s, t;
    int t1 = model_info_get_t1(arma);
    int t2 = model_info_get_t2(arma);
    int n = t2 - t1 + 1;

    const double *y = bhX[0];
    const double **X = bhX + 1;
    double **series = model_info_get_series(arma);
    double *e = series[0];
    double **de = series + 1;
    double **de_a, **de_sa, **de_m, **de_sm, **de_r;
    const double *phi, *Phi;
    const double *theta, *Theta;
    const double *beta;

    struct arma_info *ainfo;

    double ll, s2 = 0.0;
    int err = 0;

    /* retrieve ARMA-specific information */
    ainfo = model_info_get_extra_info(arma);

    /* pointers to blocks of coefficients */
    phi = coeff + ainfo->ifc;
    Phi = phi + ainfo->np;
    theta = Phi + ainfo->P;
    Theta = theta + ainfo->nq;
    beta = Theta + ainfo->Q;

    /* pointers to blocks of derivatives */
    de_a = de + ainfo->ifc;
    de_sa = de_a + ainfo->np;
    de_m = de_sa + ainfo->P;
    de_sm = de_m + ainfo->nq;
    de_r = de_sm + ainfo->Q;

#if ARMA_DEBUG
    fprintf(stderr, "arma_ll: p=%d, q=%d, P=%d, Q=%d, pd=%d\n",
	    ainfo->p, ainfo->q, ainfo->P, ainfo->Q, ainfo->pd);
#endif

    if (ma_out_of_bounds(ainfo, theta, Theta)) {
	pputs(vprn, "arma: MA estimate(s) out of bounds\n");
	fputs("arma: MA estimate(s) out of bounds\n", stderr);
	return 1;
    }

    /* update forecast errors */

    for (t=t1; t<=t2; t++) {
	int p;

	e[t] = y[t];

	/* intercept */
	if (ainfo->ifc) {
	    e[t] -= coeff[0];
	} 

	/* non-seasonal AR component */
	k = 0;
	for (i=0; i<ainfo->p; i++) {
	    if (AR_included(ainfo, i)) {
		s = t - (i + 1);
		e[t] -= phi[k++] * y[s];
	    }
	}

	/* seasonal AR component plus interactions */
	for (j=0; j<ainfo->P; j++) {
	    s = t - (j + 1) * ainfo->pd;
	    e[t] -= Phi[j] * y[s];
	    k = 0;
	    for (i=0; i<ainfo->p; i++) {
		if (AR_included(ainfo, i)) {
		    p = s - (i + 1);
		    e[t] += Phi[j] * phi[k++] * y[p];
		}
	    }
	}

	/* non-seasonal MA component */
	k = 0;
	for (i=0; i<ainfo->q; i++) {
	    if (MA_included(ainfo, i)) {
		s = t - (i + 1);
		if (s >= t1) {
		    e[t] -= theta[k] * e[s];
		}
		k++;
	    }
	}

	/* seasonal MA component plus interactions */
	for (j=0; j<ainfo->Q; j++) {
	    s = t - (j + 1) * ainfo->pd;
	    if (s >= t1) {
		e[t] -= Theta[j] * e[s];
		k = 0;
		for (i=0; i<ainfo->q; i++) {
		    if (MA_included(ainfo, i)) {
			p = s - (i + 1);
			if (p >= t1) {
			    e[t] -= Theta[j] * theta[k] * e[p];
			}
			k++;
		    }
		}
	    }
	}

	/* exogenous regressors */
	for (i=0; i<ainfo->nexo; i++) {
	    e[t] -= beta[i] * X[i][t];
	}

	s2 += e[t] * e[t];
    }

    /* get error variance and log-likelihood */

    s2 /= (double) n;

    ll = -n * (0.5 * log(s2) + LN_SQRT_2_PI_P5);
    model_info_set_ll(arma, ll, do_score);

    if (do_score) {
	int lag, xlag;
	double x;

	for (t=t1; t<=t2; t++) {

	    /* the constant term (de_0) */
	    if (ainfo->ifc) {
		de[0][t] = -1.0;
		do_MA_partials(de[0], ainfo, theta, Theta, t);
	    }

	    /* non-seasonal AR terms (de_a) */
	    k = 0;
	    for (i=0; i<ainfo->p; i++) {
		if (!AR_included(ainfo, i)) {
		    continue;
		}
		lag = i + 1;
		if (t >= lag) {
		    de_a[k][t] = -y[t-lag];
		    /* cross-partial with seasonal AR */
		    for (j=0; j<ainfo->P; j++) {
			xlag = lag + (j + 1) * ainfo->pd;
			if (t >= xlag) {
			    de_a[k][t] += Phi[j] * y[t-xlag];
			}
		    }
		    do_MA_partials(de_a[k], ainfo, theta, Theta, t);
		}
		k++;
	    }

	    /* seasonal AR terms (de_sa) */
	    for (j=0; j<ainfo->P; j++) {
		lag = (j + 1) * ainfo->pd;
		if (t >= lag) {
		    de_sa[j][t] = -y[t-lag];
		    /* cross-partial with non-seasonal AR */
		    k = 0;
		    for (i=0; i<ainfo->p; i++) {
			if (AR_included(ainfo, i)) {
			    xlag = lag + (i + 1);
			    if (t >= xlag) {
				de_sa[j][t] += phi[k] * y[t-xlag];
			    }
			    k++;
			}
		    }
		    do_MA_partials(de_sa[j], ainfo, theta, Theta, t);
		}
	    }

	    /* non-seasonal MA terms (de_m) */
	    k = 0;
	    for (i=0; i<ainfo->q; i++) {
		if (!MA_included(ainfo, i)) {
		    continue;
		}
		lag = i + 1;
		if (t >= lag) {
		    de_m[k][t] = -e[t-lag];
		    /* cross-partial with seasonal MA */
		    for (j=0; j<ainfo->Q; j++) {
			xlag = lag + (j + 1) * ainfo->pd;
			if (t >= xlag) {
			    de_m[k][t] -= Theta[j] * e[t-xlag];
			}
		    }
		    do_MA_partials(de_m[k], ainfo, theta, Theta, t);
		}
		k++;
	    }

	    /* seasonal MA terms (de_sm) */
	    for (j=0; j<ainfo->Q; j++) {
		lag = (j + 1) * ainfo->pd;
		if (t >= lag) {
		    de_sm[j][t] = -e[t-lag];
		    /* cross-partial with non-seasonal MA */
		    k = 0;
		    for (i=0; i<ainfo->q; i++) {
			if (MA_included(ainfo, i)) {
			    xlag = lag + (i + 1);
			    if (t >= xlag) {
				de_sm[j][t] -= theta[k] * e[t-xlag];
			    }
			    k++;
			}
		    }
		    do_MA_partials(de_sm[j], ainfo, theta, Theta, t);
		}
	    }

	    /* exogenous regressors (de_r) */
	    for (j=0; j<ainfo->nexo; j++) {
		de_r[j][t] = -X[j][t]; 
		do_MA_partials(de_r[j], ainfo, theta, Theta, t);
	    }

	    /* update OPG data set */
	    x = e[t] / s2; /* sqrt(s2)? does it matter? */
	    for (i=0; i<ainfo->nc; i++) {
		Z[i+1][t] = -de[i][t] * x;
	    }
	}
    }

    if (isnan(ll)) {
	err = 1;
    }

    return err;
}

/*
  Given an ARMA process $A(L)B(L) y_t = C(L)D(L) \epsilon_t$, finds the 
  roots of the four polynomials -- or just two polynomials if seasonal
  AR and MA effects, B(L) and D(L) are not present -- and attaches
  this information to the ARMA model.

  pmod: MODEL pointer to which the roots info should be attached.

  ainfo: gives various pieces of information on the ARMA model,
  including seasonal and non-seasonal AR and MA orders.

  coeff: ifc + p + q + P + Q vector of coefficients (if an intercept
  is present it is element 0 and is ignored)

  returns: zero on success, non-zero on failure
*/

static int arma_model_add_roots (MODEL *pmod, struct arma_info *ainfo,
				 const double *coeff)
{
    const double *phi = coeff + ainfo->ifc;
    const double *Phi = phi + ainfo->np;
    const double *theta = Phi + ainfo->P;
    const double *Theta = theta + ainfo->nq;

    int nr = ainfo->p + ainfo->P + ainfo->q + ainfo->Q;
    int pmax, qmax, lmax;
    double *temp = NULL, *tmp2 = NULL;
    cmplx *rptr, *roots = NULL;
    int i, k;

    pmax = (ainfo->p > ainfo->P)? ainfo->p : ainfo->P;
    qmax = (ainfo->q > ainfo->Q)? ainfo->q : ainfo->Q;
    lmax = (pmax > qmax)? pmax : qmax;

    if (pmax == 0 && qmax == 0) {
	return 0;
    }

    temp = malloc((lmax + 1) * sizeof *temp);
    tmp2 = malloc((lmax + 1) * sizeof *tmp2);
    roots = malloc(nr * sizeof *roots);

    if (temp == NULL || tmp2 == NULL || roots == NULL) {
	free(temp);
	free(tmp2);
	free(roots);
	return E_ALLOC;
    }

    temp[0] = 1.0;
    rptr = roots;

    if (ainfo->p > 0) {
	/* A(L), non-seasonal */
	k = 0;
	for (i=0; i<ainfo->p; i++) {
	    if (AR_included(ainfo, i)) {
		temp[i+1] = -phi[k++];
	    } else {
		temp[i+1] = 0;
	    }
	}
	polrt(temp, tmp2, ainfo->p, rptr);
	rptr += ainfo->p;
    }

    if (ainfo->P > 0) {
	/* B(L), seasonal */
	for (i=0; i<ainfo->P; i++) {
	    temp[i+1] = -Phi[i];
	}    
	polrt(temp, tmp2, ainfo->P, rptr);
	rptr += ainfo->P;
    }

    if (ainfo->q > 0) {
	/* C(L), non-seasonal */
	k = 0;
	for (i=0; i<ainfo->q; i++) {
	    if (MA_included(ainfo, i)) {
		temp[i+1] = theta[k++];
	    } else {
		temp[i+1] = 0;
	    }
	}  
	polrt(temp, tmp2, ainfo->q, rptr);
	rptr += ainfo->q;
    }

    if (ainfo->Q > 0) {
	/* D(L), seasonal */
	for (i=0; i<ainfo->Q; i++) {
	    temp[i+1] = Theta[i];
	}  
	polrt(temp, tmp2, ainfo->Q, rptr);
    }
    
    free(temp);
    free(tmp2);

    gretl_model_set_data(pmod, "roots", roots, GRETL_TYPE_CMPLX_ARRAY,
			 nr * sizeof *roots);

    return 0;
}

/* below: exact ML using Kalman filter apparatus */

static gretl_matrix *S = NULL;
static gretl_matrix *P = NULL;
static gretl_matrix *F = NULL;
static gretl_matrix *A = NULL;
static gretl_matrix *H = NULL;
static gretl_matrix *Q = NULL;
static gretl_matrix *E = NULL;

static gretl_matrix *Svar;
static gretl_matrix *Svar2;
static gretl_matrix *vQ;

static double *ac;
static double *mc;

static struct arma_info *kainfo;    

static int ainfo_get_r (struct arma_info *ainfo)
{
    int pmax = ainfo->p + ainfo->pd * ainfo->P;
    int qmax = ainfo->q + ainfo->pd * ainfo->Q;

    return (pmax > qmax + 1)? pmax : qmax + 1;
}

static int allocate_ac_mc (struct arma_info *ainfo)
{
    if (ainfo->P > 0) {
	int pmax = ainfo->p + ainfo->pd * ainfo->P;

	ac = malloc((pmax + 1) * sizeof *ac);
	if (ac == NULL) {
	    return E_ALLOC;
	}
    }

    if (ainfo->Q > 0) {
	int qmax = ainfo->q + ainfo->pd * ainfo->Q;

	mc = malloc((qmax + 1) * sizeof *mc);
	if (mc == NULL) {
	    return E_ALLOC;
	}
    }

    return 0;
}

static void free_ac_mc (void)
{
    if (ac != NULL) free(ac);
    if (mc != NULL) free(mc);
}

static void write_big_phi (const double *phi, 
			   const double *Phi,
			   struct arma_info *ainfo,
			   gretl_matrix *F)
{
    int pmax = ainfo->p + ainfo->pd * ainfo->P;
    double x, y;
    int i, j, k, ii;

    for (i=0; i<=pmax; i++) {
	ac[i] = 0.0;
    }

    for (j=-1; j<ainfo->P; j++) {
	x = (j < 0)? -1 : Phi[j];
	k = 0.0;
	for (i=-1; i<ainfo->p; i++) {
	    if (i < 0) {
		y = -1;
	    } else if (AR_included(ainfo, i)) {
		y = phi[k++];
	    } else {
		y = 0.0;
	    }
	    ii = (j+1) * ainfo->pd + (i+1);
	    ac[ii] -= x * y;
	}
    }

    for (i=0; i<pmax; i++) {
	gretl_matrix_set(F, 0, i, ac[i+1]);
    }
}

static void write_big_theta (const double *theta, 
			     const double *Theta,
			     struct arma_info *ainfo,
			     gretl_matrix *H)
{
    int qmax = ainfo->q + ainfo->pd * ainfo->Q;
    double x, y;
    int i, j, k, ii;

    for (i=0; i<=qmax; i++) {
	mc[i] = 0.0;
    }

    for (j=-1; j<ainfo->Q; j++) {
	x = (j < 0)? 1 : Theta[j];
	k = 0;
        for (i=-1; i<ainfo->q; i++) {
	    if (i < 0) {
		y = 1;
	    } else if (MA_included(ainfo, i)) {
		y = theta[k++];
	    } else {
		y = 0;
	    }
            ii = (j+1) * ainfo->pd + (i+1);
	    mc[ii] = x * y;
        }
    }

    for (i=1; i<=qmax; i++) {
	H->val[i] = mc[i];
    }    
}

static void condense_row (gretl_matrix *targ,
			  const gretl_matrix *src,
			  int targrow, int srcrow,
			  int n)
{
    double x;
    int i, j, k, g;
    int targcol = 0;

    for (j=0; j<n; j++) {
	for (i=j; i<n; i++) {
	    k = j * n + i;
	    g = (k % n) * n + k / n;
	    x = gretl_matrix_get(src, srcrow, k);
	    if (g != k) {
		x += gretl_matrix_get(src, srcrow, g);
	    } 
	    gretl_matrix_set(targ, targrow, targcol++, x);
	}
    }
}

static void condense_state_vcv (gretl_matrix *targ, 
				const gretl_matrix *src,
				int n)
{
    int posr = 0, posc = 0;
    int i, j;

    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    if (j >= i) {
		condense_row(targ, src, posr++, posc, n);
	    }
	    posc++;
	}
    }
}

static void kalman_matrices_init (struct arma_info *ainfo)
{
    int r = F->rows;

    gretl_matrix_zero(A);
    gretl_matrix_zero(P);

    gretl_matrix_zero(F);
    gretl_matrix_inscribe_I(F, 1, 0, r - 1);

    gretl_matrix_zero(Q);
    gretl_matrix_set(Q, 0, 0, 1.0);

    gretl_matrix_zero(H);
    gretl_vector_set(H, 0, 1.0);
}

#define ARMA_MDEBUG 0

static int write_kalman_matrices (const double *b, int idx)
{
    const double *phi = b + kainfo->ifc;
    const double *Phi = phi + kainfo->np;
    const double *theta = Phi + kainfo->P;
    const double *Theta = theta + kainfo->nq;
    const double *beta = Theta + kainfo->Q;
    double mu = (kainfo->ifc)? b[0] : 0.0;
    int i, k, err = 0;

    int rewrite_A = 0;
    int rewrite_F = 0;
    int rewrite_H = 0;

    gretl_matrix_zero(S);

    if (idx == KALMAN_ALL) {
	rewrite_A = rewrite_F = rewrite_H = 1;
    } else {
	/* called in context of calculating score */
	int pmax = kainfo->ifc + kainfo->np + kainfo->P;
	int tmax = pmax + kainfo->nq + kainfo->Q;

	if (kainfo->ifc && idx == 0) {
	    rewrite_A = 1;
	} else if (idx >= kainfo->ifc && idx < pmax) {
	    rewrite_F = 1;
	} else if (idx >= kainfo->ifc && idx < tmax) {
	    rewrite_H = 1;
	} else {
	    rewrite_A = 1;
	}
    }

#if ARMA_MDEBUG
    fprintf(stderr, "\n*** write_kalman_matrices: before\n");
    gretl_matrix_print(A, "A");
    gretl_matrix_print(P, "P");
    gretl_matrix_print(F, "F");
    gretl_matrix_print(H, "H");
#endif    

    /* See Hamilton, Time Series Analysis, ch 13, p. 375 */

    if (rewrite_A) {
	/* const and coeffs on exogenous vars */
	gretl_vector_set(A, 0, mu);
	for (i=0; i<kainfo->nexo; i++) {
	    gretl_vector_set(A, i + 1, beta[i]);
	}
    }

    if (rewrite_H) {
	/* form the H vector using theta and/or Theta */
	if (kainfo->Q > 0) {
	    write_big_theta(theta, Theta, kainfo, H);
	} else {
	    k = 0;
	    for (i=0; i<kainfo->q; i++) {
		if (MA_included(kainfo, i)) {
		    gretl_vector_set(H, i+1, theta[k++]);
		} else {
		    gretl_vector_set(H, i+1, 0.0);
		}
	    }
	}
    }

    if (rewrite_F) {
	/* form the F matrix using phi and/or Phi */
	if (kainfo->P > 0) {
	    write_big_phi(phi, Phi, kainfo, F);
	} else {
	    k = 0;
	    for (i=0; i<kainfo->p; i++) {
		if (AR_included(kainfo, i)) {
		    gretl_matrix_set(F, 0, i, phi[k++]);
		} else {
		    gretl_matrix_set(F, 0, i, 0.0);
		}
	    }
	} 

	/* form $P_{1|0}$ (MSE) matrix, as per Hamilton, ch 13, p. 378. */

	gretl_matrix_kronecker_product(F, F, Svar);
	gretl_matrix_I_minus(Svar);

	if (arma_using_vech(kainfo)) {
	    condense_state_vcv(Svar2, Svar, gretl_matrix_rows(F));
	    gretl_matrix_vectorize_h(vQ, Q);
	    err = gretl_LU_solve(Svar2, vQ);
	    if (!err) {
		gretl_matrix_unvectorize_h(P, vQ);
	    }
	} else {
	    gretl_matrix_vectorize(vQ, Q);
	    err = gretl_LU_solve(Svar, vQ);
	    if (!err) {
		gretl_matrix_unvectorize(P, vQ);
	    }
	}
    }

#if ARMA_MDEBUG
    fprintf(stderr, "\n*** after\n");
    gretl_matrix_print(A, "A");
    gretl_matrix_print(P, "P");
    gretl_matrix_print(F, "F");
    gretl_matrix_print(H, "H");
#endif    

    return err;
}

static int rewrite_kalman_matrices (kalman *K, const double *b, int i)
{
    int err = write_kalman_matrices(b, i);

    if (!err) {
	kalman_set_initial_state_vector(K, S);
	kalman_set_initial_MSE_matrix(K, P);
    }

    return err;
}

/* add covariance matrix and standard errors based on numerical
   approximation to the Hessian
*/

static void arma_hessian_vcv (MODEL *pmod, double *vcv, int k)
{
    double x;
    int t, i = 0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->uhat[t] = gretl_vector_get(E, i++);
    }

    pmod->vcv = vcv;

    for (i=0; i<k; i++) {
	x = vcv[ijton(i, i, k)];
	pmod->sderr[i] = (na(x))? NADBL : sqrt(x);
    }
}

static double *
kalman_arma_score_callback (const double *b, int i, void *data)
{
    kalman *K = (kalman *) data;
    int err;

    rewrite_kalman_matrices(K, b, i);
    err = kalman_forecast(K);

#if ARMA_DEBUG
    fprintf(stderr, "kalman_arma_score: kalman f'cast gave "
	    "err = %d, ll = %#.12g\n", err, kalman_get_loglik(K));
#endif

    return (err)? NULL : E->val;
}

/* add covariance matrix and standard errors based on Outer Product of
   Gradient
*/

static int arma_OPG_vcv (MODEL *pmod, kalman *K, double *b, 
			 double s2, int k, int T)
{
    gretl_matrix *G = NULL;
    gretl_matrix *V = NULL;
    double x;
    int i, j, s, t;
    int idx, err = 0;

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->uhat[t] = gretl_vector_get(E, s++);
    }

    G = build_OPG_matrix(b, k, T, kalman_arma_score_callback, (void *) K, &err);
    if (err) {
	goto bailout;
    }

    V = gretl_matrix_alloc(k, k);
    if (V == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    gretl_matrix_multiply_mod(G, GRETL_MOD_NONE,
			      G, GRETL_MOD_TRANSPOSE,
			      V, GRETL_MOD_NONE);

    err = gretl_invert_symmetric_matrix(V);

    if (!err) {
	for (i=0; i<k; i++) {
	    for (j=0; j<=i; j++) {
		idx = ijton(i, j, k);
		x = s2 * gretl_matrix_get(V, i, j);
		pmod->vcv[idx] = x;
		if (i == j) {
		    pmod->sderr[i] = sqrt(x);
		}
	    }
	}
    }

 bailout:

    pmod->errcode = err;

    gretl_matrix_free(G);
    gretl_matrix_free(V);
    
    return err;
}

/* in Kalman case the basic model struct is empty, so we have
   to allocate for coefficients, residuals and so on */

static int kalman_arma_model_allocate (MODEL *pmod, int k, int T,
				       double *vcv)
{
    int err = 0;

    pmod->ncoeff = k;
    pmod->full_n = T;

    err = gretl_model_allocate_storage(pmod);

    if (!err && vcv == NULL) {
	err = gretl_model_new_vcv(pmod, NULL);
    }

    return err;
}

static int kalman_arma_finish (MODEL *pmod, const int *alist,
			       struct arma_info *ainfo,
			       const double **Z, const DATAINFO *pdinfo, 
			       kalman *K, double *b, double *vcv,
			       int k, int T)
{
    double s2;
    int kopt, i;

    pmod->t1 = ainfo->t1;
    pmod->t2 = ainfo->t2;
    pmod->nobs = T;

    pmod->errcode = kalman_arma_model_allocate(pmod, k, ainfo->T, vcv);
    if (pmod->errcode) {
	return pmod->errcode;
    }

    for (i=0; i<k; i++) {
	pmod->coeff[i] = b[i];
    }

    s2 = kalman_get_arma_variance(K);
    pmod->sigma = sqrt(s2);

    pmod->lnL = kalman_get_loglik(K);
    kopt = kalman_get_options(K);

    /* rescale stuff if we are using average loglikelihood */

    if (kopt & KALMAN_AVG_LL) {
	pmod->lnL *= T;
	if (vcv != NULL) {
	    int k2 = k * (k + 1) / 2;

	    for (i=0; i<k2; i++) {
		vcv[i] /= T;
	    }
	}
    }

#if ARMA_DEBUG
    fprintf(stderr, "kalman_arma_finish: doing VCV, method %s\n",
	    (vcv != NULL)? "Hessian" : "OPG");
#endif

    if (vcv != NULL) {
	arma_hessian_vcv(pmod, vcv, k);
	gretl_model_set_int(pmod, "ml_vcv", VCV_HESSIAN);
    } else {
	arma_OPG_vcv(pmod, K, b, s2, k, T);
	gretl_model_set_int(pmod, "ml_vcv", VCV_OP);
    }

    write_arma_model_stats(pmod, alist, ainfo, Z, pdinfo);
    arma_model_add_roots(pmod, ainfo, b);

    gretl_model_set_int(pmod, "arma_flags", ARMA_EXACT);

    return pmod->errcode;
}

#if ARMA_DEBUG

static void debug_print_theta (const double *theta,
			       const double *Theta)
{
    int i, k = 0;

    fprintf(stderr, "kalman_arma_ll():\n");

    for (i=0; i<kainfo->q; i++) {
	if (MA_included(kainfo, i)) {
	    fprintf(stderr, "theta[%d] = %#.12g\n", i+1, theta[k++]);
	}
    }

    for (i=0; i<kainfo->Q; i++) {
	fprintf(stderr, "Theta[%d] = %#.12g\n", i, Theta[i]);
    }   
}

#endif

static int kalman_do_ma_check = 1;

static double kalman_arma_ll (const double *b, void *p)
{
    int offset = kainfo->ifc + kainfo->np + kainfo->P;
    const double *theta = b + offset;
    const double *Theta = theta + kainfo->nq;
    double ll = NADBL;
    kalman *K;
    int err = 0;

#if ARMA_DEBUG
    debug_print_theta(theta, Theta);
#endif

    if (kalman_do_ma_check && ma_out_of_bounds(kainfo, theta, Theta)) {
	pputs(vprn, "arma: MA estimate(s) out of bounds\n");
	return NADBL;
    }

    K = (kalman *) p;
    err = rewrite_kalman_matrices(K, b, KALMAN_ALL);
    if (!err) {
	err = kalman_forecast(K);
	ll = kalman_get_loglik(K);
    }

#if ARMA_DEBUG
    fprintf(stderr, "kalman_arma_ll: loglik = %#.12g\n", ll);
#endif

    return ll;
}

static gretl_matrix *form_arma_y_vector (const int *alist, 
					 const double **Z,
					 struct arma_info *ainfo,
					 int *err)
{
    gretl_matrix *yvec;
    const double *y;
    int s, t, T;

#if ARMA_DEBUG
    fprintf(stderr, "ainfo->t1 = %d, ainfo->t2 = %d\n",
	    ainfo->t1, ainfo->t2);
#endif

    if (ainfo->dy != NULL) {
	y = ainfo->dy;
	/* should be handled earlier? */
	for (t=ainfo->t1; t<=ainfo->t2; t++) {
	    if (na(y[t])) {
		ainfo->t1 += 1;
	    } else {
		break;
	    }
	}
    } else {
	y = Z[ainfo->yno];
    }

    T = ainfo->t2 - ainfo->t1 + 1;
    if (T == 0) {
	*err = E_DATA;
	return NULL;
    }

    yvec = gretl_column_vector_alloc(T);
    if (yvec == NULL) {
	*err = E_ALLOC;
	return NULL;
    }    

    s = 0;
    for (t=ainfo->t1; t<=ainfo->t2; t++) {
	if (na(y[t])) {
	    *err = E_DATA;
	}
	gretl_vector_set(yvec, s++, y[t]);
    }

#if ARMA_DEBUG
    gretl_matrix_print(yvec, "y");
    fprintf(stderr, "y has %d rows\n", gretl_matrix_rows(yvec));
#endif

    if (*err) {
	gretl_matrix_free(yvec);
	yvec = NULL;
    }
    
    return yvec;
}

static gretl_matrix *form_arma_x_matrix (const int *alist, 
					 const double **Z,
					 struct arma_info *ainfo)
{
    gretl_matrix *x;
    int i, xstart;
    int *xlist;
    int err = 0;

    xlist = gretl_list_new(ainfo->nexo);
    if (xlist == NULL) {
	return NULL;
    }

    xstart = arma_list_y_position(ainfo) + 1;
    for (i=xstart; i<=alist[0]; i++) {
	xlist[i - xstart + 1] = alist[i];
    }

#if ARMA_DEBUG
    printlist(alist, "alist (arma list)");
    printlist(xlist, "xlist (exog vars)");
#endif

    x = gretl_matrix_data_subset(xlist, Z, ainfo->t1, ainfo->t2, 
				 NULL, &err);
    if (err) {
	free(xlist);
	return NULL;
    }

#if ARMA_DEBUG
    gretl_matrix_print(x, "x");
    fprintf(stderr, "x has %d rows\n", gretl_matrix_rows(x));
#endif

    free(xlist);
    
    return x;
}

/* Given an estimate of the ARMA constant via OLS, convert to the form
   wanted for initializing the Kalman filter.  Note: the 'b' array
   goes: const, phi, Phi, theta, Theta, beta.
*/

static void transform_arma_const (double *b, struct arma_info *ainfo)
{
    const double *phi = b + 1;
    const double *Phi = phi + ainfo->np;
    double narfac = 1.0;
    double sarfac = 1.0;
    int i, k = 0;

#if AINIT_DEBUG
    fprintf(stderr, "transform_arma_const: initially = %g\n", b[0]);
#endif

    for (i=0; i<ainfo->p; i++) {
	if (AR_included(ainfo, i)) {
	    narfac -= phi[k++];
	}
    }

    for (i=0; i<ainfo->P; i++) {
	sarfac -= Phi[i];
    }

    b[0] /= (narfac * sarfac);
}

static int kalman_undo_y_scaling (struct arma_info *ainfo,
				  gretl_matrix *y, double *b, 
				  kalman *K)
{
    double *beta = b + 1 + ainfo->np + ainfo->P +
	ainfo->nq + ainfo->Q;
    int i, t, T = ainfo->t2 - ainfo->t1 + 1;
    int err = 0;

    b[0] *= ainfo->dyscale;

    for (i=0; i<ainfo->nexo; i++) {
	beta[i] *= ainfo->dyscale;
    }

    i = ainfo->t1;
    for (t=0; t<T; t++) {
	y->val[t] *= ainfo->dyscale;
	ainfo->dy[i++] *= ainfo->dyscale;
    }

    if (na(kalman_arma_ll(b, K))) {
	err = 1;
    }

    return err;
}

static int arma_use_hessian (gretlopt opt)
{
    if (getenv("ARMA_USE_HESSIAN") != NULL) {
	return 1;
    } else if (opt & OPT_H) {
	return 1;
    } else if (libset_get_int(ARMA_VCV) == VCV_HESSIAN) {
	return 1;
    }

    return 0;
}

static int kalman_arma (const int *alist, double *coeff, 
			const double **Z, const DATAINFO *pdinfo,
			struct arma_info *ainfo, MODEL *pmod,
			gretlopt opt, PRN *prn)
{
    kalman *K = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *x = NULL;

    int k = 1 + ainfo->nexo; /* number of exog vars plus space for const */
    int r, r2, m;

    /* BFGS apparatus */
    int maxit = 1000;
    double reltol = 1.0e-12;
    int fncount = 0;
    int grcount = 0;

    double *b;
    int i, T, err = 0;

    b = malloc(ainfo->nc * sizeof *b);
    if (b == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<ainfo->nc; i++) {
	b[i] = coeff[i];
    }

#if AINIT_DEBUG
    fputs("initial coefficients:\n", stderr);
    for (i=0; i<ainfo->nc; i++) {
	fprintf(stderr, " b[%d] = % .10E\n", i, b[i]);
    }
#endif

    y = form_arma_y_vector(alist, Z, ainfo, &err);
    if (y == NULL) {
	free(b);
	return err;
    }

    if (ainfo->nexo > 0) {
	x = form_arma_x_matrix(alist, Z, ainfo);
	if (x == NULL) {
	    free(b);
	    gretl_matrix_free(y);
	    return E_ALLOC;
	}
    }

    if (allocate_ac_mc(ainfo)) {
	free(b);
	gretl_matrix_free(y);
	gretl_matrix_free(x);
	return E_ALLOC;
    }

    r = ainfo_get_r(ainfo);
    r2 = r * r;
    m = r * (r + 1) / 2;
    T = gretl_matrix_rows(y);

    /* when should we use vech apparatus? */
    if (r > 4) {
	set_arma_use_vech(ainfo);
    }

    clear_gretl_matrix_err();

    S = gretl_column_vector_alloc(r);
    P = gretl_matrix_alloc(r, r);
    F = gretl_matrix_alloc(r, r);
    A = gretl_column_vector_alloc(k);
    H = gretl_column_vector_alloc(r);
    E = gretl_column_vector_alloc(T);
    Q = gretl_matrix_alloc(r, r);

    Svar = gretl_matrix_alloc(r2, r2);

    if (arma_using_vech(ainfo)) {
	vQ = gretl_column_vector_alloc(m);
	Svar2 = gretl_matrix_alloc(m, m);
    } else {
	vQ = gretl_column_vector_alloc(r * r);
    }

    err = get_gretl_matrix_err();
    if (err) {
	free(b);
	gretl_matrix_free(y);
	gretl_matrix_free(x);
	return E_ALLOC;
    }

    kalman_matrices_init(ainfo);

#if ARMA_DEBUG
    fprintf(stderr, "ready to estimate: ainfo specs:\n"
	    "p=%d, P=%d, q=%d, Q=%d, ifc=%d, nexo=%d, t1=%d, t2=%d\n", 
	    ainfo->p, ainfo->P, ainfo->q, ainfo->Q, ainfo->ifc, 
	    ainfo->nexo, ainfo->t1, ainfo->t2);
    fprintf(stderr, "Kalman dims: r = %d, k = %d, T = %d, ncoeff=%d\n", 
	    r, k, T, ainfo->nc);
#endif

    /* publish ainfo */
    kainfo = ainfo;

    K = kalman_new(S, P, F, A, H, Q, NULL, y, x, E, ainfo->nc, 
		   ainfo->ifc, &err);

    if (err) {
	fprintf(stderr, "kalman_new(): err = %d\n", err);
    } else {
	if (r > 3) {
	    kalman_set_nonshift(K, 1);
	} else {
	    kalman_set_nonshift(K, r);
	}

	if (T > 3072 || getenv("KALMAN_AVG_LL") != NULL) {
	    kalman_set_options(K, KALMAN_ARMA_LL | KALMAN_AVG_LL);
	} else {
	    kalman_set_options(K, KALMAN_ARMA_LL);
	}

	err = BFGS_max(b, ainfo->nc, maxit, reltol, 
		       &fncount, &grcount, kalman_arma_ll, C_LOGLIK,
		       NULL, K, opt, prn);
	if (err) {
	    fprintf(stderr, "BFGS_max returned %d\n", err);
	} 
    }

    if (!err && ainfo->dyscale != 1.0) {
	kalman_undo_y_scaling(ainfo, y, b, K);
    }

    if (err) {
	pmod->errcode = err;
    } else {
	double *hess = NULL;

	gretl_model_set_int(pmod, "fncount", fncount);
	gretl_model_set_int(pmod, "grcount", grcount);

	if (arma_use_hessian(opt)) { 
	    kalman_do_ma_check = 0;
	    hess = numerical_hessian(b, ainfo->nc, kalman_arma_ll, K, &err);
	    kalman_do_ma_check = 1;
	    if (err) {
		/* fall back to OPG */
		err = 0;
	    }
	}

	kalman_arma_finish(pmod, alist, ainfo, Z, pdinfo, 
			   K, b, hess, ainfo->nc, T);
    } 

    kalman_free(K);

    gretl_matrix_free(S);
    gretl_matrix_free(P);
    gretl_matrix_free(F);
    gretl_matrix_free(A);
    gretl_matrix_free(H);
    gretl_matrix_free(E);
    gretl_matrix_free(Q);

    gretl_matrix_free(y);
    gretl_matrix_free(x);

    gretl_matrix_free(Svar);
    gretl_matrix_free(Svar2);
    gretl_matrix_free(vQ);

    free(b);
    free_ac_mc();

    /* unpublish ainfo */
    kainfo = NULL;

    return err;
}

/* end of Kalman-specific material */

/* construct a "virtual dataset" in the form of a set of pointers into
   the main dataset: this will be passed to the bhhh_max function.
   The dependent variable is put in position 0; following this are the
   independent variables.
*/

static const double **
make_armax_X (const int *list, struct arma_info *ainfo, const double **Z)
{
    const double **X;
    int ypos, nx;
    int v, i;

    ypos = arma_list_y_position(ainfo);
    nx = list[0] - ypos;

#if ARMA_DEBUG
    fprintf(stderr, "make_armax_X: allocating %d series pointers\n",
	    nx + 1);
#endif    

    X = malloc((nx + 1) * sizeof *X);
    if (X == NULL) {
	return NULL;
    }

    /* the dependent variable */
    if (ainfo->dy != NULL) {
	X[0] = ainfo->dy;
    } else {
	X[0] = Z[list[ypos]];
    }

    /* the independent variables */
    for (i=1; i<=nx; i++) {
	v = list[i + ypos];
	X[i] = Z[v];
    }

    return X;
}

static int add_to_spec (char *targ, const char *src)
{
    if (strlen(src) + strlen(targ) > MAXLINE - 1) {
	return 1;
    } else {
	strcat(targ, src);
	return 0;
    }
}

/* for ARMAX: write the component of the NLS specification
   that takes the form (y_{t-i} - X_{t-i} \beta)
*/

static int y_Xb_at_lag (char *spec, struct arma_info *ainfo, 
			int narmax, int lag)
{
    char chunk[32];
    int i, nt;
    int err = 0;

    if (narmax == 0) {
	sprintf(chunk, "y_%d", lag);
	return add_to_spec(spec, chunk);
    }

    nt = ainfo->ifc + narmax;

    sprintf(chunk, "(y_%d-", lag);

    if (nt > 1) {
	strcat(chunk, "(");
    }

    if (ainfo->ifc) {
	strcat(chunk, "b0");
    }

    err = add_to_spec(spec, chunk);

    for (i=0; i<narmax && !err; i++) {
	if (ainfo->ifc || i > 0) {
	    err += add_to_spec(spec, "+");
	} 
	sprintf(chunk, "b%d*x%d_%d", i+1, i+1, lag);
	err += add_to_spec(spec, chunk); 
    }

    if (nt > 1) {
	err += add_to_spec(spec, "))");
    } else {
	err += add_to_spec(spec, ")");
    }

    return err;
}

static void nls_kickstart (int b0, int by1, 
			   MODEL *pmod, double ***pZ, 
			   DATAINFO *pdinfo)
{
    int list[4];

    if (b0 != 0) {
	list[0] = 3;
	list[1] = 1;
	list[2] = 0;
	list[3] = 2;
    } else {
	list[0] = 2;
	list[1] = 1;
	list[2] = 2;
    }

    *pmod = lsq(list, pZ, pdinfo, OLS, OPT_A | OPT_Z);

    if (!pmod->errcode) {
	if (b0 != 0) {
	    (*pZ)[b0][0] = pmod->coeff[0];
	    (*pZ)[by1][0] = pmod->coeff[1];
	} else {
	    (*pZ)[by1][0] = pmod->coeff[0];
	}
    }

    clear_model(pmod);
}

static int arma_get_nls_model (MODEL *amod, struct arma_info *ainfo,
			       int narmax, double ***pZ, DATAINFO *pdinfo) 
{
#if AINIT_DEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR);
    gretlopt nlsopt = OPT_A | OPT_V;
#else
    PRN *prn = NULL;
    gretlopt nlsopt = OPT_A | OPT_C;
#endif
    char fnstr[MAXLINE];
    char term[32];
    nlspec *spec;
    int *plist = NULL;
    int v, oldv = pdinfo->v;
    int b0 = 0, by1 = 0;
    int nparam, lag;
    int i, j, k, err = 0;

    spec = nlspec_new(NLS, pdinfo);
    if (spec == NULL) {
	return E_ALLOC;
    }

    nlspec_set_t1_t2(spec, 0, ainfo->t2 - ainfo->t1); /* ?? */

    nparam = ainfo->ifc + ainfo->np + ainfo->P + ainfo->nexo;

    plist = gretl_list_new(nparam);
    if (plist == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = dataset_add_scalars(nparam, pZ, pdinfo); 
    if (err) {
	goto bailout;
    }

    /* make names for the parameters; construct the param list;
       and do some rudimentary fall-back initialization */

    v = oldv;
    k = 1;

    if (ainfo->ifc) {
	(*pZ)[v][0] = gretl_mean(0, pdinfo->n - 1, (*pZ)[1]);
	strcpy(pdinfo->varname[v], "b0");
	b0 = v;
	plist[k++] = v++;
    }

    for (i=0; i<ainfo->p; i++) {
	if (AR_included(ainfo, i)) {
	    if (by1 == 0) {
		by1 = v;
		(*pZ)[v][0] = 0.1;
	    }
	    sprintf(pdinfo->varname[v], "phi%d", i+1);
	    plist[k++] = v++;
	}
    }

    for (i=0; i<ainfo->P; i++) {
	if (by1 == 0) {
	    by1 = v;
	    (*pZ)[v][0] = 0.1;
	}
	sprintf(pdinfo->varname[v], "Phi%d", i+1);
	plist[k++] = v++;
    }

    for (i=0; i<ainfo->nexo; i++) {
	sprintf(pdinfo->varname[v], "b%d", i+1);
	plist[k++] = v++;
    }

    /* construct NLS specification */

    strcpy(fnstr, "y=");

    if (ainfo->ifc) {
	strcat(fnstr, "b0");
    } else {
	strcat(fnstr, "0");
    } 

    for (i=0; i<ainfo->p && !err; i++) {
	if (AR_included(ainfo, i)) {
	    lag = i + 1;
	    sprintf(term, "+phi%d*", lag);
	    err = add_to_spec(fnstr, term);
	    if (!err) {
		err = y_Xb_at_lag(fnstr, ainfo, narmax, lag);
	    }
	}
    }

    for (j=0; j<ainfo->P && !err; j++) {
	sprintf(term, "+Phi%d*", j+1);
	strcat(fnstr, term);
	lag = (j + 1) * ainfo->pd;
	y_Xb_at_lag(fnstr, ainfo, narmax, lag);
	for (i=0; i<ainfo->p; i++) {
	    if (AR_included(ainfo, i)) {
		sprintf(term, "-phi%d*Phi%d*", i+1, j+1);
		err = add_to_spec(fnstr, term);
		if (!err) {
		    lag = (j+1) * ainfo->pd + (i+1);
		    y_Xb_at_lag(fnstr, ainfo, narmax, lag);
		}
	    }
	}
    }

    for (i=0; i<ainfo->nexo && !err; i++) {
	sprintf(term, "+b%d*x%d", i+1, i+1);
	err = add_to_spec(fnstr, term);
    }

    if (!err) {
	nls_kickstart(b0, by1, amod, pZ, pdinfo);

#if AINIT_DEBUG
	fprintf(stderr, "initting using NLS spec:\n %s\n", fnstr);
	for (i=0; i<plist[0]; i++) {
	    fprintf(stderr, "initial NLS b[%d] = %g\n",
		    i, (*pZ)[oldv+i][0]);
	}
	printlist(plist, "NLS param list");
#endif

	err = nlspec_set_regression_function(spec, fnstr, pdinfo);
    }

    if (!err) {
	err = nlspec_add_param_list(spec, plist, (const double **) *pZ,
				    pdinfo);
    }

    if (!err) {
	*amod = model_from_nlspec(spec, pZ, pdinfo, nlsopt, prn);
	err = amod->errcode;
#if AINIT_DEBUG
	if (!err) {
	    printmodel(amod, pdinfo, OPT_NONE, prn);
	}
	gretl_print_destroy(prn);
#endif
    }

 bailout:

    nlspec_destroy(spec);
    free(plist);

    return err;
}

/* compose the regression list for the case where we're initializing
   ARMA via OLS (not NLS)
*/

static int *make_ar_ols_list (struct arma_info *ainfo, int av)
{
    int * alist = gretl_list_new(av);
    int i, k, vi;

    if (alist == NULL) {
	return NULL;
    }

    alist[1] = 1;

    if (ainfo->ifc) {
	alist[2] = 0;
	k = 3;
    } else {
	alist[0] -= 1;
	k = 2;
    }

    /* allow for const and y */
    vi = 2;

    for (i=0; i<ainfo->p; i++) {
	if (AR_included(ainfo, i)) {
	    alist[k++] = vi++;
	}
    }

    for (i=0; i<ainfo->P; i++) {
	alist[k++] = vi++;
    }

    for (i=0; i<ainfo->nexo; i++) {
	alist[k++] = vi++;
    }

    return alist;
}

/* compose variable names for temporary dataset */

static void arma_init_add_varnames (struct arma_info *ainfo, 
				    int ptotal, int narmax, 
				    DATAINFO *adinfo)
{
    int i, j, k, kx, ky;
    int lag, k0 = 2;

    strcpy(adinfo->varname[1], "y");

    k = k0;
    kx = ptotal + ainfo->nexo + k0;

    for (i=0; i<ainfo->p; i++) {
	if (AR_included(ainfo, i)) {
	    lag = i + 1;
	    sprintf(adinfo->varname[k++], "y_%d", lag);
	    for (j=0; j<narmax; j++) {
		sprintf(adinfo->varname[kx++], "x%d_%d", j+1, lag);
	    }
	}
    }

    ky = ainfo->np + ainfo->P + k0;

    for (j=0; j<ainfo->P; j++) {
	lag = (j + 1) * ainfo->pd;
	k = k0 + ainfo->np + j;
	sprintf(adinfo->varname[k], "y_%d", lag);
	for (i=0; i<narmax; i++) {
	    sprintf(adinfo->varname[kx++], "x%d_%d", i+1, lag);
	}
	for (i=0; i<ainfo->p; i++) {
	    if (AR_included(ainfo, i)) {
		lag = (j + 1) * ainfo->pd + (i + 1);
		sprintf(adinfo->varname[ky++], "y_%d", lag);
		for (k=0; k<narmax; k++) {
		    sprintf(adinfo->varname[kx++], "x%d_%d", k+1, lag);
		}
	    }
	}
    }

    kx = ptotal + k0;

    for (i=0; i<ainfo->nexo; i++) {
	sprintf(adinfo->varname[kx++], "x%d", i+1);
    }
}

/* Build temporary dataset including lagged vars: if we're doing exact
   ML on an ARMAX model we need lags of the exogenous variables as
   well as lags of y_t.  Note that the auxiliary dataset has "t = 0"
   at an offset of ainfo->t1 into the "real", external dataset.
*/

static void arma_init_build_dataset (struct arma_info *ainfo, 
				     int ptotal, int narmax, 
				     const int *list,
				     const double **Z,
				     double **aZ, 
				     DATAINFO *adinfo)
{
    const double *y;
    int i, j, k, kx, ky;
    int t, s, m, k0 = 2;
    int lag, xstart;

    /* dependent variable */
    if (ainfo->dy != NULL) {
	y = ainfo->dy;
    } else {
	y = Z[ainfo->yno];
    }

    /* starting position for reading exogeneous vars */
    if (ainfo->d > 0 || ainfo->D > 0) {
	xstart = (arma_has_seasonal(ainfo))? 10 : 6;
    } else {
	xstart = (arma_has_seasonal(ainfo))? 8 : 5;
    }

    for (t=0; t<adinfo->n; t++) {
	int realt = t + ainfo->t1;
	int miss = 0;

	aZ[1][t] = y[realt];

	k = k0;
	kx = ptotal + ainfo->nexo + k0;

	for (i=0; i<ainfo->p; i++) {
	    if (!AR_included(ainfo, i)) {
		continue;
	    }
	    lag = i + 1;
	    s = realt - lag;
	    if (s < 0) {
		miss = 1;
		aZ[k++][t] = NADBL;
		for (j=0; j<narmax; j++) {
		    aZ[kx++][t] = NADBL;
		}
	    } else {
		aZ[k++][t] = y[s];
		for (j=0; j<narmax; j++) {
		    m = list[xstart + j];
		    aZ[kx++][t] = Z[m][s];
		}
	    }
	}

	ky = ainfo->np + ainfo->P + k0;

	for (j=0; j<ainfo->P; j++) {
	    lag = (j + 1) * ainfo->pd;
	    s = realt - lag;
	    k = ainfo->np + k0 + j;
	    if (s < 0) {
		miss = 1;
		aZ[k][t] = NADBL;
		for (k=0; k<narmax; k++) {
		    aZ[kx++][t] = NADBL;
		}
	    } else {
		aZ[k][t] = y[s];
		for (k=0; k<narmax; k++) {
		    m = list[xstart + k];
		    aZ[kx++][t] = Z[m][s];
		}
	    }
	    for (i=0; i<ainfo->p; i++) {
		if (!AR_included(ainfo, i)) {
		    continue;
		}
		lag = (j + 1) * ainfo->pd + (i + 1);
		s = realt - lag;
		if (s < 0) {
		    miss = 1;
		    aZ[ky++][t] = NADBL;
		    for (k=0; k<narmax; k++) {
			aZ[kx++][t] = NADBL;
		    }
		} else {
		    aZ[ky++][t] = y[s];
		    for (k=0; k<narmax; k++) {
			m = list[xstart + k];
			aZ[kx++][t] = Z[m][s];
		    }
		}
	    }
	}

	kx = ptotal + k0;

	for (i=0; i<ainfo->nexo; i++) {
	    m = list[xstart + i];
	    aZ[kx++][t] = Z[m][realt];
	}

	if (miss) {
	    adinfo->t1 = t + 1;
	}	
    }

#if AINIT_DEBUG
    fprintf(stderr, "arma init dataset:\n");
    for (i=0; i<adinfo->v; i++) {
	fprintf(stderr, "var %d '%s', obs[0] = %g\n", i, adinfo->varname[i], 
		aZ[i][0]);
    }
#endif
}

/* transcribe coeffs from the OLS or NLS model used for initializing,
   into the array that will be passed to the maximizer.
*/

static void arma_init_transcribe_coeffs (struct arma_info *ainfo,
					 MODEL *pmod, double *b)
{
    int q0 = ainfo->ifc + ainfo->np + ainfo->P;
    int Q0 = q0 + ainfo->nq;
    int i, j = 0;

    for (i=0; i<pmod->ncoeff; i++) {
	if (i == q0) {
	    /* reserve space for nonseasonal MA */
	    j += ainfo->nq;
	} 
	if (i == Q0) {
	    /* and for seasonal MA */
	    j += ainfo->Q;
	}
	b[j++] = pmod->coeff[i];
    }

    /* insert near-zeros for nonseasonal MA */
    for (i=0; i<ainfo->nq; i++) {
	b[q0 + i] = 0.0001;
    } 

    /* and also seasonal MA */
    for (i=0; i<ainfo->Q; i++) {
	b[Q0 + i] = 0.0001;
    }	
}

static void maybe_rescale_dy (struct arma_info *ainfo)
{
    double ybar = gretl_mean(0, ainfo->T - 1, ainfo->dy);

    if (fabs(ybar) > 250) {
	int t;

	for (t=0; t<ainfo->T; t++) {
	    if (!na(ainfo->dy[t])) {
		ainfo->dy[t] /= ybar;
	    }
	}
	ainfo->dyscale = ybar;
    }
}

/* Run a least squares model to get initial values for the AR
   coefficients, either OLS or NLS.  We use NLS if there is
   nonlinearity due to either (a) the presence of both a seasonal and
   a non-seasonal AR component or (b) the presence of exogenous
   variables in the context of a non-zero AR order, where estimation
   will be via exact ML.  

   In this initialization any MA coefficients are simply set to
   near-zero.
*/

static int ar_arma_init (const int *list, double *coeff, 
			 const double **Z, const DATAINFO *pdinfo,
			 struct arma_info *ainfo, PRN *prn)
{
    int an = pdinfo->t2 - ainfo->t1 + 1;
    int nmixed = ainfo->np * ainfo->P;
    int ptotal = ainfo->np + ainfo->P + nmixed;
    int av = ptotal + ainfo->nexo + 2;
    double **aZ = NULL;
    DATAINFO *adinfo = NULL;
    int *alist = NULL;
    MODEL armod;
    int narmax, nonlin = 0;
    int i, err = 0;

#if AINIT_DEBUG
    fprintf(stderr, "ar_arma_init: pdinfo->t1=%d, pdinfo->t2=%d (n=%d); "
	    "ainfo->t1=%d, ainfo->t2=%d\n",
	    pdinfo->t1, pdinfo->t2, pdinfo->n, ainfo->t1, ainfo->t2);
    fprintf(stderr, " nmixed = %d, ptotal = %d\n", nmixed, ptotal);
#endif

    if (ptotal == 0 && ainfo->nexo == 0 && !ainfo->ifc) {
	/* special case of pure MA model */
	for (i=0; i<ainfo->nq + ainfo->Q; i++) {
	    coeff[i] = 0.0001; 
	} 
	return 0;
    }

    gretl_model_init(&armod); 

    narmax = (arma_exact_ml(ainfo))? ainfo->nexo : 0;
    if (narmax > 0) {
	/* ARMAX-induced lags of exog vars */
	av += ainfo->nexo * ptotal;
    } 

    if (arma_exact_ml(ainfo) && ainfo->ifc && ainfo->dy != NULL) {
	maybe_rescale_dy(ainfo);
    }

    adinfo = create_auxiliary_dataset(&aZ, av, an);
    if (adinfo == NULL) {
	return E_ALLOC;
    }

    if (ptotal > 0 && (narmax > 0 || nmixed > 0)) {
	/* we'll have to use NLS */
	nonlin = 1;
    } else {
	/* OLS: need regression list */
	alist = make_ar_ols_list(ainfo, av);
    }

    /* add variable names to auxiliary dataset: this is required for 
       NLS, and useful for debugging when using OLS 
    */
    arma_init_add_varnames(ainfo, ptotal, narmax, adinfo);

    /* build temporary dataset */
    arma_init_build_dataset(ainfo, ptotal, narmax, list,
			    Z, aZ, adinfo);

    if (nonlin) {
#if AINIT_DEBUG
	fprintf(stderr, "arma:_init_by_ls: doing NLS\n");
#endif
	err = arma_get_nls_model(&armod, ainfo, narmax, &aZ, adinfo);
    } else {
#if AINIT_DEBUG
	printlist(alist, "'alist' in ar_arma_init (OLS)");
#endif
	armod = lsq(alist, &aZ, adinfo, OLS, OPT_A | OPT_Z);
	err = armod.errcode;
    }

#if AINIT_DEBUG
    if (!err) {
	fprintf(stderr, "LS init: ncoeff = %d, nobs = %d\n", 
		armod.ncoeff, armod.nobs);
	for (i=0; i<armod.ncoeff; i++) {
	    fprintf(stderr, " coeff[%d] = %g\n", i, armod.coeff[i]);
	}
    } else {
	fprintf(stderr, "LS init: armod.errcode = %d\n", err);
    }
#endif

    if (!err) {
	arma_init_transcribe_coeffs(ainfo, &armod, coeff);
    }

    /* handle the case where we need to translate from an
       estimate of the regression constant to the
       unconditional mean of y_t
    */
    if (!err && arma_exact_ml(ainfo) && ainfo->ifc && 
	(!nonlin || ainfo->nexo == 0)) {
	transform_arma_const(coeff, ainfo);
    }

    if (!err && prn != NULL) {
	if (nonlin) {
	    pputs(prn, "\narma initialization: using nonlinear AR model\n\n");
	} else {
	    pputs(prn, "\narma initialization: using linear AR model\n\n");
	}
    }

    /* clean up */
    free(alist);
    destroy_auxiliary_dataset(aZ, adinfo);
    clear_model(&armod);

    return err;
}

#define MINLAGS 16

static int hr_init_check (const DATAINFO *pdinfo, struct arma_info *ainfo)
{
    int nobs = pdinfo->t2 - ainfo->t1 + 1; /* ?? */
    int nlags = (ainfo->P + ainfo->Q) * pdinfo->pd;
    int ncoeff, df;
    int err = 0;

    if (nlags < MINLAGS) {
	nlags = MINLAGS;
    }

    ncoeff = nlags + ainfo->nexo + ainfo->ifc;
    nobs -= nlags;
    df = nobs - ncoeff;

    if (df < 1) {
	err = E_DF;
    }

#if AINIT_DEBUG
    fprintf(stderr, "hr_init_check: ncoeff=%d, nobs=%d, 'df'=%d\n", 
	    ncoeff, nobs, df);
#endif

    return err;
}

static int hr_transcribe_coeffs (struct arma_info *ainfo,
				 MODEL *pmod, double *b)
{
    const double *theta = NULL;
    const double *Theta = NULL;
    int j = ainfo->nexo + ainfo->ifc;
    int i, k = 0;
    int err = 0;

    if (ainfo->ifc) {
	b[0] = pmod->coeff[0];
	k = 1;
    } 

    for (i=0; i<ainfo->p; i++) {
	if (AR_included(ainfo, i)) {
	    b[k++] = pmod->coeff[j++];
	}
    }

    for (i=0; i<ainfo->P; i++) { 
	b[k++] = pmod->coeff[j];
	j += ainfo->np + 1; /* assumes ainfo->p < pd */
    }

    theta = pmod->coeff + j;

    for (i=0; i<ainfo->q; i++) {
	if (MA_included(ainfo, i)) {
	    b[k++] = pmod->coeff[j++];
	}
    }

    Theta = pmod->coeff + j;

    for (i=0; i<ainfo->Q; i++) {
	b[k++] = pmod->coeff[j];
	j += ainfo->nq + 1; /* assumes ainfo->q < pd */
    }

    j = ainfo->ifc;

    for (i=0; i<ainfo->nexo; i++) {
	b[k++] = pmod->coeff[j++];
    }

    /* check MA values? */
    if (ainfo->q > 0 || ainfo->Q > 0) {
	err = ma_out_of_bounds(ainfo, theta, Theta);
	bounds_checker_cleanup();
    }

    return err;
}

/* Hannan-Rissanen ARMA initialization via two OLS passes. In the
   first pass we run an OLS regression of y on the exogenous vars plus
   a certain (biggish) number of lags. In the second we estimate the
   ARMA model by OLS, substituting innovations and corresponding lags
   with the first-pass residuals.
*/

static int hr_arma_init (const int *list, double *coeff, 
			 const double **Z, const DATAINFO *pdinfo,
			 struct arma_info *ainfo, PRN *prn)
{
    int an = pdinfo->t2 - ainfo->t1 + 1;
    int np = ainfo->p, nq = ainfo->q;
    int nP = ainfo->P, nQ = ainfo->Q;
    int ptotal = np + nP + np * nP;
    int qtotal = nq + nQ + nq * nQ;
    int nexo = ainfo->nexo;
    int pass1lags, pass1v;

    const double *y;
    double **aZ = NULL;
    DATAINFO *adinfo = NULL;
    int *pass1list = NULL;
    int *pass2list = NULL;
    int *arlags = NULL;
    int *malags = NULL;
    MODEL armod;
    int xstart;
    int m, pos, s;
    int i, j, t;
    int err = 0;

    pass1lags = (ainfo->Q + ainfo->P) * pdinfo->pd;
    if (pass1lags < MINLAGS) {
	pass1lags = MINLAGS;
    }
    pass1v = pass1lags + nexo + 2;

    /* dependent variable */
    if (ainfo->dy != NULL) {
	y = ainfo->dy;
    } else {
	y = Z[ainfo->yno];
    }

    adinfo = create_auxiliary_dataset(&aZ, pass1v + qtotal, an);
    if (adinfo == NULL) {
	return E_ALLOC;
    }

#if AINIT_DEBUG
    fprintf(stderr, "hr_arma_init: dataset allocated: %d vars, %d obs\n", 
	    pass1v + qtotal, an);
#endif

    /* in case we bomb before estimating a model */
    gretl_model_init(&armod);

    /* Start building stuff for pass 1 */

    pass1list = gretl_list_new(pass1v);
    if (pass1list == NULL) {
	err = E_ALLOC;
	goto bailout;
    }
	
    pass1list[1] = 1;
    pass1list[2] = 0;
    for (i=2; i<pass1v; i++) {
	pass1list[i+1] = i;
    }

    /* variable names */

    strcpy(adinfo->varname[1], "y");
    for (i=0; i<nexo; i++) { 
	/* exogenous vars */
	sprintf(adinfo->varname[i+1], "x%d", i);
    }
    for (i=1; i<=pass1lags; i++) { 
	/* lags */
	sprintf(adinfo->varname[i+1+nexo], "y_%d", i);
    }

     /* Fill the dataset with the data for pass 1 */

    /* starting position for reading exogeneous vars */
    if (ainfo->d > 0 || ainfo->D > 0) {
	xstart = (arma_has_seasonal(ainfo))? 10 : 6;
    } else {
	xstart = (arma_has_seasonal(ainfo))? 8 : 5;
    }

    for (t=0; t<an; t++) {
	s = t + ainfo->t1;
	aZ[1][t] = y[s];
	for (i=0, pos=2; i<nexo; i++) {
	    m = list[xstart + i];
	    aZ[pos++][t] = Z[m][s];
	}
	for (i=1; i<=pass1lags; i++) {
	    s = t + ainfo->t1 - i;
	    aZ[pos++][t] = (s >= 0)? y[s] : NADBL;
	}
    }

    /* pass 1 proper */

    armod = lsq(pass1list, &aZ, adinfo, OLS, OPT_A);
    if (armod.errcode) {
	err = armod.errcode;
	goto bailout;
    } 

#if AINIT_DEBUG
    fprintf(stderr, "pass1 model: t1=%d, t2=%d, nobs=%d, ncoeff=%d, dfd = %d\n", 
	    armod.t1, armod.t2, armod.nobs, armod.ncoeff, armod.dfd);
#endif

    /* allocations for pass 2 */

    if (qtotal > 0) {
	malags = malloc(qtotal * sizeof *malags);
	if (malags == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0, pos=0; i<nq; i++) {
		malags[pos++] = i+1;
	    }
	    for (i=0; i<ainfo->Q; i++) {
		for (j=0; j<=nq; j++) {
		    malags[pos++] = (i+1) * pdinfo->pd + j;
		}
	    }
	}
    }

    if (ptotal > 0 && !err) {
	arlags = malloc(ptotal * sizeof *arlags);
	if (arlags == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0, pos=0; i<np; i++) {
		arlags[pos++] = i+1;
	    }
	    for (i=0; i<ainfo->P; i++) {
		for (j=0; j<=np; j++) {
		    arlags[pos++] = (i+1) * pdinfo->pd + j;
		}
	    }
	}
    }

    if (!err) {
	pass2list = gretl_list_new(2 + nexo + ptotal + qtotal);
	if (pass2list == NULL) {
	    err = E_ALLOC;
	}
    }

    /* handle error in pass2 allocations */
    if (err) {
	goto bailout;
    }

    /* stick lagged residuals into temp dataset */
    pos = pass1v;
    for (i=0; i<qtotal; i++) {
	sprintf(adinfo->varname[pos], "e_%d", malags[i]);
	for (t=0; t<an; t++) {
	    s = t - malags[i];
	    aZ[pos][t] = (s >= 0)? armod.uhat[s] : NADBL;
	}
	pos++;
    }

    /* compose pass 2 regression list */
    for (i=1, pos=1; i<=nexo+2; i++) {
	pass2list[pos++] = pass1list[i];
    }
    for (i=0; i<ptotal; i++) {
	/* FIXME? */
	if (AR_included(ainfo,i)) {
	    pass2list[pos++] = arlags[i] + nexo + 1;
	}
    }
    for (i=0; i<qtotal; i++) {
	/* FIXME? */
	if (MA_included(ainfo,i)) {
	    pass2list[pos++] = pass1v + i;
	}
    }
    
    /* now do pass2 */
    clear_model(&armod);
    armod = lsq(pass2list, &aZ, adinfo, OLS, OPT_A);

    if (armod.errcode) {
	err = armod.errcode;
    } else {
#if AINIT_DEBUG
	PRN *modprn = gretl_print_new(GRETL_PRINT_STDERR);

	printmodel(&armod, adinfo, OPT_S, modprn);
	gretl_print_destroy(modprn);
#endif
	err = hr_transcribe_coeffs(ainfo, &armod, coeff);

	if (!err && arma_exact_ml(ainfo) && 
	    ainfo->ifc && ainfo->nexo == 0) {
	    transform_arma_const(coeff, ainfo);
	}
    }

#if AINIT_DEBUG
    if (!err) {
	fprintf(stderr, "HR init:\n");
	for (i=0; i<ainfo->nc; i++) {
	    fprintf(stderr, "coeff[%d] = %g\n", i, coeff[i]);
	}
    }
#endif

 bailout:

    free(pass1list);
    free(pass2list);
    free(arlags);
    free(malags);
    destroy_auxiliary_dataset(aZ, adinfo);
    clear_model(&armod);

    if (!err && prn != NULL) {
	pputs(prn, "\narma initialization: using Hannan-Rissanen method\n\n");
    }

    return err;
}

static int user_arma_init (double *coeff, struct arma_info *ainfo, 
			   char flags, int *init_done, PRN *prn)
{
    int i, nc = n_init_vals();

    if (nc == 0) {
	return 0;
    } else if (nc < ainfo->nc) {
	pprintf(prn, "arma initialization: need %d coeffs but got %d\n",
		ainfo->nc, nc);
	return E_DATA;
    }

    if (flags & ARMA_EXACT) {
	/* initialization is handled within BFGS */
	for (i=0; i<ainfo->nc; i++) {
	    coeff[i] = 0.0;
	}	
    } else {
	const gretl_matrix *m = get_init_vals();

	pputs(prn, "\narma initialization: at user-specified values\n\n");
	for (i=0; i<ainfo->nc; i++) {
	    coeff[i] = gretl_vector_get(m, i);
	}
	free_init_vals();
    }

    *init_done = 1;

    return 0;
}

/* set up a model_info struct for passing to bhhh_max */

static model_info *
set_up_arma_model_info (struct arma_info *ainfo)
{
    double tol = libset_get_double(BHHH_TOLER);
    model_info *arma;

    if (na(tol)) {
	tol = 1.0e-6;
    }

    arma = model_info_new(ainfo->nc, ainfo->t1, ainfo->t2, ainfo->T, tol);

    if (arma == NULL) return NULL;

    model_info_set_opts(arma, PRESERVE_OPG_MODEL);
    model_info_set_n_series(arma, ainfo->nc + 1);

    /* add pointer to ARMA-specific details */
    model_info_set_extra_info(arma, ainfo);

    return arma;
}

/* retrieve results specific to bhhh procedure */

static void 
conditional_arma_model_prep (MODEL *pmod, model_info *minfo,
			     double *theta)
{
    double **series;
    int i, t;

    pmod->lnL = model_info_get_ll(minfo);

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = theta[i];
    }

    series = model_info_get_series(minfo);
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->uhat[t] = series[0][t];
    }

    pmod->sigma = NADBL; /* will be replaced */
}

static int bhhh_arma (const int *alist, double *coeff, 
		      const double **Z, const DATAINFO *pdinfo,
		      struct arma_info *ainfo, MODEL *pmod,
		      gretlopt opt, PRN *prn)
{
    model_info *minfo = NULL;
    const double **X = NULL;
    int err = 0;

    /* construct virtual dataset for dep var, real regressors */
    X = make_armax_X(alist, ainfo, Z);
    if (X == NULL) {
	pmod->errcode = E_ALLOC;
	return pmod->errcode;
    }

    /* create model_info struct to feed to bhhh_max() */
    minfo = set_up_arma_model_info(ainfo);
    if (minfo == NULL) {
	pmod->errcode = E_ALLOC;
	free(X);
	return pmod->errcode;
    }

    /* call BHHH conditional ML function (OPG regression) */
    err = bhhh_max(arma_ll, X, coeff, minfo, opt, prn);
    
    if (err) {
	fprintf(stderr, "arma: bhhh_max returned %d\n", err);
	pmod->errcode = E_NOCONV;
    } else {
	MODEL *omod = model_info_capture_OPG_model(minfo);
	double *theta = model_info_get_theta(minfo);

	conditional_arma_model_prep(omod, minfo, theta);
	write_arma_model_stats(omod, alist, ainfo, Z, pdinfo);
	arma_model_add_roots(omod, ainfo, theta);
	*pmod = *omod;
	free(omod);
    }

    free(X);
    model_info_free(minfo);

    return pmod->errcode;
}

#if 0 /* not yet */

static int bhhh_arma_init (const int *alist, double *coeff, 
			   const double **Z, const DATAINFO *pdinfo,
			   struct arma_info *ainfo)
{
    model_info *minfo = NULL;
    const double **X = NULL;
    int err = 0;

    /* construct virtual dataset for dep var, real regressors */
    X = make_armax_X(alist, ainfo, Z);
    if (X == NULL) {
	return E_ALLOC;
    }

    /* create model_info struct to feed to bhhh_max() */
    minfo = set_up_arma_model_info(ainfo);
    if (minfo == NULL) {
	free(X);
	return E_ALLOC;
    }

    /* call BHHH conditional ML function (OPG regression) */
    err = bhhh_max(arma_ll, X, coeff, minfo, OPT_NONE, NULL);
    
    if (err) {
	fprintf(stderr, "bhhh_arma_init: bhhh_max returned %d\n", err);
	err = E_NOCONV;
    } else {
	double *theta = model_info_get_theta(minfo);
	int i;
	
	for (i=0; i<ainfo->nc; i++) {
	    coeff[i] = theta[i];
	}
    }

    free(X);
    model_info_free(minfo);

    return err;
}

#endif /* not yet */

/* Should we try Hannan-Rissanen initialization of ARMA
   coefficients? */

static int prefer_hr_init (struct arma_info *ainfo)
{
    int ret = 0;

    if (ainfo->q > 1 || ainfo->Q > 0) {
	ret = 1;

	/* don't use for gappy arma (yet?) */
	if (ainfo->pqspec != NULL && *ainfo->pqspec != '\0') {
	    ret = 0;
	}

	/* unlikely to work well with small sample */
	if (ainfo->t2 - ainfo->t1 < 100) {
	    ret = 0;
	}

#if 1 
	/* not sure about this: HR catches the MA terms, but NLS
	   handles better the AR interactions
	*/
	if (ainfo->p > 0 && ainfo->P > 0) {
	    ret = 0;
	}
#endif

	/* overlapping orders screw things up */
	if ((ainfo->P > 0 && ainfo->p >= ainfo->pd) ||
	    (ainfo->Q > 0 && ainfo->q >= ainfo->pd)) {
	    ret = 0;
	}

	if (ret && arma_exact_ml(ainfo)) {
	    /* screen for cases where we'll use NLS */
	    if (ainfo->P > 0) {
		ret = 0;
	    } else if (ainfo->p + ainfo->P > 0 && ainfo->nexo > 0) {
		ret = 0;
	    }
	}
    }

#if AINIT_DEBUG
    fprintf(stderr, "prefer_hr_init? %s\n", ret? "yes" : "no");
#endif

    return ret;
}

MODEL arma_model (const int *list, const char *pqspec,
		  const double **Z, const DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn)
{
    double *coeff = NULL;
    int *alist = NULL;
    MODEL armod;
    struct arma_info ainfo;
    int init_done = 0;
    char flags = 0;
    int err = 0;

    if (!(opt & OPT_C)) {
	flags = ARMA_EXACT;
    }

    vprn = set_up_verbose_printer(opt, prn);

    arma_info_init(&ainfo, flags, pqspec, pdinfo);
    gretl_model_init(&armod); 
    gretl_model_smpl_init(&armod, pdinfo);

    alist = gretl_list_copy(list);
    if (alist == NULL) {
	err = E_ALLOC;
    } 

    if (!err) {
	err = arma_check_list(alist, opt, Z, pdinfo, &ainfo);
    }

    if (err) {
	armod.errcode = err;
	goto bailout;
    } 

    /* calculate maximum lag */
    calc_max_lag(&ainfo);

    /* adjust sample range if need be */
    if (arma_adjust_sample(pdinfo, Z, alist, &ainfo)) {
        armod.errcode = E_DATA;
	goto bailout;
    }

    /* allocate initial coefficient vector */
    coeff = malloc(ainfo.nc * sizeof *coeff);
    if (coeff == NULL) {
	armod.errcode = E_ALLOC;
	goto bailout;
    }

    /* create differenced series if needed */
    if (ainfo.d > 0 || ainfo.D > 0) {
	err = arima_difference(Z[ainfo.yno], &ainfo);
    }

    /* initialize the coefficients: there are 3 possible methods */

    /* first pass: see if the user specified some values */
    err = user_arma_init(coeff, &ainfo, flags, &init_done, vprn);
    if (err) {
	armod.errcode = err;
	goto bailout;
    }

    /* second pass: try Hannan-Rissanen if suitable */
    if (!init_done && prefer_hr_init(&ainfo)) {
	err = hr_init_check(pdinfo, &ainfo);
	if (!err) {
	    err = hr_arma_init(alist, coeff, Z, pdinfo, &ainfo, vprn);
#if AINIT_DEBUG
	    if (err) {
		fputs("*** hr_arma_init failed, will try ar_arma_init\n", stderr);
	    } else {
		fputs("*** hr_arma_init OK\n", stderr);
	    }
#endif
	}
	if (!err) {
	    init_done = 1;
	}
    }

    /* third pass: estimate pure AR model by OLS or NLS */
    if (!init_done) {
	err = ar_arma_init(alist, coeff, Z, pdinfo, &ainfo, vprn);
    }

    if (err) {
	armod.errcode = err;
	goto bailout;
    }

    if (flags & ARMA_EXACT) {
#if 0 /* not yet */
	bhhh_arma_init(alist, coeff, Z, pdinfo, &ainfo);
#endif
	kalman_arma(alist, coeff, Z, pdinfo, &ainfo, &armod, opt, vprn);
    } else {
	bhhh_arma(alist, coeff, Z, pdinfo, &ainfo, &armod, opt, vprn);
    }

 bailout:

    free(alist);
    free(coeff);
    arma_info_cleanup(&ainfo);

    /* cleanup in MA roots checker */
    bounds_checker_cleanup();

    if (vprn != prn) {
	close_down_verbose_printer(vprn);
    }
    vprn = NULL;

    return armod;
}
