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
#include "gretl_scalar.h"
#include "gretl_bfgs.h"
#include "arma_priv.h"

#include "../cephes/libprob.h"

#define ARMA_DEBUG 0

/* ln(sqrt(2*pi)) + 0.5 */
#define LN_SQRT_2_PI_P5 1.41893853320467274178

#define KALMAN_ALL 999

#include "arma_common.c"

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

static struct bchecker *bchecker_allocate (arma_info *ainfo)
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

int ma_out_of_bounds (arma_info *ainfo, const double *theta,
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
	    pprintf(ainfo->prn, "MA root %d = %g\n", i, rt);
	    err = 1;
	    break;
	}
    }

    return err;
}

void bounds_checker_cleanup (void)
{
    ma_out_of_bounds(NULL, NULL, NULL);
}

static void do_MA_partials (double *drv,
			    arma_info *ainfo,
			    const double *theta,
			    const double *Theta,
			    int t)
{
    int i, j, k, p, s;

    k = 0;
    for (i=0; i<ainfo->q; i++) {
	if (MA_included(ainfo, i)) {
	    p = i + 1;
	    if (t - p >= 0) {
		drv[0] -= theta[k] * drv[p];
	    }
	    k++;
	}
    }

    for (j=0; j<ainfo->Q; j++) {
	p = (j + 1) * ainfo->pd;
	if (t - p >= 0) {
	    drv[0] -= Theta[j] * drv[p];
	    k = 0;
	    for (i=0; i<ainfo->q; i++) {
		if (MA_included(ainfo, i)) {
		    s = p + i + 1;
		    if (t - s >= 0) {
			drv[0] -= Theta[j] * theta[k] * drv[s];
		    }
		    k++;
		}
	    }
	}
    }
}

/* for each of the arrays of derivatives, shuffle each
   value one place up */

static void push_derivs (arma_info *ainfo, double **de, int dlen)
{
    int i, j;

    for (i=0; i<ainfo->n_aux; i++) {
	for (j=dlen-1; j>0; j--) {
	    de[i][j] = de[i][j-1];
	}
	de[i][0] = 0.0;
    }
}

static void zero_derivs (arma_info *ainfo, double **de, int dlen)
{
    int i, j;

    for (i=0; i<ainfo->n_aux; i++) {
	for (j=0; j<dlen; j++) {
	    de[i][j] = 0.0;
	}
    }
}

static int arma_analytical_score (arma_info *ainfo,
				  const double *y,
				  const double **X,
				  const double *phi,
				  const double *Phi,
				  const double *theta,
				  const double *Theta,
				  double s2,
				  gretl_matrix *G)
{
    /* forecast errors */
    const double *e = ainfo->e;
    /* pointers to blocks of derivatives (workspace) */
    double **de = ainfo->aux;
    double **de_a =    de + ainfo->ifc;
    double **de_sa = de_a + ainfo->np;
    double **de_m = de_sa + ainfo->P;
    double **de_sm = de_m + ainfo->nq;
    double **de_r = de_sm + ainfo->Q; 
    int dlen = 1 + ainfo->q + ainfo->pd * ainfo->Q;
    double x, Gsi;
    int t, gt;
    int i, j, k, p, s;

    zero_derivs(ainfo, de, dlen);

    for (t=ainfo->t1, gt=0; t<=ainfo->t2; t++, gt++) {

	/* the constant term (de_0) */
	if (ainfo->ifc) {
	    de[0][0] = -1.0;
	    do_MA_partials(de[0], ainfo, theta, Theta, t);
	}

	/* non-seasonal AR terms (de_a) */
	k = 0;
	for (i=0; i<ainfo->p; i++) {
	    if (!AR_included(ainfo, i)) {
		continue;
	    }
	    p = i + 1;
	    if (t - p >= 0) {
		de_a[k][0] = -y[t-p];
		/* cross-partial with seasonal AR */
		for (j=0; j<ainfo->P; j++) {
		    s = p + (j + 1) * ainfo->pd;
		    if (t - s >= 0) {
			de_a[k][0] += Phi[j] * y[t-s];
		    }
		}
		do_MA_partials(de_a[k], ainfo, theta, Theta, t);
	    }
	    k++;
	}

	/* seasonal AR terms (de_sa) */
	for (j=0; j<ainfo->P; j++) {
	    p = (j + 1) * ainfo->pd;
	    if (t - p >= 0) {
		de_sa[j][0] = -y[t-p];
		/* cross-partial with non-seasonal AR */
		k = 0;
		for (i=0; i<ainfo->p; i++) {
		    if (AR_included(ainfo, i)) {
			s = p + i + 1;
			if (t - s >= 0) {
			    de_sa[j][0] += phi[k] * y[t-s];
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
	    p = i + 1;
	    if (t - p >= 0) {
		de_m[k][0] = -e[t-p];
		/* cross-partial with seasonal MA */
		for (j=0; j<ainfo->Q; j++) {
		    s = p + (j + 1) * ainfo->pd;
		    if (t - s >= 0) {
			de_m[k][0] -= Theta[j] * e[t-s];
		    }
		}
		do_MA_partials(de_m[k], ainfo, theta, Theta, t);
	    }
	    k++;
	}

	/* seasonal MA terms (de_sm) */
	for (j=0; j<ainfo->Q; j++) {
	    p = (j + 1) * ainfo->pd;
	    if (t - p >= 0) {
		de_sm[j][0] = -e[t-p];
		/* cross-partial with non-seasonal MA */
		k = 0;
		for (i=0; i<ainfo->q; i++) {
		    if (MA_included(ainfo, i)) {
			s = p + i + 1;
			if (t - s >= 0) {
			    de_sm[j][0] -= theta[k] * e[t-s];
			}
			k++;
		    }
		}
		do_MA_partials(de_sm[j], ainfo, theta, Theta, t);
	    }
	}

	/* exogenous regressors (de_r) */
	for (j=0; j<ainfo->nexo; j++) {
	    de_r[j][0] = -X[j][t]; 
	    do_MA_partials(de_r[j], ainfo, theta, Theta, t);
	}

	/* update gradient matrix */
	x = e[t] / s2; /* sqrt(s2)? does it matter? */
	for (i=0; i<ainfo->nc; i++) {
	    Gsi = -de[i][0] * x;
	    gretl_matrix_set(G, gt, i, Gsi);
	}

	push_derivs(ainfo, de, dlen);
    }

    return 0;
}

static int conditional_arma_forecast_errors (arma_info *ainfo,
					     const double *y,
					     const double **X,
					     double b0,
					     const double *phi,
					     const double *Phi,
					     const double *theta,
					     const double *Theta,
					     const double *beta,
					     double *s2)
{
    double *e = ainfo->e;
    int i, j, k, s, t, p;

    *s2 = 0.0;

    for (t=ainfo->t1; t<=ainfo->t2; t++) {
	e[t] = y[t];

	/* intercept */
	if (ainfo->ifc) {
	    e[t] -= b0;
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
		if (s >= ainfo->t1) {
		    e[t] -= theta[k] * e[s];
		}
		k++;
	    }
	}

	/* seasonal MA component plus interactions */
	for (j=0; j<ainfo->Q; j++) {
	    s = t - (j + 1) * ainfo->pd;
	    if (s >= ainfo->t1) {
		e[t] -= Theta[j] * e[s];
		k = 0;
		for (i=0; i<ainfo->q; i++) {
		    if (MA_included(ainfo, i)) {
			p = s - (i + 1);
			if (p >= ainfo->t1) {
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

	*s2 += e[t] * e[t];
    }

    return 0;
}

/* Calculate ARMA log-likelihood.  This function is passed to the
   bhhh_max() routine as a callback. */

static double bhhh_arma_callback (double *coeff, 
				  gretl_matrix *G, 
				  void *data,
				  int do_score,
				  int *err)
{
    arma_info *ainfo = (arma_info *) data;
    /* pointers to blocks of data */
    const double *y = ainfo->X[0];
    const double **X = ainfo->X + 1;
    /* pointers to blocks of coefficients */
    const double *phi =   coeff + ainfo->ifc;
    const double *Phi =     phi + ainfo->np;
    const double *theta =   Phi + ainfo->P;
    const double *Theta = theta + ainfo->nq;
    const double *beta =  Theta + ainfo->Q;
    double ll, s2 = 0.0;

    *err = 0;

#if ARMA_DEBUG
    fprintf(stderr, "arma_ll: p=%d, q=%d, P=%d, Q=%d, pd=%d\n",
	    ainfo->p, ainfo->q, ainfo->P, ainfo->Q, ainfo->pd);
#endif

    if (ma_out_of_bounds(ainfo, theta, Theta)) {
	pputs(ainfo->prn, "arma: MA estimate(s) out of bounds\n");
	fputs("arma: MA estimate(s) out of bounds\n", stderr);
	*err = E_NOCONV;
	return NADBL;
    }

    conditional_arma_forecast_errors(ainfo, y, X, coeff[0],
				     phi, Phi, theta, Theta,
				     beta, &s2);

    /* error variance and log-likelihood */
    s2 /= (double) ainfo->T;
    ll = -ainfo->T * (0.5 * log(s2) + LN_SQRT_2_PI_P5);

    if (isnan(ll)) {
	*err = E_NAN;
    }

    if (do_score) {
	ainfo->ll = ll;
	arma_analytical_score(ainfo, y, X, 
			      phi, Phi, theta, Theta,
			      s2, G);
    }

    return ll;
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

static int arma_model_add_roots (MODEL *pmod, arma_info *ainfo,
				 const double *coeff)
{
    const double *phi =   coeff + ainfo->ifc;
    const double *Phi =     phi + ainfo->np;
    const double *theta =   Phi + ainfo->P;
    const double *Theta = theta + ainfo->nq;
    int nr = ainfo->p + ainfo->P + ainfo->q + ainfo->Q;
    int pmax, qmax, lmax;
    double *temp = NULL, *tmp2 = NULL;
    cmplx *rptr, *roots = NULL;
    int i, k, cerr = 0;

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
	cerr = polrt(temp, tmp2, ainfo->p, rptr);
	rptr += ainfo->p;
    }

    if (!cerr && ainfo->P > 0) {
	/* B(L), seasonal */
	for (i=0; i<ainfo->P; i++) {
	    temp[i+1] = -Phi[i];
	}    
	cerr = polrt(temp, tmp2, ainfo->P, rptr);
	rptr += ainfo->P;
    }

    if (!cerr && ainfo->q > 0) {
	/* C(L), non-seasonal */
	k = 0;
	for (i=0; i<ainfo->q; i++) {
	    if (MA_included(ainfo, i)) {
		temp[i+1] = theta[k++];
	    } else {
		temp[i+1] = 0;
	    }
	}  
	cerr = polrt(temp, tmp2, ainfo->q, rptr);
	rptr += ainfo->q;
    }

    if (!cerr && ainfo->Q > 0) {
	/* D(L), seasonal */
	for (i=0; i<ainfo->Q; i++) {
	    temp[i+1] = Theta[i];
	}  
	cerr = polrt(temp, tmp2, ainfo->Q, rptr);
    }
    
    free(temp);
    free(tmp2);

    if (cerr) {
	free(roots);
    } else {
	gretl_model_set_data(pmod, "roots", roots, GRETL_TYPE_CMPLX_ARRAY,
			     nr * sizeof *roots);
    }

    return 0;
}

/* below: exact ML using Kalman filter apparatus */

static gretl_matrix *S;
static gretl_matrix *P;
static gretl_matrix *F;
static gretl_matrix *A;
static gretl_matrix *H;
static gretl_matrix *Q;
static gretl_matrix *E;
static gretl_matrix *Svar;

static gretl_matrix *Svar2;
static gretl_matrix *vQ;

static arma_info *kainfo;   

static void allocate_kalman_matrices (int r, int r2,
				      int k, int T)
{
    S = gretl_column_vector_alloc(r);
    P = gretl_matrix_alloc(r, r);
    F = gretl_matrix_alloc(r, r);
    A = gretl_column_vector_alloc(k);
    H = gretl_column_vector_alloc(r);
    E = gretl_column_vector_alloc(T);
    Q = gretl_matrix_alloc(r, r);
    Svar = gretl_matrix_alloc(r2, r2);
}

static void free_kalman_matrices (void)
{
    gretl_matrix_replace(&S, NULL);
    gretl_matrix_replace(&P, NULL);
    gretl_matrix_replace(&F, NULL);
    gretl_matrix_replace(&A, NULL);
    gretl_matrix_replace(&H, NULL);
    gretl_matrix_replace(&E, NULL);
    gretl_matrix_replace(&Q, NULL);
    gretl_matrix_replace(&Svar, NULL);
}

static int ainfo_get_r (arma_info *ainfo)
{
    int pmax = ainfo->p + ainfo->pd * ainfo->P;
    int qmax = ainfo->q + ainfo->pd * ainfo->Q;

    return (pmax > qmax + 1)? pmax : qmax + 1;
}

static int allocate_ac_mc (arma_info *ainfo)
{
    int m = (ainfo->P > 0) + (ainfo->Q > 0);
    int err = 0;

    if (m > 0) {
	double *ac = NULL, *mc = NULL;
	int n, i = 0;

	ainfo->aux = doubles_array_new(m, 0);
	if (ainfo->aux == NULL) {
	    return E_ALLOC;
	}

	if (ainfo->P > 0) {
	    n = 1 + ainfo->p + ainfo->pd * ainfo->P;
	    ac = malloc(n * sizeof *ac);
	    if (ac == NULL) {
		err = E_ALLOC;
	    } else {
		ainfo->aux[i++] = ac;
	    }
	}

	if (!err && ainfo->Q > 0) {
	    n = 1 + ainfo->q + ainfo->pd * ainfo->Q;
	    mc = malloc(n * sizeof *mc);
	    if (mc == NULL) {
		err = E_ALLOC;
	    } else {
		ainfo->aux[i++] = mc;
	    }
	}

	if (err) {
	    doubles_array_free(ainfo->aux, m);
	} else {
	    ainfo->n_aux = m;
	}
    }

    return err;
}

static void write_big_phi (const double *phi, 
			   const double *Phi,
			   arma_info *ainfo,
			   gretl_matrix *F)
{
    int pmax = ainfo->p + ainfo->pd * ainfo->P;
    double *ac = ainfo->aux[0];
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
			     arma_info *ainfo,
			     gretl_matrix *H)
{
    int qmax = ainfo->q + ainfo->pd * ainfo->Q;
    int i = (ainfo->P > 0)? 1 : 0;
    double *mc = ainfo->aux[i];
    double x, y;
    int j, k, ii;

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

static void kalman_matrices_init (arma_info *ainfo)
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
    const double *phi =       b + kainfo->ifc;
    const double *Phi =     phi + kainfo->np;
    const double *theta =   Phi + kainfo->P;
    const double *Theta = theta + kainfo->nq;
    const double *beta =  Theta + kainfo->Q;
    double mu = (kainfo->ifc)? b[0] : 0.0;
    int rewrite_A = 0;
    int rewrite_F = 0;
    int rewrite_H = 0;
    int i, k, err = 0;

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

static const double *kalman_arma_llt_callback (const double *b, int i, 
					       void *data)
{
    kalman *K = (kalman *) data;
    int err;

    rewrite_kalman_matrices(K, b, i);
    err = kalman_forecast(K, NULL);

#if ARMA_DEBUG
    fprintf(stderr, "kalman_arma_llt: kalman f'cast gave "
	    "err = %d, ll = %#.12g\n", err, kalman_get_loglik(K));
#endif

    return (err)? NULL : E->val;
}

/* add covariance matrix and standard errors based on Outer Product of
   Gradient
*/

static int arma_OPG_vcv (MODEL *pmod, kalman *K, double *b, 
			 double s2, int k, int T,
			 PRN *prn)
{
    gretl_matrix *G = NULL;
    gretl_matrix *V = NULL;
    double rcond;
    int s, t;
    int err = 0;

    for (t=pmod->t1, s=0; t<=pmod->t2; t++, s++) {
	pmod->uhat[t] = gretl_vector_get(E, s);
    }

    G = build_score_matrix(b, k, T, kalman_arma_llt_callback, 
			   K, &err);
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

    rcond = gretl_symmetric_matrix_rcond(V, &err);
    if (!err && rcond < 1.0E-10) {
	pprintf(prn, "OPG: rcond = %g; will try Hessian\n", rcond);
	err = 1;
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(V);
    }

    if (!err) {
	gretl_matrix_multiply_by_scalar(V, s2);
	err = gretl_model_write_vcv(pmod, V);
    }

 bailout:

    gretl_matrix_free(G);
    gretl_matrix_free(V);
    
    return err;
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

static double kalman_arma_ll (const double *b, void *data)
{
    kalman *K = (kalman *) data;
    int offset = kainfo->ifc + kainfo->np + kainfo->P;
    const double *theta = b + offset;
    const double *Theta = theta + kainfo->nq;
    double ll = NADBL;
    int err = 0;

#if ARMA_DEBUG
    debug_print_theta(theta, Theta);
#endif

    if (kalman_do_ma_check && ma_out_of_bounds(kainfo, theta, Theta)) {
	pputs(kalman_get_printer(K), "arma: MA estimate(s) out of bounds\n");
	return NADBL;
    }

    err = rewrite_kalman_matrices(K, b, KALMAN_ALL);
    if (!err) {
	err = kalman_forecast(K, NULL);
	ll = kalman_get_loglik(K);
    }

#if ARMA_DEBUG
    fprintf(stderr, "kalman_arma_ll: loglik = %#.12g\n", ll);
#endif

    return ll;
}

static void maybe_rescale_hessian (int kopt, double *vcv, 
				   int k, int T)
{
    if (kopt & KALMAN_AVG_LL) {
	int i, k2 = k * (k + 1) / 2;

	for (i=0; i<k2; i++) {
	    vcv[i] /= T;
	}
    }
}

static int arima_ydiff_only (arma_info *ainfo)
{
    if ((ainfo->d > 0 || ainfo->D > 0) &&
	ainfo->nexo > 0 && !(ainfo->flags & ARMA_XDIFF)) {
	return 1;
    } else {
	return 0;
    }
}

static int arma_use_hessian (gretlopt opt)
{
    int ret = 1;

    if (opt & OPT_G) {
	ret = 0;
    } else if (libset_get_int(ARMA_VCV) == VCV_OP) {
	ret = 0;
    }

    return ret;
}

static int kalman_arma_finish (MODEL *pmod, arma_info *ainfo,
			       const double **Z, const DATAINFO *pdinfo, 
			       kalman *K, double *b, 
			       gretlopt opt, PRN *prn)
{
    double s2;
    int do_hess = arma_use_hessian(opt);
    int kopt, i, k = ainfo->nc;
    int err;

    pmod->t1 = ainfo->t1;
    pmod->t2 = ainfo->t2;
    pmod->nobs = ainfo->T;
    pmod->ncoeff = ainfo->nc;
    pmod->full_n = pdinfo->n;
    
    /* in the Kalman case the basic model struct is empty, so we 
       have to allocate for coefficients, residuals and so on
    */

    err = gretl_model_allocate_storage(pmod);
    if (err) {
	return err;
    }

    for (i=0; i<k; i++) {
	pmod->coeff[i] = b[i];
    }

    s2 = kalman_get_arma_variance(K);
    pmod->sigma = sqrt(s2);

    pmod->lnL = kalman_get_loglik(K);
    kopt = kalman_get_options(K);

    /* rescale if we're using average loglikelihood */
    if (kopt & KALMAN_AVG_LL) {
	pmod->lnL *= ainfo->T;
    }

#if ARMA_DEBUG
    fprintf(stderr, "kalman_arma_finish: doing VCV, method %s\n",
	    (do_hess)? "Hessian" : "OPG");
#endif

    if (!do_hess) {
	/* try OPG first, hessian as fallback */
	err = arma_OPG_vcv(pmod, K, b, s2, k, ainfo->T, prn);
	if (err) {
	    err = 0;
	    do_hess = 1;
	} else {
	    gretl_model_set_vcv_info(pmod, VCV_ML, VCV_OP);
	    pmod->opt |= OPT_G;
	}
    }	

    if (do_hess) { 
	double *vcv;

	kalman_do_ma_check = 0;
	vcv = numerical_hessian(b, ainfo->nc, kalman_arma_ll, K, &err);
	kalman_do_ma_check = 1;
	if (!err) {
	    maybe_rescale_hessian(kopt, vcv, k, ainfo->T);
	    arma_hessian_vcv(pmod, vcv, k);
	    gretl_model_set_vcv_info(pmod, VCV_ML, VCV_HESSIAN);
	}
    }

    if (!err) {
	write_arma_model_stats(pmod, ainfo, Z, pdinfo);
	arma_model_add_roots(pmod, ainfo, b);
	gretl_model_set_int(pmod, "arma_flags", ARMA_EXACT);
	if (ainfo->flags & ARMA_LBFGS) {
	    pmod->opt |= OPT_L;
	}
	if (arima_ydiff_only(ainfo)) {
	    pmod->opt |= OPT_Y;
	}
    }

    return err;
}

static gretl_matrix *form_arma_y_vector (arma_info *ainfo,
					 const double **Z,
					 int *err)
{
    gretl_matrix *yvec;
    const double *y;
    int s, t;

#if ARMA_DEBUG
    fprintf(stderr, "ainfo->t1 = %d, ainfo->t2 = %d\n",
	    ainfo->t1, ainfo->t2);
#endif

    if (ainfo->y != NULL) {
	y = ainfo->y;
    } else {
	y = Z[ainfo->yno];
    }

    yvec = gretl_column_vector_alloc(ainfo->T);
    if (yvec == NULL) {
	*err = E_ALLOC;
	return NULL;
    }    

    s = 0;
    for (t=ainfo->t1; t<=ainfo->t2; t++) {
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

static gretl_matrix *form_arma_X_matrix (arma_info *ainfo,
					 const double **Z,
					 int *err)
{
    gretl_matrix *X;

#if ARMA_DEBUG
    printlist(ainfo->xlist, "xlist (exog vars)");
#endif

    X = gretl_matrix_data_subset(ainfo->xlist, Z, ainfo->t1, ainfo->t2, 
				 M_MISSING_ERROR, err);

#if ARMA_DEBUG
    gretl_matrix_print(X, "X");
    fprintf(stderr, "X has %d rows\n", gretl_matrix_rows(X));
#endif

    return X;
}

static int kalman_undo_y_scaling (arma_info *ainfo,
				  gretl_matrix *y, double *b, 
				  kalman *K)
{
    double *beta = b + 1 + ainfo->np + ainfo->P +
	ainfo->nq + ainfo->Q;
    int i, t, T = ainfo->t2 - ainfo->t1 + 1;
    int err = 0;

    b[0] *= ainfo->yscale;

    for (i=0; i<ainfo->nexo; i++) {
	beta[i] *= ainfo->yscale;
    }

    i = ainfo->t1;
    for (t=0; t<T; t++) {
	y->val[t] *= ainfo->yscale;
	ainfo->y[i++] *= ainfo->yscale;
    }

    if (na(kalman_arma_ll(b, K))) {
	err = 1;
    }

    return err;
}

static void free_arma_X_matrix (arma_info *ainfo, gretl_matrix *X)
{
    if (X == ainfo->dX) {
	gretl_matrix_free(ainfo->dX);
	ainfo->dX = NULL;
    } else {
	gretl_matrix_free(X);
    }
}

#define KALMAN_ARMA_INITH 0 /* try initializing inverse Hessian */
#define INITH_USE_OPG 0     /* try using OPG for that purpose */

#if KALMAN_ARMA_INITH
# if INITH_USE_OPG

static gretl_matrix *kalman_arma_init_H (double *b, int k, int T,
					 kalman *K)
{
    gretl_matrix *G, *H = NULL;
    int err = 0;

    G = build_score_matrix(b, k, T, kalman_arma_llt_callback, 
			   K, &err);

    if (!err) {
	H = gretl_matrix_alloc(k, k);
	if (H == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	gretl_matrix_multiply_mod(G, GRETL_MOD_NONE,
				  G, GRETL_MOD_TRANSPOSE,
				  H, GRETL_MOD_NONE);
	err = gretl_invert_symmetric_matrix(H);
	if (err) {
	    fprintf(stderr, "kalman_arma_init_H: H is not pd\n");
	    gretl_matrix_free(H);
	    H = NULL;
	}
    }

    gretl_matrix_free(G);

    return H;
}

# else /* just use scaled identity matrix */

static gretl_matrix *kalman_arma_init_H (double *b, int k, int T,
					 kalman *K)
{
    gretl_matrix *H;

    H = gretl_identity_matrix_new(k);
    if (H != NULL) {
	gretl_matrix_divide_by_scalar(H, T);
    } 

    return H;
}

# endif /* INITH variants */
#endif /* KALMAN_ARMA_INITH */

static int kalman_arma (double *coeff, 
			const double **Z, const DATAINFO *pdinfo,
			arma_info *ainfo, MODEL *pmod,
			gretlopt opt)
{
    kalman *K = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    int k = 1 + ainfo->nexo; /* number of exog vars plus space for const */
    int r, r2, m;
    int fncount = 0, grcount = 0;
    double *b;
    int i, err = 0;
    
    b = malloc(ainfo->nc * sizeof *b);
    if (b == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<ainfo->nc; i++) {
	b[i] = coeff[i];
    }

#if ARMA_DEBUG
    fputs("initial coefficients:\n", stderr);
    for (i=0; i<ainfo->nc; i++) {
	fprintf(stderr, " b[%d] = % .10E\n", i, b[i]);
    }
#endif

    y = form_arma_y_vector(ainfo, Z, &err);

    if (!err && ainfo->nexo > 0) {
	if (ainfo->dX != NULL) {
	    X = ainfo->dX;
	} else {
	    X = form_arma_X_matrix(ainfo, Z, &err);
	}
    }

    if (!err) {
	err = allocate_ac_mc(ainfo);
    }

    if (err) {
	goto bailout;
    }

    r = ainfo_get_r(ainfo);
    r2 = r * r;
    m = r * (r + 1) / 2;

    /* when should we use vech apparatus? */
    if (r > 4) {
	set_arma_use_vech(ainfo);
    }

    clear_gretl_matrix_err();

    allocate_kalman_matrices(r, r2, k, ainfo->T);

    if (arma_using_vech(ainfo)) {
	vQ = gretl_column_vector_alloc(m);
	Svar2 = gretl_matrix_alloc(m, m);
    } else {
	vQ = gretl_column_vector_alloc(r * r);
    }

    err = get_gretl_matrix_err();
    if (err) {
	goto bailout;
    }

    kalman_matrices_init(ainfo);

#if ARMA_DEBUG
    fprintf(stderr, "ready to estimate: ainfo specs:\n"
	    "p=%d, P=%d, q=%d, Q=%d, ifc=%d, nexo=%d, t1=%d, t2=%d\n", 
	    ainfo->p, ainfo->P, ainfo->q, ainfo->Q, ainfo->ifc, 
	    ainfo->nexo, ainfo->t1, ainfo->t2);
    fprintf(stderr, "Kalman dims: r = %d, k = %d, T = %d, ncoeff=%d\n", 
	    r, k, ainfo->T, ainfo->nc);
#endif

    /* publish ainfo */
    kainfo = ainfo;

    K = kalman_new(S, P, F, A, H, Q, NULL, y, X, E, &err);

    if (err) {
	fprintf(stderr, "kalman_new(): err = %d\n", err);
    } else {
	int maxit, save_lbfgs = 0;
	double toler;
	gretl_matrix *A0 = NULL;

	kalman_attach_printer(K, ainfo->prn);

	if (r > 3) {
	    kalman_set_nonshift(K, 1);
	} else {
	    kalman_set_nonshift(K, r);
	}

	if (getenv("KALMAN_AVG_LL") != NULL) {
	    kalman_set_options(K, KALMAN_ARMA_LL | KALMAN_AVG_LL);
	} else {
	    kalman_set_options(K, KALMAN_ARMA_LL);
	}

	save_lbfgs = libset_get_bool(USE_LBFGS);

	if (save_lbfgs) {
	    ainfo->flags |= ARMA_LBFGS;
	} else if (opt & OPT_L) {
	    libset_set_bool(USE_LBFGS, 1);
	    ainfo->flags |= ARMA_LBFGS;
	}

	BFGS_defaults(&maxit, &toler, ARMA);

#if KALMAN_ARMA_INITH
	A0 = kalman_arma_init_H(b, ainfo->nc, ainfo->T, K);
#endif

	err = BFGS_max(b, ainfo->nc, maxit, toler, 
		       &fncount, &grcount, kalman_arma_ll, C_LOGLIK,
		       NULL, K, A0, opt, ainfo->prn);

	gretl_matrix_free(A0);

	if (err) {
	    fprintf(stderr, "BFGS_max returned %d\n", err);
	} 

	if (!save_lbfgs && (opt & OPT_L)) {
	    libset_set_bool(USE_LBFGS, 0);
	}
    }

    if (!err && ainfo->yscale != 1.0) {
	kalman_undo_y_scaling(ainfo, y, b, K);
    }

    if (!err) {
	gretl_model_set_int(pmod, "fncount", fncount);
	gretl_model_set_int(pmod, "grcount", grcount);
	err = kalman_arma_finish(pmod, ainfo, 
				 Z, pdinfo, K, b, 
				 opt, ainfo->prn);
    } 

 bailout:

    if (err) {
	pmod->errcode = err;
    }

    kalman_free(K);
    free_kalman_matrices();

    gretl_matrix_free(y);
    free_arma_X_matrix(ainfo, X);

    gretl_matrix_replace(&Svar2, NULL);
    gretl_matrix_replace(&vQ, NULL);

    free(b);

    /* unpublish ainfo */
    kainfo = NULL;

    return err;
}

/* end of Kalman-specific material */

static int user_arma_init (double *coeff, arma_info *ainfo, int *init_done)
{
    PRN *prn = ainfo->prn;
    int i, nc = n_init_vals();

    if (nc == 0) {
	return 0;
    } else if (nc < ainfo->nc) {
	pprintf(prn, "arma initialization: need %d coeffs but got %d\n",
		ainfo->nc, nc);
	return E_DATA;
    }

    if (ainfo->flags & ARMA_EXACT) {
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

/* construct a "virtual dataset" in the form of a set of pointers into
   the main dataset: this will be passed to the bhhh_max function.
   The dependent variable is put in position 0; following this are the
   independent variables.
*/

static const double **make_armax_X (arma_info *ainfo, const double **Z)
{
    const double **X;
    int *list = ainfo->alist;
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
    if (ainfo->y != NULL) {
	X[0] = ainfo->y;
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

/* add extra OP-related info to the arma info struct */

static int set_up_arma_OPG_info (arma_info *ainfo, 
				 const double **Z,
				 const DATAINFO *pdinfo)
{
    /* array length needed for derivatives */
    int nd = 1 + ainfo->q + ainfo->pd * ainfo->Q;
    /* number of derivatives */
    int k = ainfo->nc;
    int err = 0;

    /* construct virtual dataset for dep var, real regressors */
    ainfo->X = make_armax_X(ainfo, Z);
    if (ainfo->X == NULL) {
	err = E_ALLOC;
    }  

    if (!err) {
	/* allocate gradient matrix */
	ainfo->G = gretl_zero_matrix_new(ainfo->T, k);
	if (ainfo->G == NULL) {
	    err = E_ALLOC;
	}
    }    

    if (!err && !(ainfo->flags & ARMA_EXACT)) {
	/* allocate covariance matrix */
	ainfo->V = gretl_matrix_alloc(k, k);
	if (ainfo->V == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* forecast errors array */
	ainfo->e = malloc((ainfo->t2 + 1) * sizeof *ainfo->e);
	if (ainfo->e == NULL) {
	    err = E_ALLOC;
	} else {
	    int t;

	    for (t=0; t<=ainfo->t2; t++) {
		ainfo->e[t] = 0.0;
	    }
	}
    }

    if (!err) {
	/* derivatives arrays */
	ainfo->aux = doubles_array_new0(k, nd); 
	if (ainfo->aux == NULL) {
	    err = E_ALLOC;
	} else {
	    ainfo->n_aux = k;
	}
    }

    return err;
}

/* retrieve results specific to bhhh procedure */

static int 
conditional_arma_model_prep (MODEL *pmod, arma_info *ainfo,
			     double *theta)
{
    int i, t, err;

    pmod->t1 = ainfo->t1;
    pmod->t2 = ainfo->t2;
    pmod->nobs = pmod->t2 - pmod->t1 + 1;
    pmod->ncoeff = ainfo->nc;

    err = gretl_model_allocate_storage(pmod);
    if (err) {
	return err;
    }

    pmod->lnL = ainfo->ll;
    pmod->sigma = NADBL; /* will be replaced */

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = theta[i];
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->uhat[t] = ainfo->e[t];
    }

    err = gretl_model_write_vcv(pmod, ainfo->V);

    return err;
}

static int bhhh_arma (double *theta, 
		      const double **Z, const DATAINFO *pdinfo,
		      arma_info *ainfo, MODEL *pmod,
		      gretlopt opt)
{
    gretlopt bhhh_opt = OPT_NONE;
    double tol = libset_get_double(BHHH_TOLER);
    int iters, err = 0;

    err = set_up_arma_OPG_info(ainfo, Z, pdinfo);
    if (err) {
	pmod->errcode = err;
	return err;
    }

    if (opt & OPT_V) {
	bhhh_opt |= OPT_V;
    }

    err = bhhh_max(theta, ainfo->nc, ainfo->G,
		   bhhh_arma_callback, tol, &iters,
		   ainfo, ainfo->V, bhhh_opt, ainfo->prn);
    
    if (err) {
	fprintf(stderr, "arma: bhhh_max returned %d\n", err);
    } else {
	pmod->full_n = pdinfo->n;
	err = conditional_arma_model_prep(pmod, ainfo, theta);
    }

    if (!err) {
	gretl_model_set_int(pmod, "iters", iters);
	write_arma_model_stats(pmod, ainfo, Z, pdinfo);
	arma_model_add_roots(pmod, ainfo, theta);
    }

    if (err && !pmod->errcode) {
	pmod->errcode = err;
    }

    return pmod->errcode;
}

/* Should we try Hannan-Rissanen initialization of ARMA
   coefficients? */

static int prefer_hr_init (arma_info *ainfo)
{
    int ret = 0;

    if (ainfo->q > 1 || ainfo->Q > 0) {
	ret = 1;
	if (ainfo->pqspec != NULL && *ainfo->pqspec != '\0') {
	    /* don't use for gappy arma (yet?) */
	    ret = 0;
	} else if (ainfo->flags & ARMA_XDIFF) {
	    /* don't use for ARIMAX (yet?) */
	    ret = 0;
	} else if (ainfo->t2 - ainfo->t1 < 100) {
	    /* unlikely to work well with small sample */
	    ret = 0;
	} else if (ainfo->p > 0 && ainfo->P > 0) {
	    /* not sure about this: HR catches the MA terms, but NLS
	       handles the AR interactions better
	    */
	    ret = 0;
	} else if ((ainfo->P > 0 && ainfo->p >= ainfo->pd) ||
		   (ainfo->Q > 0 && ainfo->q >= ainfo->pd)) {
	    /* overlapping orders screw things up */
	    ret = 0;
	} else if (ret && arma_exact_ml(ainfo)) {
	    /* screen for cases where we'll use NLS */
	    if (ainfo->P > 0) {
		ret = 0;
	    } else if (ainfo->p + ainfo->P > 0 && ainfo->nexo > 0) {
		ret = 0;
	    }
	}
    }

#if ARMA_DEBUG
    fprintf(stderr, "prefer_hr_init? %s\n", ret? "yes" : "no");
#endif

    return ret;
}

static int arma_via_OLS (arma_info *ainfo, const double *coeff, 
			 const double **Z, const DATAINFO *pdinfo,
			 MODEL *pmod)
{
    int err = 0;

    ainfo->flags |= ARMA_LS;

    err = arma_by_ls(coeff, Z, pdinfo, ainfo, pmod);

    if (!err) {
	pmod->t1 = ainfo->t1;
	pmod->t2 = ainfo->t2;
	pmod->full_n = pdinfo->n;
	write_arma_model_stats(pmod, ainfo, Z, pdinfo);
	arma_model_add_roots(pmod, ainfo, pmod->coeff);
	gretl_model_set_int(pmod, "arma_flags", ARMA_LS);
    }

    return err;
}

/* Set flag to indicate differencing of exogenous regressors, in the
   case of an ARIMAX model using native exact ML -- unless this is
   forbidden by OPT_Y.
*/

static void maybe_set_xdiff_flag (arma_info *ainfo, gretlopt opt)
{
    if ((ainfo->flags & ARMA_EXACT) &&
	(ainfo->d > 0 || ainfo->D > 0) &&
	ainfo->nexo > 0 && !(opt & OPT_Y)) {
	ainfo->flags |= ARMA_XDIFF;
    }
}

MODEL arma_model (const int *list, const char *pqspec,
		  const double **Z, const DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn)
{
    double *coeff = NULL;
    MODEL armod;
    arma_info ainfo;
    int init_done = 0;
    int err = 0;

    arma_info_init(&ainfo, opt, pqspec, pdinfo);
    ainfo.prn = set_up_verbose_printer(opt, prn);

    gretl_model_init(&armod); 

    err = incompatible_options(opt, OPT_C | OPT_L);

    if (!err) {
	ainfo.alist = gretl_list_copy(list);
	if (ainfo.alist == NULL) {
	    err = E_ALLOC;
	}
    } 

    if (!err) {
	err = arma_check_list(&ainfo, Z, pdinfo, opt);
    }

    if (!err) {
	/* calculate maximum lag */
	maybe_set_xdiff_flag(&ainfo, opt);
	calc_max_lag(&ainfo);
    }

    if (!err) {
	/* adjust sample range if need be */
	err = arma_adjust_sample(&ainfo, Z, pdinfo);
    }

    if (!err) {
	/* allocate initial coefficient vector */
	coeff = malloc(ainfo.nc * sizeof *coeff);
	if (coeff == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && (ainfo.d > 0 || ainfo.D > 0)) {
	err = arima_difference(&ainfo, Z, pdinfo);
    }

    if (err) {
	goto bailout;
    }

    /* initialize the coefficients: there are 3 possible methods */

    /* first pass: see if the user specified some values */
    err = user_arma_init(coeff, &ainfo, &init_done);
    if (err) {
	goto bailout;
    }

    if (!(ainfo.flags & ARMA_EXACT) && ainfo.q == 0 && ainfo.Q == 0) {
	/* pure AR model can be estimated via least squares */
	const double *b = (init_done)? coeff : NULL;

	err = arma_via_OLS(&ainfo, b, Z, pdinfo, &armod);
	goto bailout;
    }

    /* second pass: try Hannan-Rissanen, if suitable */
    if (!init_done && prefer_hr_init(&ainfo)) {
	hr_arma_init(coeff, Z, pdinfo, &ainfo, &init_done);
    }

    /* third pass: estimate pure AR model by OLS or NLS */
    if (!err && !init_done) {
	err = ar_arma_init(coeff, Z, pdinfo, &ainfo, &armod);
    }

    if (!err) {
	clear_model_xpx(&armod);
	if (ainfo.flags & ARMA_EXACT) {
	    kalman_arma(coeff, Z, pdinfo, &ainfo, &armod, opt);
	} else {
	    bhhh_arma(coeff, Z, pdinfo, &ainfo, &armod, opt);
	}
    }

 bailout:

    if (err && !armod.errcode) {
	armod.errcode = err;
    }

    if (armod.errcode) {
	if (opt & OPT_U) {
	    armod.opt |= OPT_U; /* continue on error */
	}
    } else {
	gretl_model_smpl_init(&armod, pdinfo);
    }

    free(coeff);
    arma_info_cleanup(&ainfo);

    /* cleanup in MA roots checker */
    bounds_checker_cleanup();

    if (ainfo.prn != prn) {
	close_down_verbose_printer(ainfo.prn);
    }

    return armod;
}
