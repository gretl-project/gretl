/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2006 Allin Cottrell and Riccardo "Jack" Lucchetti
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

/* Native gretl code for ARMA estimation.  Much of the code here
   was contributed by Riccardo "Jack" Lucchetti, the rest is due
   to Allin Cottrell; thanks also to Stephen Moshier for cephes.
*/

#include "../cephes/polrt.c"

#include "libgretl.h"
#include "bhhh_max.h"
#include "libset.h"
#include "kalman.h"
#include "matrix_extra.h"

#define ARMA_DEBUG 0

/* ln(sqrt(2*pi)) + 0.5 */
#define LN_SQRT_2_PI_P5 1.41893853320467274178

#define KALMAN_ALL 999

#include "arma_common.c"

static PRN *errprn;

/* check whether the MA estimates have gone out of bounds in the
   course of iteration */

static int 
ma_out_of_bounds (struct arma_info *ainfo, const double *theta,
		  const double *Theta)
{
    static double *temp;
    static double *tmp2;
    static cmplx *roots;
    static int qmax;

    double re, im, rt;
    int i, j, k, qtot;
    int tzero = 1, Tzero = 1;
    int err = 0, cerr = 0;

    if (ainfo == NULL) {
	/* signal for cleanup */
	free(temp);
	temp = NULL;
	free(tmp2);
	tmp2 = NULL;
	free(roots);
	roots = NULL;
	qmax = 0;
	return 0;
    }

    for (i=0; i<ainfo->q && tzero; i++) {
	if (theta[i] != 0.0) {
	    tzero = 0;
	}    
    }  

    for (i=0; i<ainfo->Q && Tzero; i++) {
	if (Theta[i] != 0.0) {
	    Tzero = 0;
	}    
    }  
    
    if (tzero && Tzero) {
	return 0;
    }

    if (temp == NULL) {
	/* not allocated yet */
	qmax = ainfo->q + ainfo->Q * ainfo->pd;

	temp  = malloc((qmax + 1) * sizeof *temp);
	tmp2  = malloc((qmax + 1) * sizeof *tmp2);
	roots = malloc(qmax * sizeof *roots);

	if (temp == NULL || tmp2 == NULL || roots == NULL) {
	    free(temp);
	    temp = NULL;
	    free(tmp2);
	    tmp2 = NULL;
	    free(roots);
	    roots = NULL;
	    return 1;
	}
    }

    temp[0] = 1.0;

    /* initialize to non-seasonal MA or zero */
    for (i=0; i<qmax; i++) {
        if (i < ainfo->q) {
            temp[i+1] = theta[i];
        } else {
            temp[i+1] = 0.0;
        }
    }

    /* add seasonal MA and interaction */
    if (Tzero) {
	qtot = ainfo->q;
    } else {
	qtot = qmax;
	for (i=0; i<ainfo->Q; i++) {
	    k = (i + 1) * ainfo->pd;
	    temp[k] += Theta[i];
	    for (j=0; j<ainfo->q; j++) {
		int m = k + j + 1;

		temp[m] += Theta[i] * theta[j];
	    }
	}
    }

    cerr = polrt(temp, tmp2, qtot, roots);
    if (cerr) {
	fprintf(stderr, "ma_out_of_bounds: polrt returned %d\n", cerr);
	return 0; /* ?? */
    }

    for (i=0; i<qtot; i++) {
	re = roots[i].r;
	im = roots[i].i;
	rt = re * re + im * im;
	if (rt > DBL_EPSILON && rt <= 1.0) {
	    pprintf(errprn, "MA root %d = %g\n", i, rt);
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
    int i, j, s, p;

    for (i=0; i<ainfo->q; i++) {
	s = t - (i + 1);
	if (s >= 0) {
	    drv[t] -= theta[i] * drv[s];
	}
    }

    for (i=0; i<ainfo->Q; i++) {
	s = t - ainfo->pd * (i + 1);
	if (s >= 0) {
	    drv[t] -= Theta[i] * drv[s];
	    for (j=0; j<ainfo->q; j++) {
		p = s - (j + 1);
		if (p >= 0) {
		    drv[t] -= Theta[i] * theta[j] * drv[p];
		}
	    }
	}
    }
}

/* Calculate ARMA log-likelihood.  This function is passed to the
   bhhh_max() routine as a "callback". */

static int arma_ll (double *coeff, 
		    const double **bhX, double **Z, 
		    model_info *arma,
		    int do_score)
{
    int i, j, s, t;
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
    Phi = phi + ainfo->p;
    theta = Phi + ainfo->P;
    Theta = theta + ainfo->q;
    beta = Theta + ainfo->Q;

    /* pointers to blocks of derivatives */
    de_a = de + ainfo->ifc;
    de_sa = de_a + ainfo->p;
    de_m = de_sa + ainfo->P;
    de_sm = de_m + ainfo->q;
    de_r = de_sm + ainfo->Q;

#if ARMA_DEBUG
    fprintf(stderr, "arma_ll: p=%d, q=%d, P=%d, Q=%d, pd=%d\n",
	    ainfo->p, ainfo->q, ainfo->P, ainfo->Q, ainfo->pd);
#endif

    if (ma_out_of_bounds(ainfo, theta, Theta)) {
	pputs(errprn, "arma: MA estimate(s) out of bounds\n");
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
	for (i=0; i<ainfo->p; i++) {
	    s = t - (i + 1);
	    e[t] -= phi[i] * y[s];
	}

	/* seasonal AR component plus interactions */
	for (i=0; i<ainfo->P; i++) {
	    s = t - ainfo->pd * (i + 1);
	    e[t] -= Phi[i] * y[s];
	    for (j=0; j<ainfo->p; j++) {
		p = s - (j + 1);
		e[t] += Phi[i] * phi[j] * y[p];
	    }
	}

	/* non-seasonal MA component */
	for (i=0; i<ainfo->q; i++) {
	    s = t - (i + 1);
	    if (s >= t1) {
		e[t] -= theta[i] * e[s];
	    }
	}

	/* seasonal MA component plus interactions */
	for (i=0; i<ainfo->Q; i++) {
	    s = t - ainfo->pd * (i + 1);
	    if (s >= t1) {
		e[t] -= Theta[i] * e[s];
		for (j=0; j<ainfo->q; j++) {
		    p = s - (j + 1);
		    if (p >= t1) {
			e[t] -= Theta[i] * theta[j] * e[p];
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
	    for (j=0; j<ainfo->p; j++) {
		lag = j + 1;
		if (t >= lag) {
		    de_a[j][t] = -y[t-lag];
		    /* cross-partial with seasonal AR */
		    for (i=0; i<ainfo->P; i++) {
			xlag = lag + ainfo->pd * (i + 1);
			if (t >= xlag) {
			    de_a[j][t] += Phi[i] * y[t-xlag];
			}
		    }
		    do_MA_partials(de_a[j], ainfo, theta, Theta, t);
		}
	    }

	    /* seasonal AR terms (de_sa) */
	    for (j=0; j<ainfo->P; j++) {
		lag = ainfo->pd * (j + 1);
		if (t >= lag) {
		    de_sa[j][t] = -y[t-lag];
		    /* cross-partial with non-seasonal AR */
		    for (i=0; i<ainfo->p; i++) {
			xlag = lag + (i + 1);
			if (t >= xlag) {
			    de_sa[j][t] += phi[i] * y[t-xlag];
			}
		    }
		    do_MA_partials(de_sa[j], ainfo, theta, Theta, t);
		}
	    }

	    /* non-seasonal MA terms (de_m) */
	    for (j=0; j<ainfo->q; j++) {
		lag = j + 1;
		if (t >= lag) {
		    de_m[j][t] = -e[t-lag];
		    /* cross-partial with seasonal MA */
		    for (i=0; i<ainfo->Q; i++) {
			xlag = lag + ainfo->pd * (i + 1);
			if (t >= xlag) {
			    de_m[j][t] -= Theta[i] * e[t-xlag];
			}
		    }
		    do_MA_partials(de_m[j], ainfo, theta, Theta, t);
		}
	    }

	    /* seasonal MA terms (de_sm) */
	    for (j=0; j<ainfo->Q; j++) {
		lag = ainfo->pd * (j + 1);
		if (t >= lag) {
		    de_sm[j][t] = -e[t-lag];
		    /* cross-partial with non-seasonal MA */
		    for (i=0; i<ainfo->q; i++) {
			xlag = lag + (i + 1);
			if (t >= xlag) {
			    de_sm[j][t] -= theta[i] * e[t-xlag];
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
    const double *Phi = phi + ainfo->p;
    const double *theta = Phi + ainfo->P;
    const double *Theta = theta + ainfo->q;

    int nr = ainfo->p + ainfo->P + ainfo->q + ainfo->Q;
    int pmax, qmax, lmax;
    double *temp = NULL, *temp2 = NULL;
    cmplx *rptr, *roots = NULL;
    int i;

    pmax = (ainfo->p > ainfo->P)? ainfo->p : ainfo->P;
    qmax = (ainfo->q > ainfo->Q)? ainfo->q : ainfo->Q;
    lmax = (pmax > qmax)? pmax : qmax;

    if (pmax == 0 && qmax == 0) {
	return 0;
    }

    temp  = malloc((lmax + 1) * sizeof *temp);
    temp2 = malloc((lmax + 1) * sizeof *temp2);
    roots = malloc(nr * sizeof *roots);

    if (temp == NULL || temp2 == NULL || roots == NULL) {
	free(temp);
	free(temp2);
	free(roots);
	return E_ALLOC;
    }

    temp[0] = 1.0;
    rptr = roots;

    if (ainfo->p > 0) {
	/* A(L), non-seasonal */
	for (i=0; i<ainfo->p; i++) {
	    temp[i+1] = -phi[i];
	}
	polrt(temp, temp2, ainfo->p, rptr);
	rptr += ainfo->p;
    }

    if (ainfo->P > 0) {
	/* B(L), seasonal */
	for (i=0; i<ainfo->P; i++) {
	    temp[i+1] = -Phi[i];
	}    
	polrt(temp, temp2, ainfo->P, rptr);
	rptr += ainfo->P;
    }

    if (ainfo->q > 0) {
	/* C(L), non-seasonal */
	for (i=0; i<ainfo->q; i++) {
	    temp[i+1] = theta[i];
	}  
	polrt(temp, temp2, ainfo->q, rptr);
	rptr += ainfo->q;
    }

    if (ainfo->Q > 0) {
	/* D(L), seasonal */
	for (i=0; i<ainfo->Q; i++) {
	    temp[i+1] = Theta[i];
	}  
	polrt(temp, temp2, ainfo->Q, rptr);
    }
    
    free(temp);
    free(temp2);

    gretl_model_set_data(pmod, "roots", roots, MODEL_DATA_CMPLX_ARRAY,
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
    int i, j, k;

    for (i=0; i<=pmax; i++) {
	ac[i] = 0.0;
    }

    for (i=0; i<=ainfo->P; i++) {
	x = (i == 0)? -1 : Phi[i-1];
	for (j=0; j<=ainfo->p; j++) {
	    y = (j == 0)? -1 : phi[j-1];
	    k = j + ainfo->pd * i;
	    ac[k] -= x * y;
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
    int i, j, k;

    for (i=0; i<=qmax; i++) {
	mc[i] = 0.0;
    }

    for (i=0; i<=ainfo->Q; i++) {
	x = (i == 0)? 1 : Theta[i-1];
        for (j=0; j<=ainfo->q; j++) {
	    y = (j == 0)? 1 : theta[j-1];
            k = j + ainfo->pd * i;
	    mc[k] = x * y;
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
    const double *Phi = phi + kainfo->p;
    const double *theta = Phi + kainfo->P;
    const double *Theta = theta + kainfo->q;
    const double *beta = Theta + kainfo->Q;
    double mu = (kainfo->ifc)? b[0] : 0.0;
    int i, err = 0;

    int rewrite_A = 0;
    int rewrite_F = 0;
    int rewrite_H = 0;

    gretl_matrix_zero(S);

    if (idx == KALMAN_ALL) {
	rewrite_A = rewrite_F = rewrite_H = 1;
    } else {
	int pmax = kainfo->ifc + kainfo->p + kainfo->P;
	int tmax = pmax + kainfo->q + kainfo->Q;

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
	    for (i=0; i<kainfo->q; i++) {
		gretl_vector_set(H, i + 1, theta[i]);
	    }
	}
    }

    if (rewrite_F) {

	/* form the F matrix using phi and/or Phi */
	if (kainfo->P > 0) {
	    write_big_phi(phi, Phi, kainfo, F);
	} else {
	    for (i=0; i<kainfo->p; i++) {
		gretl_matrix_set(F, 0, i, phi[i]);
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
    int t, i = 0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->uhat[t] = gretl_vector_get(E, i++);
    }

    pmod->vcv = vcv;

    for (i=0; i<k; i++) {
	pmod->sderr[i] = sqrt(vcv[ijton(i, i, k)]);
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
    int i;

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

#if ARMA_DEBUG
    fprintf(stderr, "kalman_arma_finish: doing VCV, method %s\n",
	    (vcv != NULL)? "Hessian" : "OPG");
#endif

    if (vcv != NULL) {
	arma_hessian_vcv(pmod, vcv, k);
    } else {
	arma_OPG_vcv(pmod, K, b, s2, k, T);
    }

    write_arma_model_stats(pmod, alist, ainfo, Z, pdinfo);
    arma_model_add_roots(pmod, ainfo, b);

    gretl_model_set_int(pmod, "arma_flags", ARMA_EXACT);

    return pmod->errcode;
}

static int kalman_do_ma_check = 1;

static double kalman_arma_ll (const double *b, void *p)
{
    int offset = kainfo->ifc + kainfo->p + kainfo->P;
    const double *theta = b + offset;
    const double *Theta = theta + kainfo->q;
    double ll = NADBL;
    kalman *K;
    int err = 0;

#if ARMA_DEBUG
    int i;
    fprintf(stderr, "kalman_arma_ll():\n");
    for (i=0; i<kainfo->q; i++) {
	fprintf(stderr, "theta[%d] = %#.12g\n", i, theta[i]);
    }
    for (i=0; i<kainfo->Q; i++) {
	fprintf(stderr, "Theta[%d] = %#.12g\n", i, Theta[i]);
    }    
#endif

    if (kalman_do_ma_check && ma_out_of_bounds(kainfo, theta, Theta)) {
	pputs(errprn, "arma: MA estimate(s) out of bounds\n");
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

    x = gretl_matrix_data_subset(xlist, Z, ainfo->t1, ainfo->t2, NULL);
    if (x == NULL) {
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
   wanted for initializing the Kalman filter
*/

static void transform_arma_const (double *b, struct arma_info *ainfo)
{
    const double *phi = b + 1;
    const double *Phi = phi + ainfo->p;
    double narfac = 1.0;
    double sarfac = 1.0;
    int i;

    for (i=0; i<ainfo->p; i++) {
	narfac -= phi[i];
    }

    for (i=0; i<ainfo->P; i++) {
	sarfac -= Phi[i];
    }

    b[0] /= (narfac * sarfac);
}

static int kalman_arma (const int *alist, double *coeff, 
			const double **Z, const DATAINFO *pdinfo,
			struct arma_info *ainfo, MODEL *pmod,
			PRN *prn)
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
    double *hess = NULL;
    int i, T, err = 0;

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
	fprintf(stderr, " b[%d] = %g\n", i, b[i]);
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
	if (Svar2 == NULL) {
	    err = E_ALLOC;
	}
    } else {
	vQ = gretl_column_vector_alloc(r * r);
    }

    if (err || S == NULL || P == NULL || F == NULL || A == NULL ||
	H == NULL || E == NULL || Q == NULL || Svar == NULL || 
	vQ == NULL) {
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

    K = kalman_new(S, P, F, A, H, Q, NULL, y, x, E, ainfo->nc, ainfo->ifc, &err);

    if (err) {
	fprintf(stderr, "kalman_new(): err = %d\n", err);
    } else {
	if (r > 3) {
	    kalman_set_nonshift(K, 1);
	} else {
	    kalman_set_nonshift(K, r);
	}
	kalman_use_ARMA_ll(K);
	err = BFGS_max(b, ainfo->nc, maxit, reltol, 
		       &fncount, &grcount, kalman_arma_ll, C_LOGLIK,
		       NULL, K, (prn != NULL)? OPT_V : OPT_NONE, prn);
	if (err) {
	    fprintf(stderr, "BFGS_max returned %d\n", err);
	} else if (1 || getenv("ARMA_HESSIAN")) { /* Hello! Testing */
	    kalman_do_ma_check = 0;
	    hess = numerical_hessian(b, ainfo->nc, kalman_arma_ll, K, &err);
	    kalman_do_ma_check = 1;
	}
    }

    if (err) {
	pmod->errcode = err;
    } else {
	gretl_model_set_int(pmod, "fncount", fncount);
	gretl_model_set_int(pmod, "grcount", grcount);
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

/* for ARMAX: write the component of the NLS specification
   that takes the form (y_{t-i} - X_{t-i} \beta)
*/

static void y_Xb_at_lag (char *spec, struct arma_info *ainfo, 
			 int narmax, int lag)
{
    char term[32];
    int i, nt;

    if (narmax == 0) {
	sprintf(term, "y_%d", lag);
	strcat(spec, term);
	return;
    }

    nt = ainfo->ifc + narmax;

    sprintf(term, "(y_%d-", lag);
    strcat(spec, term);

    if (nt > 1) {
	strcat(spec, "(");
    }

    if (ainfo->ifc) {
	strcat(spec, "b0");
    }

    for (i=1; i<=narmax; i++) {
	if (ainfo->ifc || i > 1) {
	    strcat(spec, "+");
	} 
	sprintf(term, "b%d*x%d_%d", i, i, lag);
	strcat(spec, term);
    }

    if (nt > 1) {
	strcat(spec, "))");
    } else {
	strcat(spec, ")");
    }
}

static int arma_get_nls_model (MODEL *amod, struct arma_info *ainfo,
			       int narmax, double ***pZ, DATAINFO *pdinfo) 
{
#if ARMA_DEBUG
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
    int nparam;
    int i, j, k, err = 0;

    spec = nlspec_new(NLS, pdinfo);
    if (spec == NULL) {
	return E_ALLOC;
    }

    nlspec_set_t1_t2(spec, 0, ainfo->t2 - ainfo->t1); /* ?? */

    nparam = ainfo->ifc + ainfo->p + ainfo->P + ainfo->nexo;

    plist = gretl_list_new(nparam);
    if (plist == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = dataset_add_scalars(nparam, pZ, pdinfo); 
    if (err) {
	goto bailout;
    }

    /* construct names for the parameters, and param list;
       also do some initialization -- but FIXME in that regard
    */

    v = oldv;
    k = 1;
    if (ainfo->ifc) {
	(*pZ)[v][0] = gretl_mean(0, pdinfo->n - 1, (*pZ)[1]); /* ? */
	strcpy(pdinfo->varname[v], "b0");
	plist[k++] = v++;
    }
    for (i=1; i<=ainfo->p; i++) {
	if (i == 1) {
	    (*pZ)[v][0] = 0.1; /* ? */
	}
	sprintf(pdinfo->varname[v], "phi%d", i);
	plist[k++] = v++;
    }
    for (i=1; i<=ainfo->P; i++) {
	sprintf(pdinfo->varname[v], "Phi%d", i);
	plist[k++] = v++;
    }
    for (i=1; i<=ainfo->nexo; i++) {
	sprintf(pdinfo->varname[v], "b%d", i);
	plist[k++] = v++;
    }

    /* construct NLS specification */

    strcpy(fnstr, "y=");

    if (ainfo->ifc) {
	strcat(fnstr, "b0");
    } else {
	strcat(fnstr, "0");
    } 

    for (i=0; i<=ainfo->p; i++) {
	if (i > 0) {
	    sprintf(term, "+phi%d*", i);
	    strcat(fnstr, term);
	    y_Xb_at_lag(fnstr, ainfo, narmax, i);
	}
	for (j=0; j<=ainfo->P; j++) {
	    if (i == 0 && j > 0) {
		sprintf(term, "+Phi%d*", j);
		strcat(fnstr, term);
		y_Xb_at_lag(fnstr, ainfo, narmax, j * ainfo->pd);
	    }
	    if (i > 0 && j > 0) {
		sprintf(term, "-phi%d*Phi%d*", i, j);
		strcat(fnstr, term);
		y_Xb_at_lag(fnstr, ainfo, narmax, j * ainfo->pd + i);
	    }
	}
    } 

    for (i=1; i<=ainfo->nexo; i++) {
	sprintf(term, "+b%d*x%d", i, i);
	strcat(fnstr, term);
    }	

#if ARMA_DEBUG
    fprintf(stderr, "initting using NLS spec:\n %s\n", fnstr);
#endif

    err = nlspec_set_regression_function(spec, fnstr, pdinfo);

    if (!err) {
	err = nlspec_add_param_list(spec, plist, (const double **) *pZ,
				    pdinfo);
    }

    if (!err) {
	*amod = model_from_nlspec(spec, pZ, pdinfo, nlsopt, prn);
	err = amod->errcode;
#if ARMA_DEBUG
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

static int *make_ar_ols_list (struct arma_info *ainfo, int av, int ptotal)
{
    int * alist = gretl_list_new(av);
    int i, offset;

    if (alist == NULL) {
	return NULL;
    }

    alist[1] = 1;

    if (ainfo->ifc) {
	alist[2] = 0;
	offset = 2;
    } else {
	alist[0] -= 1;
	offset = 1;
    }

    for (i=1; i<=ainfo->p; i++) {
	alist[offset + i] = 1 + i;
    }

    for (i=1; i<=ainfo->P; i++) {
	alist[offset + ainfo->p + i] = ainfo->p + 1 + i;
    }

    offset += ptotal;  

    for (i=1; i<=ainfo->nexo; i++) {
	alist[offset + i] = ptotal + 1 + i;
    }

    return alist;
}

/* Run a least squares model to get initial values for the AR
   coefficients, either OLS or NLS.  We use NLS if there is
   nonlinearity due to either (a) the presence of both a seasonal and
   a non-seasonal AR component or (b) the presence of exogenous
   variables in the context of a non-zero AR order, where estimation
   will be via exact ML.  In this initialization any MA coefficients
   are simply set to zero.
*/

static int ar_arma_init (const int *list, double *coeff, 
			 const double **Z, const DATAINFO *pdinfo,
			 struct arma_info *ainfo, PRN *prn)
{
    int an = pdinfo->t2 - ainfo->t1 + 1;
    int np = ainfo->p, nq = ainfo->q;
    int nmixed = ainfo->p * ainfo->P;
    int ptotal = ainfo->p + ainfo->P + nmixed;
    int av = ptotal + ainfo->nexo + 2;
    const double *y;
    double **aZ = NULL;
    DATAINFO *adinfo = NULL;
    int *alist = NULL;
    MODEL armod;
    int nonlin = 0;
    int narmax = 0;
    int xstart, lag;
    int axi = 0, ayi = 0;
    int i, j, k, t;
    int err = 0;

#if ARMA_DEBUG
    fprintf(stderr, "ar_arma_init: pdinfo->t1=%d, pdinfo->t2=%d (n=%d); "
	    "ainfo->t1=%d, ainfo->t2=%d\n",
	    pdinfo->t1, pdinfo->t2, pdinfo->n, ainfo->t1, ainfo->t2);
#endif

    /* dependent variable */
    if (ainfo->dy != NULL) {
	y = ainfo->dy;
    } else {
	y = Z[ainfo->yno];
    }

    if (ptotal == 0 && ainfo->nexo == 0 && !ainfo->ifc) {
	/* special case of pure MA model */
	for (i=0; i<nq; i++) {
	    coeff[i] = 0.0; 
	} 
	return 0;
    }

    if (arma_exact_ml(ainfo)) {
	narmax = ainfo->nexo;
	if (narmax > 0) {
	    /* ARMAX-induced lags of exog vars */
	    av += ainfo->nexo * ptotal;
	}
    }

    gretl_model_init(&armod); 

    adinfo = create_new_dataset(&aZ, av, an, 0);
    if (adinfo == NULL) {
	return E_ALLOC;
    }

    if (ptotal > 0 && (narmax > 0 || nmixed > 0)) {
	/* have to use NLS */
	nonlin = 1;
    } else {
	/* OLS: need regression list */
	alist = make_ar_ols_list(ainfo, av, ptotal);
    }

    /* starting position for reading exogeneous vars */
    if (ainfo->d > 0 || ainfo->D > 0) {
	xstart = (arma_has_seasonal(ainfo))? 10 : 6;
    } else {
	xstart = (arma_has_seasonal(ainfo))? 8 : 5;
    }

    /* construct the variable names */

    strcpy(adinfo->varname[1], "y");

    axi = ptotal + ainfo->nexo + 2;

    for (i=1; i<=ainfo->p; i++) {
	sprintf(adinfo->varname[i+1], "y_%d", i);
	for (j=1; j<=narmax; j++) {
	    sprintf(adinfo->varname[axi++], "x%d_%d", j, i);
	}
    }

    ayi = ainfo->p + ainfo->P + 2;

    for (i=1; i<=ainfo->P; i++) {
	k = ainfo->p + 1 + i;
	sprintf(adinfo->varname[k], "y_%d", ainfo->pd * i);
	for (j=1; j<=narmax; j++) {
	    sprintf(adinfo->varname[axi++], "x%d_%d", j, ainfo->pd * i);
	}
	for (j=1; j<=ainfo->p; j++) {
	    sprintf(adinfo->varname[ayi++], "y_%d", ainfo->pd * i + j);
	    for (k=1; k<=narmax; k++) {
		sprintf(adinfo->varname[axi++], "x%d_%d", k, ainfo->pd * i + j);
	    }
	}
    }

    axi = ptotal + 2;

    for (i=1; i<=ainfo->nexo; i++) {
	sprintf(adinfo->varname[axi++], "x%d", i);
    }

    /* Build temporary dataset including lagged vars: if we're doing
       exact ML on an ARMAX model we need lags of the exogenous
       variables as well as lags of y_t.  Note that the auxiliary
       dataset has "t = 0" at an offset of ainfo->t1 into the "real",
       external dataset.
    */

    for (t=0; t<an; t++) {
	int realt = t + ainfo->t1;
	int miss = 0;
	int s, m;

	aZ[1][t] = y[realt];

	axi = ptotal + ainfo->nexo + 2;

	for (i=1; i<=ainfo->p; i++) {
	    s = realt - i;
	    if (s < 0) miss = 1;
	    aZ[i+1][t] = (s >= 0)? y[s] : NADBL;
	    for (j=1; j<=narmax; j++) {
		m = list[xstart + j - 1];
		aZ[axi++][t] = (s >= 0)? Z[m][s] : NADBL;
	    }
	}

	ayi = ainfo->p + ainfo->P + 2;

	for (i=1; i<=ainfo->P; i++) {
	    lag = ainfo->pd * i;
	    s = realt - lag;
	    if (s < 0) miss = 1;
	    k = ainfo->p + 1 + i;
	    aZ[k][t] = (s >= 0)? y[s] : NADBL;
	    for (k=1; k<=narmax; k++) {
		m = list[xstart + k - 1];
		aZ[axi++][t] = (s >= 0)? Z[m][s] : NADBL;
	    }
	    for (j=1; j<=ainfo->p; j++) {
		lag = ainfo->pd * i + j;
		s = realt - lag;
		if (s < 0) miss = 1;
		aZ[ayi++][t] = (s >= 0)? y[s] : NADBL;
		for (k=1; k<=narmax; k++) {
		    m = list[xstart + k - 1];
		    aZ[axi++][t] = (s >= 0)? Z[m][s] : NADBL;
		}
	    }
	}

	axi = ptotal + 2;

	for (i=1; i<=ainfo->nexo; i++) {
	    m = list[xstart + i - 1];
	    aZ[axi++][t] = Z[m][realt];
	}

	if (miss) {
	    adinfo->t1 = t + 1;
	}	
    }

    if (arma_has_seasonal(ainfo)) {
	np += ainfo->P;
	nq += ainfo->Q;
    }

#if ARMA_DEBUG
    fprintf(stderr, "arma init dataset:\n");
    for (i=0; i<adinfo->v; i++) {
	fprintf(stderr, "var %d '%s', obs[0] = %g\n", i, adinfo->varname[i], 
		aZ[i][0]);
    }
#endif

    if (nonlin) {
#if ARMA_DEBUG
	fprintf(stderr, "arma:_init_by_ls: doing NLS\n");
#endif
	err = arma_get_nls_model(&armod, ainfo, narmax, &aZ, adinfo);
    } else {
#if ARMA_DEBUG
	printlist(alist, "'alist' in ar_arma_init (OLS)");
#endif
	armod = lsq(alist, &aZ, adinfo, OLS, OPT_A | OPT_Z);
	err = armod.errcode;
    }

    if (!err) {
	j = 0;
	for (i=0; i<armod.ncoeff; i++) {
	    if (i == np + ainfo->ifc) {
		j += nq; /* reserve space for MA coeffs */
	    }
	    coeff[j++] = armod.coeff[i];
	}
	for (i=0; i<nq; i++) {
	    /* insert near-zeros for MA coeffs */
	    coeff[i + np + ainfo->ifc] = 0.0001;
	} 
    }

    if (!err && arma_exact_ml(ainfo) && ainfo->ifc) {
	if (!nonlin || ainfo->nexo == 0) {
	    /* handle the case where we need to translate from an
	       estimate of the regression constant to the
	       unconditional mean of y_t
	    */
	    transform_arma_const(coeff, ainfo);
	}
    }

#if ARMA_DEBUG
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

    if (!err && prn != NULL) {
	if (nonlin) {
	    pputs(prn, "\narma initialization: using nonlinear AR model\n\n");
	} else {
	    pputs(prn, "\narma initialization: using linear AR model\n\n");
	}
    }

    /* clean up */
    free(alist);
    destroy_dataset(aZ, adinfo);
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

#if ARMA_DEBUG
    fprintf(stderr, "hr_init_check: ncoeff=%d, nobs=%d, 'df'=%d\n", 
	    ncoeff, nobs, df);
#endif

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
    double *theta = NULL;
    double *Theta = NULL;
    MODEL armod;
    int xstart;
    int m, pos, s;
    int i, j, t;
    int err = 0;

    pass1lags = (nQ + nP) * pdinfo->pd;
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

    adinfo = create_new_dataset(&aZ, pass1v + qtotal, an, 0);
    if (adinfo == NULL) {
	return E_ALLOC;
    }

#if ARMA_DEBUG
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

#if ARMA_DEBUG
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
	    for (i=0; i<nQ; i++) {
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
	    for (i=0; i<nP; i++) {
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
	pass2list[pos++] = arlags[i] + nexo + 1;
    }
    for (i=0; i<qtotal; i++) {
	pass2list[pos++] = pass1v + i;
    }
    
    /* now do pass2 */
    clear_model(&armod);
    armod = lsq(pass2list, &aZ, adinfo, OLS, OPT_A);

    if (armod.errcode) {
	err = armod.errcode;
    } else {
	/* rearrange coefficients */
	int pos2 = nexo + ainfo->ifc;

	if (ainfo->ifc) {
	    coeff[0] = armod.coeff[0];
	    pos = 1;
	} else {
	    pos = 0;
	}

#if ARMA_DEBUG
	PRN *prn = gretl_print_new(GRETL_PRINT_STDERR);
	printmodel(&armod, adinfo, OPT_S, prn);
	gretl_print_destroy(prn);
#endif

	for (i=0; i<np; i++) { /* phi */
#if ARMA_DEBUG
	    fprintf(stderr, "phi[%d] = coeff[%d] = %g\n", i+1, pos2, armod.coeff[pos2]);
#endif
	    coeff[pos++] = armod.coeff[pos2++];
	}
	for (i=0; i<nP; i++) { /* Phi */
#if ARMA_DEBUG
	    fprintf(stderr, "Phi[%d] = coeff[%d] = %g\n", i+1, pos2, armod.coeff[pos2]);
#endif
	    coeff[pos++] = armod.coeff[pos2];
	    pos2 += np + 1;
	}
	for (i=0; i<nq; i++) { /* theta */
#if ARMA_DEBUG
	    fprintf(stderr, "theta[%d] = coeff[%d] = %g\n", i+1, pos2, armod.coeff[pos2]);
#endif
	    if (i == 0) {
		theta = armod.coeff + pos2;
	    }
	    coeff[pos++] = armod.coeff[pos2++];
	}
	for (i=0; i<nQ; i++) { /* Theta */
#if ARMA_DEBUG
	    fprintf(stderr, "Theta[%d] = coeff[%d] = %g\n", i+1, pos2, armod.coeff[pos2]);
#endif
	    if (i == 0) {
		Theta = armod.coeff + pos2;
	    }
	    coeff[pos++] = armod.coeff[pos2];
	    pos2 += nq + 1;
	}
	for (i=0; i<nexo; i++) { /* exog */
#if ARMA_DEBUG
	    fprintf(stderr, "beta[%d] = coeff[%d] = %g\n", i+1, i+1, armod.coeff[i+1]);
#endif
	    coeff[pos++] = armod.coeff[i+1];
	}
    }

    bailout:

    free(pass1list);
    free(pass2list);
    free(arlags);
    free(malags);
    destroy_dataset(aZ, adinfo);
    clear_model(&armod);

    /* check MA values? */
    if (!err && (ainfo->q > 0 || ainfo->Q > 0)) {
	err = ma_out_of_bounds(ainfo, theta, Theta);
	bounds_checker_cleanup();
    }

    if (!err && prn != NULL) {
	pputs(prn, "\narma initialization: using Hannan-Rissanen method\n\n");
    }

    return err;
}

static int user_arma_init (double *coeff, struct arma_info *ainfo, 
			   char flags, int *init_done, PRN *prn)
{
    int i, nc = n_init_vals();

    if ((nc > 0) && (nc < ainfo->nc)) {
	pprintf(prn, "arma initialization: need %d coeffs but got %d\n",
		ainfo->nc, nc);
	return E_DATA;
    }

    if (!(flags & ARMA_EXACT)) {
	/* in EXACT case, this is handled within BFGS */
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
    double tol = get_bhhh_toler();
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
		      PRN *prn)
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
    err = bhhh_max(arma_ll, X, coeff, minfo, prn);
    
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

/* Should we try Hannan-Rissanen initialization of ARMA
   coefficients? */

static int prefer_hr_init (struct arma_info *ainfo)
{
    int ret = 0;

    /* unlikely to work well with small sample */
    if (ainfo->t2 - ainfo->t1 < 100) {
	return 0;
    }

    if (ainfo->q > 1 || ainfo->Q > 0) {
	ret = 1;
	if (arma_exact_ml(ainfo)) {
	    if (ainfo->P > 0) {
		ret = 0;
	    } else if (ainfo->p + ainfo->P > 0 && ainfo->nexo > 0) {
		ret = 0;
	    }
	}
    }

    return ret;
}

MODEL arma_model (const int *list, const double **Z, const DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn)
{
    double *coeff = NULL;
    int *alist = NULL;
    PRN *aprn = NULL;
    MODEL armod;
    struct arma_info ainfo;
    int init_done = 0;
    char flags = 0;
    int err = 0;

    if (!(opt & OPT_C)) {
	flags = ARMA_EXACT;
    }

    if (opt & OPT_V) {
	aprn = prn;
	errprn = prn;
    } else {
	errprn = NULL;
    }

    arma_info_init(&ainfo, flags, pdinfo);
    gretl_model_init(&armod); 
    gretl_model_smpl_init(&armod, pdinfo);

    alist = gretl_list_copy(list);
    if (alist == NULL) {
	armod.errcode = E_ALLOC;
	goto bailout;
    }

    err = arma_check_list(alist, opt, Z, pdinfo, &ainfo);
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
    err = user_arma_init(coeff, &ainfo, flags, &init_done, errprn);
    if (err) {
	armod.errcode = err;
	goto bailout;
    }

    /* second pass: try Hannan-Rissanen? */
    if (!init_done && prefer_hr_init(&ainfo)) {
	err = hr_init_check(pdinfo, &ainfo);
	if (!err) {
	    err = hr_arma_init(alist, coeff, Z, pdinfo, &ainfo, errprn);
#if ARMA_DEBUG
	    if (err) {
		fputs("hr_arma_init failed, will try ar_arma_init\n", stderr);
	    }
#endif
	}
	if (!err) {
	    init_done = 1;
	}
    }

    /* third pass: estimate pure AR model by OLS or NLS */
    if (!init_done) {
	err = ar_arma_init(alist, coeff, Z, pdinfo, &ainfo, errprn);
    }

    if (err) {
	armod.errcode = err;
	goto bailout;
    }

    if (flags & ARMA_EXACT) {
	kalman_arma(alist, coeff, Z, pdinfo, &ainfo, &armod, aprn);
    } else {
	bhhh_arma(alist, coeff, Z, pdinfo, &ainfo, &armod, aprn);
    }

 bailout:

    free(alist);
    free(coeff);
    free(ainfo.dy);

    /* cleanup in MA roots checker */
    bounds_checker_cleanup();

    errprn = NULL;

    return armod;
}
