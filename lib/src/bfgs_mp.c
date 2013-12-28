/* BFGS, and some related GMM utilities, in multiple precision */

#include <gmp.h>

static void free_triangular_mp_array (mpf_t **m, int n)
{
    if (m != NULL) {
	int i;

	for (i=0; i<n; i++) {
	    free(m[i]);
	}
	free(m);
    }
}

static mpf_t **triangular_mp_array_new (int n)
{
    mpf_t **m = malloc(n * sizeof *m);
    int i, j;

    if (m != NULL) {
	for (i=0; i<n; i++) {
	    m[i] = NULL;
	}
	for (i=0; i<n; i++) {
	    m[i] = malloc((i + 1) * sizeof **m);
	    if (m[i] == NULL) {
		free_triangular_mp_array(m, n);
		return NULL;
	    }
	}
 
	for (i=0; i<n; i++) {
	    for (j=0; j<i; j++) {
		mpf_init(m[i][j]);
	    }
	    mpf_init_set_d(m[i][i], 1.0);
	}
    }

    return m;
}

static void mp_workspace_free (mpf_t *m, int n)
{
    if (m != NULL) {
	int i;

	for (i=0; i<n; i++) {
	    mpf_clear(m[i]);
	}
	free(m);
    }
}

static mpf_t *mp_workspace_new (int n)
{
    mpf_t *m = malloc(n * sizeof *m);
    int i;

    if (m != NULL) {
	for (i=0; i<n; i++) {
	    mpf_init(m[i]);
	}
    }

    return m;
}

static int mp_richardson (double *b, mpf_t *g, int n,
			  BFGS_CRIT_FUNC func, void *data)
{
    mpf_t df[RSTEPS];
    mpf_t num, den, tmp, prod;
    mpf_t pf1, pf2;
    double eps = 1.0e-4;
    double d = 0.0001;
    double v = 2.0;
    double h, p4m;
    double bi0, f1, f2;
    int r = RSTEPS;
    int i, k, m;
    int err = 0;

    for (i=0; i<RSTEPS; i++) {
	mpf_init(df[i]);
    }

    mpf_init(num);
    mpf_init(den);
    mpf_init(tmp);
    mpf_init(prod);
    mpf_init(pf1);
    mpf_init(pf2);

    for (i=0; i<n; i++) {
	bi0 = b[i];
	h = d * b[i] + eps * (b[i] == 0.0);
	for (k=0; k<r; k++) {
	    b[i] = bi0 - h;
	    f1 = func(b, data);
	    b[i] = bi0 + h;
	    f2 = func(b, data);
	    if (na(f1) || na(f2)) {
		b[i] = bi0;
		err = 1;
		break;
	    }
	    mpf_set_d(pf1, f1);
	    mpf_set_d(pf2, f2);
	    mpf_sub(num, pf2, pf1);
	    mpf_set_d(den, 2.0 * h);
	    mpf_div(df[k], num, den);
	    h /= v;
	}
	if (err) {
	    break;
	}
	b[i] = bi0;
	p4m = 4.0;
	for (m=0; m<r-1; m++) {
	    for (k=0; k<r-m-1; k++) {
		mpf_set_d(tmp, p4m);
		mpf_mul(prod, df[k+1], tmp);
		mpf_sub(num, prod, df[k]);
		mpf_sub_ui(den, tmp, 1);
		mpf_div(df[k], num, den);
	    }
	    p4m *= 4.0;
	}
	mpf_set(g[i], df[0]);
    }

    for (i=0; i<RSTEPS; i++) {
	mpf_clear(df[i]);
    }

    mpf_clear(num);
    mpf_clear(den);
    mpf_clear(tmp);
    mpf_clear(prod);    
    mpf_clear(pf1);
    mpf_clear(pf2);    

    return err;
}

static int mp_simple_gradient (double *b, mpf_t *g, int n,
			       BFGS_CRIT_FUNC func, void *data,
			       int *redo)
{
    mpf_t num, den, pf1, pf2;
    double bi0, f1, f2;
    double h = 1.0e-8;
    int i, err = 0;

    mpf_init(num);
    mpf_init(den);
    mpf_init(pf1);
    mpf_init(pf2);

    for (i=0; i<n; i++) {
	bi0 = b[i];
	b[i] = bi0 - h;
	if (bi0 != 0.0 && fabs((bi0 - b[i]) / bi0) < B_RELMIN) {
	    fprintf(stderr, "mp_simple_gradient: switching to Richardson\n");
	    *redo = 1;
	    break;
	}
	f1 = func(b, data);
	b[i] = bi0 + h;
	f2 = func(b, data);
	b[i] = bi0;
	if (na(f1) || na(f2)) {
	    err = 1;
	    break;
	}
	mpf_set_d(pf1, f1);
	mpf_set_d(pf2, f2);
	mpf_sub(num, pf2, pf1);
	mpf_set_d(den, 2.0 * h);
	mpf_div(g[i], num, den);
    }

    mpf_clear(num);
    mpf_clear(den);
    mpf_clear(pf1);
    mpf_clear(pf2);

    return err;
}

static int mp_numeric_gradient (double *b, mpf_t *g, int n,
				BFGS_CRIT_FUNC func, void *data)
{
    int err = 0;

    if (libset_get_bool(BFGS_RSTEP)) {
	err = mp_richardson(b, g, n, func, data);
    } else {
	int redo = 0;

	err = mp_simple_gradient(b, g, n, func, data, &redo);
	if (redo) {
	    err = mp_richardson(b, g, n, func, data);
	}
    }

    return err;
}

static int mp_steplen (int n, int *pndelta, double *b, 
		       mpf_t *X, mpf_t *t, 
		       double *pf, BFGS_CRIT_FUNC cfunc, void *data, 
		       mpf_t g0, double f0, int *pfcount,
		       mpf_t *psteplen)
{
    mpf_t steplen, d, prod, tmp;
    double xi, f1 = *pf;
    int i, crit_ok = 0, fcount = 0;
    int ndelta;

    mpf_init_set_d(steplen, 1.0);
    mpf_init(d);
    mpf_init(prod);
    mpf_init(tmp);

    /* Below: iterate so long as (a) we haven't achieved an acceptable
       value of the criterion and (b) there is still some prospect
       of doing so.
    */    

    do {
	ndelta = n;
	crit_ok = 0;
	for (i=0; i<n; i++) {
	    mpf_mul(prod, steplen, t[i]);
	    mpf_add(tmp, X[i], prod);
	    b[i] = mpf_get_d(tmp);
	    xi = mpf_get_d(X[i]);
	    if (coeff_unchanged(b[i], xi)) {
		ndelta--;
	    }
	}
	if (ndelta > 0) {
	    f1 = cfunc(b, data);
	    mpf_set_d(tmp, acctol);
	    mpf_mul(prod, steplen, tmp);
	    mpf_mul(d, g0, prod);
	    mpf_neg(d, d);
	    fcount++;
#if 0
	    /* find the optimal steplength by quadratic interpolation; 
	       inspired by Kelley (1999), "Iterative Methods for Optimization", 
	       especially section 3.2.1. 
	    */
	    if (xna(f1)) {
		/* function goes into NA zone, presumably outside the 
		   admissible parameter space; hence, try a much smaller 
		   step. FIXME execution can come back here indefinitely.
		*/
		mpf_set_d(tmp, STEPFRAC);
		mpf_mul(steplen, steplen, tmp);
	    } else if (f1 < f0 + mpf_get_d(d)) {
		/* function computes, but goes down: try quadratic approx */
		/* steplen *= 0.5 * g0 / (f0 - f1 + g0); */
		mpf_t num, den;

		mpf_init(num);
		mpf_init(den);

		mpf_div_ui(num, g0, 2);
		mpf_set_d(tmp, f0 - f1);
		mpf_add(den, tmp, g0);
		mpf_div(tmp, num, den);
		mpf_mul(steplen, steplen, tmp);

		mpf_clear(num);
		mpf_clear(den);

		if (mpf_cmp_d(steplen, 1.0e-12) < 0) {
		    /* safeguarding */
		    fprintf(stderr, "safeguarding: f0 - f1 = %g, g0 = %g\n",
			    f0 - f1, mpf_get_d(g0));
		    mpf_set_d(steplen, 1.0e-12);
		    crit_ok = 1;
		}
	    } else {
		crit_ok = 1;
	    }
#else
	    /* find the optimal steplength by successive powers of STEPFRAC */
	    crit_ok = !na(f1) && (f1 >= f0 + mpf_get_d(d));
	    if (!crit_ok) {
		/* calculated criterion no good: try smaller step */
		mpf_set_d(tmp, STEPFRAC);
		mpf_mul(steplen, steplen, tmp);
	    }
#endif
	}
    } while (ndelta > 0 && !crit_ok);

    *pndelta = ndelta;
    *pfcount += fcount;
    *pf = f1;

    mpf_set(*psteplen, steplen);

    mpf_clear(steplen);
    mpf_clear(d);
    mpf_clear(prod);
    mpf_clear(tmp);

    return 0;
}

static void mp_H_to_I (mpf_t **H, int n)
{
    int i, j;

    for (i=0; i<n; i++) {
	for (j=0; j<i; j++) {
	    mpf_set_d(H[i][j], 0.0);
	}
	mpf_set_d(H[i][i], 1.0);
    }
}

static void mp_H_clear (mpf_t **H, int n)
{
    int i, j;

    for (i=0; i<n; i++) {
	for (j=0; j<i; j++) {
	    mpf_clear(H[i][j]);
	}
	mpf_clear(H[i][i]);
    }
}

static void mp_print_iter_info (int iter, double crit, int type, int k, 
				const double *b, mpf_t *mpg, 
				mpf_t *sl, PRN *prn)
{
    double g[k];
    int i;

    for (i=0; i<k; i++) {
	g[i] = mpf_get_d(mpg[i]);
    }

    print_iter_info(iter, crit, type, k, b, g, mpf_get_d(*sl), prn);
}

int mp_BFGS (double *b, int n, int maxit, double reltol,
	     int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	     int crittype, void *data, gretlopt opt, PRN *prn)
{
    int done;
    mpf_t *wspace = NULL;
    mpf_t **H = NULL;
    mpf_t *g, *t, *X, *c;
    int verbskip, verbose = (opt & OPT_V);
    int fcount, gcount, ndelta = 0;
    int show_activity = 0;
    mpf_t sumgrad, gradnorm;
    mpf_t s, steplen;
    mpf_t prod, tmp, D1, D2;
    double f, f0, fmax, gradmax;
    int i, j, ilast, iter;
    int err = 0;

    mpf_init(sumgrad);
    mpf_init(gradnorm);
    mpf_init(s);
    mpf_init(steplen);
    mpf_init(prod);
    mpf_init(tmp);
    mpf_init(D1);
    mpf_init(D2);

    optim_get_user_values(b, n, &maxit, &reltol, &gradmax, opt, prn);

    wspace = mp_workspace_new(4 * n);
    H = triangular_mp_array_new(n);

    if (wspace == NULL || H == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    g = wspace;
    t = g + n;
    X = t + n;
    c = X + n;

    f = cfunc(b, data);

    if (na(f)) {
	gretl_errmsg_set(_("BFGS: initial value of objective function is not finite"));
	err = E_NAN;
	goto bailout;
    }

    f0 = fmax = f;
    iter = ilast = fcount = gcount = 1;

    mp_numeric_gradient(b, g, n, cfunc, data);

    if (maxit == 0) {
	goto skipcalc;
    }

    verbskip = libset_get_int("bfgs_verbskip");
    show_activity = show_activity_func_installed();

    do {
	if (bfgs_print_iter(verbose, verbskip, iter)) {
	    mp_print_iter_info(iter, f, crittype, n, b, g, &steplen, prn);
	}

	if (show_activity && (iter % 10 == 0)) {
	    show_activity_callback();
	}

	if (iter > 1 && ilast == gcount) {
	     /* restart: set curvature matrix to I */
	    mp_H_to_I(H, n);
	}

	for (i=0; i<n; i++) {
	    /* copy coefficients to X, gradient to c */
	    mpf_set_d(X[i], b[i]);
	    mpf_set(c[i], g[i]);
	}

	mpf_set_d(gradnorm, 0.0);
	mpf_set_d(sumgrad, 0.0);

	for (i=0; i<n; i++) {
	    mpf_set_d(s, 0.0);
	    for (j=0; j<=i; j++) {
		mpf_mul(prod, H[i][j], g[j]);
		mpf_add(s, s, prod);
	    }
	    for (j=i+1; j<n; j++) {
		mpf_mul(prod, H[j][i], g[j]);
		mpf_add(s, s, prod);	
	    }
	    mpf_set(t[i], s);
	    mpf_mul(prod, s, g[i]);
	    mpf_add(sumgrad, sumgrad, prod);
	    mpf_set_d(tmp, b[i]);
	    mpf_mul(prod, tmp, g[i]);
	    mpf_abs(prod, prod);
	    mpf_add(gradnorm, gradnorm, prod);
	}

	mpf_div_ui(gradnorm, gradnorm, (unsigned long) n);
	mpf_sqrt(gradnorm, gradnorm);

	if (mpf_sgn(sumgrad) > 0) { 
	    /* heading in the right direction */
	    mp_steplen(n, &ndelta, b, X, t, &f, cfunc, data, 
		       sumgrad, fmax, &fcount, &steplen);

	    done = fabs(fmax - f) <= reltol * (fabs(fmax) + reltol);

	    /* prepare to stop if relative change is small enough */
	    if (done) {
		ndelta = 0;
		fmax = f;
	    }

	    if (ndelta > 0) {
		/* making progress */
		fmax = f;
		mp_numeric_gradient(b, g, n, cfunc, data);
		gcount++;
		iter++;
		mpf_set_d(D1, 0.0);
		for (i=0; i<n; i++) {
		    mpf_mul(t[i], t[i], steplen);
		    mpf_sub(c[i], c[i], g[i]);
		    mpf_mul(prod, t[i], c[i]);
		    mpf_add(D1, D1, prod);
		}
		if (mpf_sgn(D1) > 0) {
		    mpf_set_d(D2, 0.0);
		    for (i=0; i<n; i++) {
			mpf_set_d(s, 0.0);
			for (j=0; j<=i; j++) {
			    mpf_mul(prod, H[i][j], c[j]);
			    mpf_add(s, s, prod);
			}
			for (j=i+1; j<n; j++) {
			    mpf_mul(prod, H[j][i], c[j]);
			    mpf_add(s, s, prod);
			}
			mpf_set(X[i], s);
			mpf_mul(prod, s, c[i]);
			mpf_add(D2, D2, prod);
		    }
		    mpf_div(tmp, D2, D1);
		    mpf_add_ui(D2, tmp, 1);
		    for (i=0; i<n; i++) {
			for (j=0; j<=i; j++) {
			    mpf_set(tmp, D2);
			    mpf_mul(prod, t[i], t[j]);
			    mpf_mul(tmp, tmp, prod);
			    mpf_mul(prod, X[i], t[j]);
			    mpf_sub(tmp, tmp, prod);
			    mpf_mul(prod, t[i], X[j]);
			    mpf_sub(tmp, tmp, prod);
			    mpf_div(tmp, tmp, D1);
			    mpf_add(H[i][j], H[i][j], tmp);
			}
		    }
		} else {
		    /* D1 <= 0.0 */
		    ilast = gcount;
		}
	    } else if (ilast < gcount) {
		ndelta = n;
		ilast = gcount;
	    }
	} else if (mpf_sgn(sumgrad) == 0) {
	    fprintf(stderr, "gradient is exactly zero!\n");
	    break;
	} else {
	    /* heading in the wrong direction */
	    if (ilast == gcount) {
		/* we just did a reset, so don't reset again; instead set 
		   ndelta = 0 so that we exit the main loop
		*/
		ndelta = 0;
		if (gcount == 1) {
		    err = E_NOCONV;
		}
	    } else {
		/* reset for another attempt */
		ilast = gcount;
		ndelta = n;
	    }
	}

	if (iter >= maxit) {
	    break;
	}

	if (gcount - ilast > 2 * n) {
	    /* periodic restart of curvature computation */
	    ilast = gcount;
	}

    } while (ndelta > 0 || ilast < gcount);

    if (iter >= maxit) {
	fprintf(stderr, _("stopped after %d iterations\n"), iter);
	err = E_NOCONV;
    } else if (mpf_cmp_d(gradnorm, gradmax) > 0) {
	err = E_NOCONV;
    } else if (fmax < f0) {
	/* allow a small sloppiness factor here? */
	double rdiff;

	rdiff = (f0 == 0.0)? -fmax : fabs((f0 - fmax) / f0);
	if (rdiff > 1.0e-12) {
	    fprintf(stderr, "failed to match initial value of objective function:\n"
		    "   f0=%.18g\n fmax=%.18g\n", f0, fmax);
	    err = E_NOCONV;
	}
    }

    if (!err && mpf_cmp_d(gradnorm, GRAD_TOLER) > 0) {
	gretl_warnmsg_sprintf(_("norm of gradient = %g"), mpf_get_d(gradnorm));
	set_gretl_warning(W_GRADIENT);
    }

 skipcalc:

    *fncount = fcount;
    *grcount = gcount;

    if (verbose) {
	mp_print_iter_info(-1, f, crittype, n, b, g, &steplen, prn);
	pputc(prn, '\n');
    }

 bailout:

    mpf_clear(sumgrad);
    mpf_clear(gradnorm);
    mpf_clear(s);
    mpf_clear(steplen);
    mpf_clear(prod);
    mpf_clear(tmp);
    mpf_clear(D1);
    mpf_clear(D2);

    mp_workspace_free(wspace, 4*n);

    if (H != NULL) {
	mp_H_clear(H, n);
	free_triangular_mp_array(H, n);
    }

    return err;
}

typedef struct gretl_mp_matrix_ gretl_mp_matrix;

struct gretl_mp_matrix_ {
    int rows;
    int cols;
    mpf_t *val;
};

static gretl_mp_matrix *gretl_mp_matrix_alloc (int r, int c)
{
    gretl_mp_matrix *m = malloc(sizeof *m);

    if (m != NULL) {
	int i, n = r * c;

	m->val = malloc(n * sizeof *m->val);
	if (m->val == NULL) {
	    free(m);
	    m = NULL;
	} else {
	    m->rows = r;
	    m->cols = c;
	    for (i=0; i<n; i++) {
		mpf_init(m->val[i]);
	    }
	}
    }

    return m;
}

static void gretl_mp_matrix_free (gretl_mp_matrix *m)
{
    if (m != NULL) {
	int i, n = m->rows * m->cols;

	for (i=0; i<n; i++) {
	    mpf_clear(m->val[i]);
	}
	free(m->val);
	free(m);
    }
}

static void gretl_mp_matrix_set (gretl_mp_matrix *m,
				 int i, int j,
				 mpf_t *val)
{
    int k = j * m->rows + i;

    mpf_set(m->val[k], *val);
}

static void gretl_mp_matrix_get (gretl_mp_matrix *m,
				 int i, int j,
				 mpf_t *val)
{
    int k = j * m->rows + i;

    mpf_set(*val, m->val[k]);
}

static int mp_columnwise_product (const gretl_matrix *A,
				  const gretl_matrix *B,
				  const gretl_matrix *S,
				  gretl_mp_matrix *C)
{
    int k = A->cols;
    int n = B->cols;
    int T = A->rows;
    mpf_t mpx, mpy, prod;
    int i, j, t, p;

    mpf_init(mpx);
    mpf_init(mpy);
    mpf_init(prod);

    p = 0;
    for (i=0; i<k; i++) {
	for (j=0; j<n; j++) {
	    if (S == NULL || gretl_matrix_get(S, i, j) != 0) {
		for (t=0; t<T; t++) {
		    mpf_set_d(mpx, gretl_matrix_get(A, t, i));
		    mpf_set_d(mpy, gretl_matrix_get(B, t, j));
		    mpf_mul(prod, mpx, mpy);
		    gretl_mp_matrix_set(C, t, p, &prod);
		}
		p++;
	    }
	}
    }

    mpf_clear(mpx);
    mpf_clear(mpy);
    mpf_clear(prod);    
	    
    return 0;
}

double mp_gmm_criterion (const gretl_matrix *A,
			 const gretl_matrix *B,
			 const gretl_matrix *S,
			 const gretl_matrix *W,
			 gretl_matrix *dsum,
			 gretl_matrix *dtmp,
			 int noc, int *err)
{
    static gretl_mp_matrix *sum;
    static gretl_mp_matrix *C;
    double crit = NADBL;
    mpf_t tmp, cti, mpc;
    mpf_t prod, wp;
    int T, k = noc;
    int i, j, p, t;

    if (A == NULL) {
	/* clean-up signal */
	gretl_mp_matrix_free(sum);
	gretl_mp_matrix_free(C);
	sum = NULL;
	C = NULL;
	return 0;
    }

    T = A->rows;

    if (sum == NULL) {
	/* starting from scratch */
	sum = gretl_mp_matrix_alloc(k, 1);
	if (sum == NULL) {
	    *err = E_ALLOC;
	    return NADBL;
	}
	C = gretl_mp_matrix_alloc(T, k);
	if (C == NULL) {
	    gretl_mp_matrix_free(sum);
	    *err = E_ALLOC;
	    return NADBL;
	}
    }  

    mpf_init(tmp);
    mpf_init(cti);
    mpf_init(mpc);
    mpf_init(prod);
    mpf_init(wp);

    mp_columnwise_product(A, B, S, C);

    for (i=0; i<T*k; i++) {
	dtmp->val[i] = mpf_get_d(C->val[i]);
    }

    for (i=0; i<k; i++) {
	mpf_set_d(tmp, 0.0);
	for (t=0; t<T; t++) {
	    gretl_mp_matrix_get(C, t, i, &cti);
	    mpf_add(tmp, tmp, cti);
	}
	mpf_set(sum->val[i], tmp);
	dsum->val[i] = mpf_get_d(tmp);
    }

    p = 0;
    for (j=0; j<k; j++) {
	mpf_set_d(tmp, 0.0);
	for (i=0; i<k; i++) {
	    mpf_set_d(wp, W->val[p++]);
	    mpf_mul(prod, sum->val[i], wp);
	    mpf_add(tmp, tmp, prod);
	}
	mpf_mul(prod, tmp, sum->val[j]); 
	mpf_add(mpc, mpc, prod);
    }

    mpf_neg(mpc, mpc);
    crit = mpf_get_d(mpc);

    mpf_clear(tmp);
    mpf_clear(cti);
    mpf_clear(mpc);
    mpf_clear(prod);
    mpf_clear(wp);

    return crit;
}

void mp_gmm_cleanup (void)
{
    mp_gmm_criterion(NULL, NULL, NULL, NULL,
		     NULL, NULL, 0, NULL);
}
