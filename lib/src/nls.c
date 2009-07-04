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

/* Nonlinear least squares for libgretl, using minpack; also Maximum
   Likelihood and GMM estimation using BFGS.  However, much of the
   GMM-specific material is in the companion file, gmm.c
*/

#include "libgretl.h" 
#include "libset.h"
#include "usermat.h"
#include "matrix_extra.h"
#include "gretl_func.h"
#include "nlspec.h"
#include "cmd_private.h"
#include "gretl_scalar.h"
#include "gretl_bfgs.h"
#include "tsls.h"

#include "gretl_f2c.h"
#include "../../minpack/minpack.h"  

#define NLS_DEBUG 0
#define ML_DEBUG 0

enum {
    ANALYTIC_DERIVS = 1 << 0,
    NLS_AUTOREG     = 1 << 1
} nls_flags;

struct parm_ {
    char name[VNAMELEN];  /* name of parameter */
    int type;             /* type of parameter (scalar or vector) */
    char *deriv;          /* string representation of derivative of regression
			     function with respect to param (or NULL) */
    int dnum;             /* ID number of variable holding the derivative */
    int nc;               /* number of individual coefficients associated
			     with the parameter */
    char dname[VNAMELEN]; /* name of variable holding the derivative */
    gretl_matrix *vec;    /* pointer to vector parameter */
    gretl_matrix *dmat;   /* pointer to matrix derivative */
};

#define scalar_param(s,i) (s->params[i].type == GRETL_TYPE_DOUBLE)
#define matrix_deriv(s,i) (s->params[i].dmat != NULL)
#define scalar_deriv(s,i) (gretl_is_scalar(s->params[i].dname))

#define numeric_mode(s) (!(s->flags & ANALYTIC_DERIVS))
#define analytic_mode(s) (s->flags & ANALYTIC_DERIVS)

/* file-scope global variables */

static nlspec private_spec;
static integer one = 1;

static char *adjust_saved_nlfunc (char *s);

static void set_numeric_mode (nlspec *s)
{
    s->flags &= ~ANALYTIC_DERIVS;
}

static void set_analytic_mode (nlspec *s)
{
    s->flags |= ANALYTIC_DERIVS;
}

static void destroy_genrs_array (GENERATOR **genrs, int n)
{
    int i;

    for (i=0; i<n; i++) {
	destroy_genr(genrs[i]);
    }

    free(genrs);
}

static int check_lhs_vec (nlspec *s)
{
    int v = gretl_vector_get_length(s->lvec);

    if (v != s->nobs) {
	fprintf(stderr, "LHS vector should be of length %d, is %d\n", 
		s->nobs, v);
	return 1;
    }

    return 0;
}

static int check_derivative_matrix (int i, gretl_matrix *m,
				    nlspec *s)
{
    int r, c, v;

    if (m == NULL) {
	fprintf(stderr, "param %d, got NULL matrix derivative\n", i);
	return 1;
    }

    r = gretl_matrix_rows(m);
    c = gretl_matrix_cols(m);
    v = s->params[i].nc;
    
    if (c != v || (r != 1 && r != s->nobs)) {
	fprintf(stderr, "matrix deriv for param %d is %d x %d: WRONG\n", 
		i, r, c);
	fprintf(stderr, " should be %d x %d, or %d x %d\n", s->nobs,
		v, 1, v);
	return 1;
    }

    return 0;
}

static int nls_dynamic_check (nlspec *s, char *formula)
{
    GENERATOR *genr;
    double *y;
    int t, err = 0;

    /* back up the dependent variable */
    y = copyvec((*s->Z)[s->dv], s->dinfo->n);
    if (y == NULL) {
	return E_ALLOC;
    }

    /* compile the formula for the dependent variable and see
       if it is autoregressive */
    strcpy(formula, s->nlfunc);
    adjust_saved_nlfunc(formula);

    genr = genr_compile(formula, s->Z, s->dinfo, OPT_P, &err);

    if (!err && genr_is_autoregressive(genr)) {
	s->flags |= NLS_AUTOREG;
    }

    /* restore the dependent variable */
    for (t=0; t<s->dinfo->n; t++) {
	(*s->Z)[s->dv][t] = y[t];
    }    

    destroy_genr(genr);
    free(y);

    return err;
}

/* we "compile" the required equations first, so we can subsequently
   execute the compiled versions for maximum efficiency 
*/

static int nls_genr_setup (nlspec *s)
{
    GENERATOR **genrs;
    gretl_matrix *m;
    char formula[MAXLINE];
    int i, j, ngen, np;
    int v, err = 0;

    s->ngenrs = 0;
    np = (analytic_mode(s))? s->nparam : 0;
    ngen = s->naux + np;
    if (s->nlfunc != NULL) {
	ngen++;
    }

#if NLS_DEBUG
    fprintf(stderr, "nls_genr_setup: current v = %d, n_gen = %d\n", 
	    s->dinfo->v, ngen);
#endif

    genrs = malloc(ngen * sizeof *genrs);
    if (genrs == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<ngen; i++) {
	genrs[i] = NULL;
    }

    j = 0;

    for (i=0; i<ngen && !err; i++) {
	char *dname = NULL;

	if (i < s->naux) {
	    /* auxiliary variables */
	    strcpy(formula, s->aux[i]);
	} else if (i == s->naux) {
	    /* residual or likelihood function */
	    if (*s->lhname != '\0') {
		sprintf(formula, "%s = %s", s->lhname, s->nlfunc);
	    } else {
		if (s->ci == NLS && dataset_is_time_series(s->dinfo)) {
		    err = nls_dynamic_check(s, formula);
		}
		sprintf(formula, "$nl_y = %s", s->nlfunc);
	    }
	} else {
	    /* derivatives/gradients */
	    sprintf(s->params[j].dname, "$nl_x%d", i);
	    if (scalar_param(s, j)) {
		sprintf(formula, "%s = %s", s->params[j].dname, 
			s->params[j].deriv);
	    } else {
		sprintf(formula, "matrix %s = %s",  s->params[j].dname,
			s->params[j].deriv);
	    }
	    dname = s->params[j].dname;
	    j++;
	}

	if (!err) {
	    genrs[i] = genr_compile(formula, s->Z, s->dinfo, OPT_P, &err);
	}

	if (err) {
	    fprintf(stderr, "genr_compile: genrs[%d] = %p, err = %d\n", i, 
		    (void *) genrs[i], err);
	    fprintf(stderr, "formula: '%s'\n", formula);
	    fprintf(stderr, "%s\n", gretl_errmsg_get());
	    break;
	}

	/* see if the formula actually works, and flush out NAs
	   while we're at it
	*/
	genr_set_na_check(genrs[i]);
	err = execute_genr(genrs[i], s->Z, s->dinfo, OPT_S, s->prn);
	genr_unset_na_check(genrs[i]);

	if (err) {
	    fprintf(stderr, "execute_genr: genrs[%d] = %p, err = %d\n", i, 
		    (void *) genrs[i], err);
	    fprintf(stderr, "formula: '%s'\n", formula);
	    break;
	}

	v = genr_get_output_varnum(genrs[i]);
	m = genr_get_output_matrix(genrs[i]);

	if (v == 0 && m == NULL) {
	    /* not a series, not a matrix: should be scalar */
	    if (genr_is_print(genrs[i])) {
		continue;
	    } else if (i >= s->naux && (dname == NULL || !gretl_is_scalar(dname))) {
		fprintf(stderr, "nls_genr_setup: bad type: %s\n", formula);
		err = E_TYPES;
		break;
	    }
	}

	if (i == s->naux) {
	    s->lhv = v;
	    s->lvec = m;
	    if (s->lvec != NULL) {
		err = check_lhs_vec(s);
	    }
	} else if (j > 0) {
	    int k = j - 1;

	    s->params[k].dnum = v;
	    s->params[k].dmat = m;

	    if (m != NULL || !scalar_param(s, k)) {
		err = check_derivative_matrix(k, m, s);
	    } 
	}

#if NLS_DEBUG
	fprintf(stderr, " formula '%s'\n", formula);
	fprintf(stderr, " v = %d, m = %p\n", v, (void *) m);
	if (v > 0) {
	    fprintf(stderr, " first value: Z[%d][%d] = %g\n", 
		    v, s->t1, (*s->Z)[v][s->t1]);
	}
#endif
    }

    if (err) {
	destroy_genrs_array(genrs, ngen);
    } else {
	s->ngenrs = ngen;
	s->genrs = genrs;
    }

    return err;
}

/* If i == 0 we're calculating the function; if i > 0 we're calculating
   a derivative.  Either way, we recalculate any auxiliary variables
   first.
*/

static int nls_auto_genr (nlspec *s, int i)
{
    int j;

    s->generr = 0;

#if NLS_DEBUG
    fprintf(stderr, "nls_auto_genr: input i = %d\n", i);
#endif

    if (s->genrs == NULL) {
	s->generr = nls_genr_setup(s);
	if (s->generr) {
	    fprintf(stderr, " nls_genr_setup failed\n");
	}
	return s->generr;
    }

    for (j=0; j<s->naux; j++) {
#if NLS_DEBUG
	fprintf(stderr, " generating aux var %d:\n %s\n", j, s->aux[j]);
#endif
	s->generr = execute_genr(s->genrs[j], s->Z, s->dinfo, OPT_S, s->prn);
	if (s->generr) {
	    return s->generr;
	}
    }

    if (i == 0 && s->nlfunc == NULL) {
	/* we're done */
	return s->generr;
    }

    j = s->naux + i;
#if NLS_DEBUG
    fprintf(stderr, " j = naux+i = %d+%d = %d: executing genr[%d] at %p:\n", 
	    s->naux, i, j, j, (void *) s->genrs[j]);
    fprintf(stderr, " Z[%d] = %p\n", s->dinfo->v-1, (void *) (*s->Z)[s->dinfo->v-1]);
#endif
    s->generr = execute_genr(s->genrs[j], s->Z, s->dinfo, OPT_S, s->prn);

    /* make sure we have a correct pointer to matrix deriv */
    if (!s->generr && i > 0 && matrix_deriv(s, i-1)) {
	gretl_matrix *m = genr_get_output_matrix(s->genrs[j]);

	s->generr = check_derivative_matrix(i-1, m, s);
    }

#if NLS_DEBUG
    if (s->generr) {
	int v = genr_get_output_varnum(s->genrs[j]);

	fprintf(stderr, " varnum = %d, err = %d\n", v, s->generr);
	errmsg(s->generr, s->prn);
    } 
#endif

    return s->generr;    
}

/* wrappers for the above to enhance comprehensibility below */

int nl_calculate_fvec (nlspec *s)
{
    return nls_auto_genr(s, 0);
}

static int nls_calculate_deriv (nlspec *s, int i)
{
    return nls_auto_genr(s, i + 1);
}

/* end wrappers */

void nlspec_destroy_arrays (nlspec *s)
{
    free(s->params);
    s->params = NULL;
    s->nparam = 0;

    free(s->coeff);
    s->coeff = NULL;
    s->ncoeff = 0;
}

/* add a vector of coefficients to the model specification */

static int push_vec_coeffs (nlspec *s, gretl_matrix *m, int k)
{
    double *coeff;
    int i, nc = s->ncoeff;

    if (k == 0) {
	return E_DATA;
    }
    
    coeff = realloc(s->coeff, (nc + k) * sizeof *coeff);
    if (coeff == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<k; i++) {
	coeff[nc + i] = m->val[i];
    }

#if NLS_DEBUG
    fprintf(stderr, "added %d coeffs\n", k);
#endif

    s->coeff = coeff;
    s->ncoeff = nc + k;

    return 0;
}

/* add a scalar coefficient to the model specification */

static int push_scalar_coeff (nlspec *s, double x)
{
    double *coeff;
    int nc = s->ncoeff;
    
    coeff = realloc(s->coeff, (nc + 1) * sizeof *coeff);
    if (coeff == NULL) {
	return E_ALLOC;
    }

    coeff[nc] = x;

#if NLS_DEBUG
    fprintf(stderr, "added coeff[%d] = %g\n", nc, x);
#endif

    s->coeff = coeff;
    s->ncoeff = nc + 1;

    return 0;
}

/* add a parameter to the model specification: this may be
   either a scalar or a vector
*/

static int nlspec_push_param (nlspec *s, const char *name, char *deriv)
{
    parm *params, *p;
    int np = s->nparam;
    int err;
    
    params = realloc(s->params, (np + 1) * sizeof *params);
    if (params == NULL) {
	return E_ALLOC;
    }

    p = &params[np];

    p->name[0] = '\0';
    strncat(p->name, name, VNAMELEN - 1);
    p->type = GRETL_TYPE_DOUBLE;
    p->deriv = deriv;
    p->dnum = 0;
    p->nc = 1;
    p->dname[0] = '\0';
    p->vec = NULL;
    p->dmat = NULL;

#if NLS_DEBUG
    fprintf(stderr, "added param[%d] = '%s'\n", np, p->name);
#endif

    s->params = params;
    s->nparam = np + 1;

    if (gretl_is_scalar(name)) {
	err = push_scalar_coeff(s, gretl_scalar_get_value(name));
    } else {
	gretl_matrix *m = get_matrix_by_name(name);
	int k = gretl_vector_get_length(m);

#if NLS_DEBUG
	fprintf(stderr, "vector param: m = %p, k = %d\n", (void *) m, k);
#endif

	p->type = GRETL_TYPE_MATRIX;
	p->vec = m;
	p->nc = k;
	err = push_vec_coeffs(s, m, k);
	if (!err) {
	    s->nvec += 1;
	}
    }

    return err;
}

/* Do we have a valid name (the name of a scalar or vector?):
   if so return 0, else return E_DATATYPE.
*/

static int check_param_name (const char *name)
{
    if (gretl_is_scalar(name)) {
	return 0;
    } else {
	gretl_matrix *m = get_matrix_by_name(name);

	if (gretl_vector_get_length(m) == 0) {
	    return E_DATATYPE;
	} else {
	    return 0;
	}
    }
}

/* For case where analytical derivatives are not given, the user
   must supply a line like:

   params b0 b1 b2 b3

   specifying the parameters to be estimated.  Here we parse such a
   list and add the parameter info to the spec.  The terms in the list
   must be pre-existing scalar variables or vectors.
*/

static int 
nlspec_add_params_from_line (nlspec *s, const char *str,
			     const double **Z, const DATAINFO *pdinfo)
{
    int i, nf = count_fields(str);
    int err = 0;

    if (s->params != NULL || nf == 0) {
	return E_DATA;
    }

#if NLS_DEBUG
    fprintf(stderr, "nlspec_add_params_from_line:\n "
	    "line = '%s', nf = %d\n", str, nf);
#endif

    for (i=0; i<nf && !err; i++) {
	char *name = gretl_word_strdup(str, &str);

	if (name == NULL) {
	    err = E_ALLOC;
	} else {
	    err = check_param_name(name);
	}

	if (!err) {
	    err = nlspec_push_param(s, name, NULL);
	}

	free(name);
    }

    if (err) {
	nlspec_destroy_arrays(s);
    } 

    return err;
}

/**
 * nlspec_add_param_list:
 * @spec: nls specification.
 * @np: number of parameters.
 * @vals: array of initial parameter values.
 * @names: array of parameter names.
 *
 * Adds to @spec a list of (scalar) parameters to be estimated.
 * For an example of use see arma.c in the gretl plugin directory.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int nlspec_add_param_list (nlspec *spec, int np, double *vals,
			   char **names, double ***pZ,
			   DATAINFO *pdinfo)
{
    int i, err = 0;

    if (spec->params != NULL || np == 0) {
	return E_DATA;
    }

    for (i=0; i<np && !err; i++) {
	err = gretl_scalar_add(names[i], vals[i]);
	if (!err) {
	    err = nlspec_push_param(spec, names[i], NULL);
	}
    }

    if (err) {
	nlspec_destroy_arrays(spec);
    } 

    return err;
}

/* update the 'external' values of scalars or matrices using
   the values produced by the optimizer.
*/

int update_coeff_values (const double *b, nlspec *s)
{
    int i, j, k = 0;
    parm *p;

    for (i=0; i<s->nparam; i++) {
	p = &s->params[i];
	if (scalar_param(s, i)) {
	    gretl_scalar_set_value(p->name, b[k++]);
	} else {
	    gretl_matrix *m = get_matrix_by_name(p->name);

	    if (m == NULL) {
		fprintf(stderr, "Couldn't find location for coeff %d\n", k);
		return E_DATA;
	    } else {
		if (m != p->vec) {
		    fprintf(stderr, "*** coeff_address: by name, '%s' is at %p; "
			    "stored addr = %p\n", p->name,
			    (void *) m, (void *) p->vec);
		    p->vec = m;
		}
		for (j=0; j<p->nc; j++) {
		    gretl_vector_set(m, j, b[k++]);
		}
	    }
	}
    }

    return 0;
}

#define NLS_SKIP_MISSING 0 /* not properly tested yet */

#if NLS_SKIP_MISSING

static int nls_make_trimmed_dataset (nlspec *spec, int t1, int t2)
{
    DATAINFO *dinfo = NULL;
    double ***pZ = NULL;
    double **origZ = *spec->Z;
    int nvar = spec->dinfo->v;
    int nobs = 0;
    int i, t, s;

    spec->missmask = malloc(spec->dinfo->n + 1);
    if (spec->missmask == NULL) {
	return E_ALLOC;
    }

    spec->missmask[spec->dinfo->n] = '\0';
    memset(spec->missmask, '0', spec->dinfo->n);

    for (t=t1; t<=t2; t++) {
	if (na(origZ[spec->lhv][t])) {
	    spec->missmask[t] = '1';
	} else {
	    nobs++;
	}
    }

    if (nobs < spec->ncoeff) {
	return E_DF;
    }

    pZ = malloc(sizeof *pZ);
    if (pZ == NULL) {
	return E_ALLOC;
    }

    dinfo = create_auxiliary_dataset(pZ, nvar, nobs);
    if (dinfo == NULL) {
	return E_ALLOC;
    }

    for (i=1; i<nvar; i++) {
	strcpy(dinfo->varname[i], spec->dinfo->varname[i]);
	s = 0;
	for (t=t1; t<=t2; t++) {
	    if (!na(origZ[spec->lhv][t])) {
		(*pZ)[i][s++] = origZ[i][t];
	    }
	}
    }

    spec->real_t1 = t1;
    spec->real_t2 = t2;

    spec->Z = pZ;
    spec->dinfo = dinfo;
    spec->t1 = 0;
    spec->t2 = nobs - 1;
    spec->nobs = nobs;

#if NLS_DEBUG
    fprintf(stderr, "s->t1 = %d, s->t2 = %d, s->nobs = %d, nvar = %d\n",
	    spec->t1, spec->t2, spec->nobs, nvar);
#endif

    return 0;
}

#endif

/* Adjust starting and ending points of sample if need be, to avoid
   missing values; abort if there are missing values within the
   (possibly reduced) sample range.  For this purpose we generate the
   nls residual variable.  
*/

static int nl_missval_check (nlspec *s)
{
    int t1 = s->t1, t2 = s->t2;
    int t, v;
    int err = 0;

#if NLS_DEBUG
    fprintf(stderr, "nl_missval_check: calling nl_calculate_fvec\n");
#endif

    /* calculate the function (NLS residual, MLE likelihood) */
    err = nl_calculate_fvec(s);
    if (err) {
	return err;
    }

    if (s->lvec != NULL) {
	/* vector result */
	goto nl_miss_exit;
    }

    /* ID number of LHS variable */
    v = s->lhv;

#if NLS_DEBUG
    fprintf(stderr, " checking var %d (%s)\n",
	    v, s->dinfo->varname[v]);
    fprintf(stderr, "  before trimming: spec->t1 = %d, spec->t2 = %d\n",
	    s->t1, s->t2);
#endif

    for (t1=s->t1; t1<=s->t2; t1++) {
	if (!na((*s->Z)[v][t1])) {
	    break;
	}
    }

    for (t2=s->t2; t2>=t1; t2--) {
	if (!na((*s->Z)[v][t2])) {
	    break;
	}
    }

    if (t2 - t1 + 1 == 0) {
	fprintf(stderr, "nl_missval_check: no valid data\n");
	return E_DATA;
    }

    if (t2 - t1 + 1 < s->ncoeff) {
	return E_DF;
    }

    for (t=t1; t<=t2; t++) {
	if (na((*s->Z)[v][t])) {
	    fprintf(stderr, "  after setting t1=%d, t2=%d, "
		    "got NA for var %d at obs %d\n", t1, t2, v, t);
#if NLS_SKIP_MISSING
	    return nls_make_trimmed_dataset(s, t1, t2);
#else
	    return E_MISSDATA;
#endif
	}
    }  

 nl_miss_exit:

    s->t1 = s->real_t1 = t1;
    s->t2 = s->real_t2 = t2;
    s->nobs = t2 - t1 + 1;

#if NLS_DEBUG
    fprintf(stderr, "  after: spec->t1 = %d, spec->t2 = %d, spec->nobs = %d\n\n",
	    s->t1, s->t2, s->nobs);
#endif

    return err;
}

/* the next two functions are used in the context of BFGS */

static double get_mle_ll (const double *b, void *p)
{
    nlspec *s = (nlspec *) p;
    double x;
    int t, k;

    update_coeff_values(b, s);

    if (nl_calculate_fvec(s)) {
	return NADBL;
    }

    s->crit = 0.0;

    if (s->lvec != NULL) {
	k = gretl_vector_get_length(s->lvec);
	for (t=0; t<k; t++) {
	    x = s->lvec->val[t];
	    if (na(x)) {
		return NADBL;
	    }
	    s->crit += x;
	}
    } else {
	k = s->lhv;
	for (t=s->t1; t<=s->t2; t++) {
	    x = (*s->Z)[k][t];
	    if (na(x)) {
		return NADBL;
	    }
	    s->crit += x;
	}
    }

    return s->crit;
}

/* this function is used in the context of the minpack callback, and
   also for checking derivatives in the MLE case
*/

static int nl_function_calc (double *f, void *p)
{
    nlspec *s = (nlspec *) p;
    const double *y;
    int t, err = 0;

#if NLS_DEBUG > 1
    fprintf(stderr, "\n*** nl_function_calc called\n");
#endif

    /* calculate residual given current parameter estimates */
    if ((err = nl_calculate_fvec(s))) {
	return err;
    }

    s->crit = 0.0;

    if (s->lvec != NULL) {
	y = s->lvec->val;
    } else {
	y = (*s->Z)[s->lhv] + s->t1;
    }

    /* transcribe from dataset to array f */

    for (t=0; t<s->nobs; t++) {
	if (na(y[t])) {
	    fprintf(stderr, "nl_calculate_fvec: produced NA at obs %d\n", t);
	    return 1;
	}

	f[t] = y[t];

#if NLS_DEBUG > 1
	fprintf(stderr, "fvec[%d] = %.14g\n", t,  f[t]);
#endif

	if (s->ci == MLE) {
	    s->crit += f[t];
	} else {
	    s->crit += f[t] * f[t];
	}
    }

    s->iters += 1;

    if (s->ci == NLS && (s->opt & OPT_V)) {
	pprintf(s->prn, _("iteration %2d: SSR = %.8g\n"), s->iters, s->crit);
    }

    return 0;
}

static int get_nls_derivs (int k, int T, int offset, double *g, double **G,
			   void *p)
{
    nlspec *spec = (nlspec *) p;
    double *gi;
    double x;
    int j, t;
    int err = 0;

    if (g != NULL) {
	gi = g;
    } else if (G != NULL) {
	gi = (*G) + offset;
    } else {
	return 1;
    }

#if NLS_DEBUG
    fprintf(stderr, "get_nls_derivs: T = %d, offset = %d\n", T, offset);
#endif

    for (j=0; j<spec->nparam; j++) {

	if (nls_calculate_deriv(spec, j)) {
	    return 1;
	}

#if NLS_DEBUG
	fprintf(stderr, "param[%d]: done nls_calculate_deriv\n", j);
#endif

	if (matrix_deriv(spec, j)) {
	    gretl_matrix *m = spec->params[j].dmat;
	    int i;

#if NLS_DEBUG
	    fprintf(stderr, "param[%d]: matrix at %p (%d x %d)\n", j, 
		    (void *) m, m->rows, m->cols);
#endif

	    for (i=0; i<m->cols; i++) {
		x = gretl_matrix_get(m, 0, i);
		for (t=0; t<T; t++) {
		    if (t > 0 && m->rows > 0) {
			x = gretl_matrix_get(m, t, i);
		    }
		    gi[t] = (spec->ci == MLE)? x : -x;
#if NLS_DEBUG > 1
		    fprintf(stderr, " set g[%d] = M(%d,%d) = %.14g\n", 
			    t, t, i, gi[t]);
#endif
		}
		if (g != NULL) {
		    gi += T;
		} else {
		    gi = *(++G) + offset;
		}
	    }
	} else if (scalar_deriv(spec, j)) {
	    /* transcribe from scalar var to array g */
	    x = gretl_scalar_get_value(spec->params[j].dname);
	    for (t=0; t<T; t++) {
		gi[t] = (spec->ci == MLE)? x : -x;
	    }
	    if (j < spec->nparam - 1) {
		/* advance the writing position */
		if (g != NULL) {
		    gi += T;
		} else {
		    gi = *(++G) + offset;
		}
	    }	    
	} else {
	    /* derivative is series */
	    int s, v = spec->params[j].dnum;

	    if (v == 0) {
		/* FIXME? */
		v = spec->params[j].dnum = spec->dinfo->v - 1;
	    }

#if NLS_DEBUG
	    fprintf(stderr, "param[%d]: dnum = %d\n", j, v);
#endif

	    /* transcribe from dataset to array g */
	    for (t=0; t<T; t++) {
		s = t + spec->t1;
		x = (*spec->Z)[v][s];
		gi[t] = (spec->ci == MLE)? x : -x;
#if NLS_DEBUG > 1
		fprintf(stderr, " set g[%d] = s->Z[%d][%d] = %.14g\n", 
			t, v, s, gi[t]);
#endif
	    }
	    if (j < spec->nparam - 1) {
		/* advance the writing position */
		if (g != NULL) {
		    gi += T;
		} else {
		    gi = *(++G) + offset;
		}
	    }
	}
    }

    return err;
}

/* analytical derivatives, used in the context of BFGS */

static int get_mle_gradient (double *b, double *g, int n, 
			     BFGS_CRIT_FUNC llfunc,
			     void *p)
{
    nlspec *spec = (nlspec *) p;
    gretl_matrix *m;
    double x;
    int i, j, k, t;
    int err = 0;

    update_coeff_values(b, spec);

#if ML_DEBUG
    fprintf(stderr, "get_mle_gradient\n");
#endif

    i = 0;

    for (j=0; j<spec->nparam; j++) {

	if (nls_calculate_deriv(spec, j)) {
	    return 1;
	}

#if ML_DEBUG
	fprintf(stderr, "mle: param %d (%s): done nls_calculate_deriv\n", 
		j, spec->params[j].name);
#endif

	if (matrix_deriv(spec, j)) {
	    m = spec->params[j].dmat;
#if ML_DEBUG > 1
	    gretl_matrix_print(m, "deriv matrix");
#endif
	    for (k=0; k<m->cols; k++) {
		g[i] = 0.0;
		for (t=0; t<m->rows; t++) {
		    x = gretl_matrix_get(m, t, k);
		    if (na(x)) {
			fprintf(stderr, "NA in gradient calculation\n");
			err = 1;
		    } else {
			g[i] += x;
		    }
		}
#if ML_DEBUG > 1
		fprintf(stderr, "set gradient g[%d] = %g\n", i, g[i]);
#endif
		i++;
	    }
	} else if (scalar_deriv(spec, j)) {
	    x = gretl_scalar_get_value(spec->params[j].dname);
	    if (na(x)) {
		fprintf(stderr, "NA in gradient calculation\n");
		err = 1;
	    } else {
		g[i] = x;
	    }
	} else {
	    /* the derivative must be a series */
	    int v = spec->params[j].dnum;

#if ML_DEBUG > 1
	    fprintf(stderr, "param[%d], dnum = %d\n", j, v);
#endif
	    g[i] = 0.0;

	    for (t=spec->t1; t<=spec->t2; t++) {
		x = (*spec->Z)[v][t];
#if ML_DEBUG > 1
		fprintf(stderr, "s->Z[%d][%d] = %g\n", v, t, x);
#endif
		if (na(x)) {
		    fprintf(stderr, "NA in gradient calculation\n");
		    err = 1;
		} else {
		    g[i] += x;
		}
	    }
#if ML_DEBUG > 1
	    fprintf(stderr, "set gradient g[%d] = %g\n", i, g[i]);
#endif
	    i++;
	} 
    }

    return err;
}

/* Compute auxiliary statistics and add them to the NLS 
   model struct */

static void add_stats_to_model (MODEL *pmod, nlspec *spec,
				const double **Z)
{
    int dv = spec->dv;
    double d, tss;
    int s, t;

    pmod->ess = spec->crit;
    pmod->sigma = sqrt(pmod->ess / (pmod->nobs - spec->ncoeff));
    
    pmod->ybar = gretl_mean(pmod->t1, pmod->t2, Z[dv]);
    pmod->sdy = gretl_stddev(pmod->t1, pmod->t2, Z[dv]);

    s = (spec->missmask != NULL)? 0 : pmod->t1;

    tss = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	d = Z[dv][s++] - pmod->ybar;
	tss += d * d;
    } 

    if (tss == 0.0) {
	pmod->rsq = pmod->adjrsq = NADBL;
    } else {
	pmod->rsq = 1.0 - pmod->ess / tss;
	pmod->adjrsq = 1.0 - (1.0 - pmod->rsq) * 
	    ((double) (pmod->nobs - 1) / (pmod->nobs - pmod->ncoeff));
    }
}

static int QML_vcv (nlspec *spec, gretl_matrix *V)
{
    gretl_matrix *Hinv = NULL;
    gretl_matrix *tmp = NULL;
    int k = V->rows;
    double x;
    int i, j;
    
    Hinv = gretl_matrix_alloc(k, k);
    tmp = gretl_matrix_alloc(k, k);

    if (Hinv == NULL || tmp == NULL) {
	free(Hinv);
	free(tmp);
	return E_ALLOC;
    }

    /* expand Hinv */
    for (i=0; i<k; i++) {
	for (j=0; j<=i; j++) {
	    x = spec->hessvec[ijton(i, j, k)];
	    gretl_matrix_set(Hinv, i, j, x);
	    gretl_matrix_set(Hinv, j, i, x);
	}
    }
	    
    /* form sandwich: V <- H^{-1} V H^{-1} */
    gretl_matrix_copy_values(tmp, V);
    gretl_matrix_qform(Hinv, GRETL_MOD_NONE, tmp,
		       V, GRETL_MOD_NONE);

    gretl_matrix_free(Hinv);
    gretl_matrix_free(tmp);

    return 0;
}

static double *mle_score_callback (const double *b, int i, void *p)
{
    nlspec *s = (nlspec *) p;

    update_coeff_values(b, s);
    nl_calculate_fvec(s);

    if (s->lvec != NULL) {
	return s->lvec->val;
    } else {
	return (*s->Z)[s->lhv] + s->t1;
    }
}

/* build the G matrix, given a final set of coefficient
   estimates, b, and a function for calculating the score
   vector, scorefun
*/

gretl_matrix *build_OPG_matrix (double *b, int k, int T,
				BFGS_SCORE_FUNC scorefun,
				void *data, int *err)
{
    double h = 1e-8;
#if ALT_OPG
    double d = 1.0e-4;
#endif
    gretl_matrix *G;
    const double *x;
    double bi0, x0;
    int i, t;

    G = gretl_matrix_alloc(k, T);
    if (G == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    gretl_matrix_zero(G);

    for (i=0; i<k; i++) {
	bi0 = b[i];
#if ALT_OPG
	h = d * bi0 + d * (b[i] == 0.0);
#endif
	b[i] = bi0 - h;
	x = scorefun(b, i, data);
	if (x == NULL) {
	    *err = E_NAN;
	    goto bailout;
	}
	for (t=0; t<T; t++) {
	    gretl_matrix_set(G, i, t, x[t]);
	}
	b[i] = bi0 + h;
	x = scorefun(b, i, data);
	if (x == NULL) {
	    *err = E_NAN;
	    goto bailout;
	}
	for (t=0; t<T; t++) {
	    x0 = gretl_matrix_get(G, i, t);
	    gretl_matrix_set(G, i, t, (x[t] - x0) / (2.0 * h));
	}
	b[i] = bi0;
#if NLS_DEBUG
	fprintf(stderr, "b[%d]: using %#.12g and %#.12g\n", i, bi0 - h, bi0 + h);
#endif
    }

#if NLS_DEBUG
    gretl_matrix_print(G, "Numerically estimated score");
#endif

 bailout:

    if (*err) {
	gretl_matrix_free(G);
	G = NULL;
    }

    return G;
}

/* add variance matrix based on OPG, (GG')^{-1}, or on QML sandwich,
   H^{-1} GG' H^{-1}
*/

static int mle_build_vcv (MODEL *pmod, nlspec *spec, int *vcvopt)
{
    gretl_matrix *G = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *m;
    double x = 0.0;
    int k = spec->ncoeff;
    int T = spec->nobs;
    int i, j, v, s, t;
    int err = 0;

    V = gretl_matrix_alloc(k, k);
    if (V == NULL) {
	return E_ALLOC;
    }

    if (numeric_mode(spec)) {
	G = build_OPG_matrix(spec->coeff, k, T, mle_score_callback, 
			     (void *) spec, &err);
	if (err) {
	    gretl_matrix_free(V);
	    return err;
	}
    } else {
	G = gretl_matrix_alloc(k, T);
	if (G == NULL) {
	    gretl_matrix_free(V);
	    return E_ALLOC;
	}
	j = 0;
	for (i=0; i<spec->nparam; i++) {
	    if (matrix_deriv(spec, i)) {
		m = spec->params[i].dmat;
		for (s=0; s<m->cols; s++) {
		    x = gretl_matrix_get(m, 0, s);
		    for (t=0; t<T; t++) {
			if (t > 0 && m->rows > 0) {
			    x = gretl_matrix_get(m, t, s);
			}
			gretl_matrix_set(G, j, t, x);
		    }
		    j++;
		}
	    } else if (scalar_deriv(spec, i)) {
		x = gretl_scalar_get_value(spec->params[i].dname);
		for (t=0; t<T; t++) {
		    gretl_matrix_set(G, j, t, x);
		}
		j++;
	    } else {
		v = spec->params[i].dnum;
		for (t=0; t<T; t++) {
		    x = (*spec->Z)[v][t + spec->t1];
		    gretl_matrix_set(G, j, t, x);
		}
		j++;
	    } 
	}		
    }

    gretl_matrix_multiply_mod(G, GRETL_MOD_NONE,
			      G, GRETL_MOD_TRANSPOSE,
			      V, GRETL_MOD_NONE);

    if ((spec->opt & OPT_R) && spec->hessvec != NULL) {
	/* robust option -> QML */
	err = QML_vcv(spec, V);
	*vcvopt = VCV_QML;
    } else {
	/* plain OPG */
	err = gretl_invert_symmetric_matrix(V);
	*vcvopt = VCV_OP;
    }

    if (!err) {
	err = gretl_model_write_vcv(pmod, V);
    }

    gretl_matrix_free(G);
    gretl_matrix_free(V);

    return err;
}

static int mle_add_vcv (MODEL *pmod, nlspec *spec)
{
    int i, k = spec->ncoeff;
    int vcvopt = VCV_OP;
    double x;
    int err = 0;

    if ((spec->opt & OPT_H) && spec->hessvec != NULL) {
	/* vcv based on Hessian */
	int n = (k * (k + 1)) / 2;

	for (i=0; i<n; i++) {
	    pmod->vcv[i] = spec->hessvec[i];
	}
	for (i=0; i<k; i++) {
	    x = pmod->vcv[ijton(i, i, k)];
	    pmod->sderr[i] = sqrt(x);
	}
	vcvopt = VCV_HESSIAN;
    } else {
	/* either OPG or QML */
	err = mle_build_vcv(pmod, spec, &vcvopt);
    }

    if (!err) {
	gretl_model_set_vcv_info(pmod, VCV_ML, vcvopt);
    }    

    return err;
}

/* NLS: add coefficient covariance matrix and standard errors 
   based on GNR */

static int add_nls_std_errs_to_model (MODEL *pmod)
{
    int i, k;

    if (pmod->vcv == NULL && makevcv(pmod, pmod->sigma)) {
	return E_ALLOC;
    }

    for (i=0; i<pmod->ncoeff; i++) {
	k = ijton(i, i, pmod->ncoeff);
	if (pmod->vcv[k] == 0.0) {
	    pmod->sderr[i] = 0.0;
	} else if (pmod->vcv[k] > 0.0) {
	    pmod->sderr[i] = sqrt(pmod->vcv[k]);
	} else {
	    pmod->sderr[i] = NADBL;
	}
    }

    return 0;
}

/* transcribe coefficient estimates into model struct */

static void add_coeffs_to_model (MODEL *pmod, double *coeff)
{
    int i;

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = coeff[i];
    }
}

/* convert (back) from '%s - (%s)' to '%s = %s' */

static char *adjust_saved_nlfunc (char *s)
{
    char *p = strchr(s, '-');

    if (p != NULL) {
	*p = '=';
    }

    p = strchr(s, '(');
    if (p != NULL) {
	shift_string_left(p, 1);
    }

    p = strrchr(s, ')');
    if (p != NULL) {
	*p = '\0';
    }  

    return s;
}

 /* Attach additional spcification info to make it possible to
    reconstruct the model from the saved state.  We need this
    for the "Modify model" option in the gretl GUI, and may
    also make use of it for bootstrapping NLS models.
 */

static int nl_model_add_spec_info (MODEL *pmod, nlspec *spec)
{
    const char *cmd = gretl_command_word(spec->ci);
    PRN *prn;
    char *buf;
    int i, err = 0;

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
    if (err) {
	return err;
    }

    pputs(prn, cmd);

    if (pmod->depvar != NULL) {
	pprintf(prn, " %s\n", pmod->depvar);
    } else {
	pputc(prn, '\n');
    }

    for (i=0; i<spec->naux; i++) {
	pprintf(prn, "%s\n", spec->aux[i]);
    }

    if (spec->ci == GMM) {
	/* orthog, weights */
	nlspec_print_gmm_info(spec, prn);
    }

    if (numeric_mode(spec)) {
	pprintf(prn, "params");
	for (i=0; i<spec->nparam; i++) {
	    pprintf(prn, " %s", spec->params[i].name);
	}
	pputc(prn, '\n');
    } else {
	for (i=0; i<spec->nparam; i++) {
	    pprintf(prn, "deriv %s = %s\n", spec->params[i].name,
		    spec->params[i].deriv);
	}
    }

    pprintf(prn, "end %s\n", cmd);

    buf = gretl_print_steal_buffer(prn);
    gretl_model_set_string_as_data(pmod, "nlinfo", buf);
    gretl_print_destroy(prn);

    return 0;
}

/* Called for all of NLS, MLE, GMM: attach to the model struct
   the names of the parameters and some other string info
*/

static int 
add_param_names_to_model (MODEL *pmod, nlspec *spec, const DATAINFO *pdinfo)
{
    char pname[VNAMELEN];
    int nc = pmod->ncoeff;
    int i, j, k, m, n;
    int err = 0;

    pmod->params = strings_array_new(nc);
    if (pmod->params == NULL) {
	return 1;
    }

    pmod->nparams = nc;

    if (spec->ci == NLS) {
	pmod->depvar = gretl_strdup(spec->nlfunc);
	if (pmod->depvar != NULL) {
	    adjust_saved_nlfunc(pmod->depvar);
	} else {
	    err = E_ALLOC;
	}
    } else if (spec->nlfunc != NULL) {
	n = strlen(spec->nlfunc) + strlen(spec->lhname);
	pmod->depvar = malloc(n + 4);
	if (pmod->depvar != NULL) {
	    sprintf(pmod->depvar, "%s = %s", spec->lhname, spec->nlfunc);
	} else {
	    err = E_ALLOC;
	}
    } 

    if (err) {
	free(pmod->params);
	return err;
    }    

    i = 0;
    for (j=0; j<spec->nparam; j++) {
	if (scalar_param(spec, j)) {
	    pmod->params[i++] = gretl_strdup(spec->params[j].name);
	} else {
	    m = spec->params[j].nc;
	    sprintf(pname, "%d", m + 1);
	    n = VNAMELEN - strlen(pname) - 3;
	    for (k=0; k<m; k++) {
		sprintf(pname, "%.*s[%d]", n, spec->params[j].name, k + 1);
		pmod->params[i++] = gretl_strdup(pname);
	    }
	}
    }

    if (spec->ci == NLS || gretl_in_gui_mode()) {
	err = nl_model_add_spec_info(pmod, spec);
    }

    return err;
}

static void 
add_fit_resid_to_model (MODEL *pmod, nlspec *spec, double *uhat, 
			const double **Z, const DATAINFO *pdinfo,
			int perfect)
{
    int t, j = 0;

    if (perfect) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    pmod->uhat[t] = 0.0;
	    pmod->yhat[t] = Z[spec->dv][t];
	}
    } else {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    pmod->uhat[t] = uhat[j];
	    pmod->yhat[t] = Z[spec->dv][t] - uhat[j];
	    j++;
	}
    }

    if (perfect || (spec->flags & NLS_AUTOREG)) {
	pmod->rho = pmod->dw = NADBL;
    }
}

/* this may be used later for generating out-of-sample forecasts --
   see nls_forecast() in forecast.c
*/

static int transcribe_nls_function (MODEL *pmod, const char *s)
{
    char *formula;
    int err = 0;

    /* skip "depvar - " */
    s += strcspn(s, " ") + 3;

    formula = gretl_strdup(s);
    if (s != NULL) {
	gretl_model_set_string_as_data(pmod, "nl_regfunc", formula); 
    } else {
	err = E_ALLOC;
    }

    return err;
}

#if NLS_DEBUG > 1
static void 
print_GNR_dataset (const int *list, double **gZ, DATAINFO *gdinfo)
{
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
    int t1 = gdinfo->t1;

    fprintf(stderr, "gdinfo->t1 = %d, gdinfo->t2 = %d\n",
	    gdinfo->t1, gdinfo->t2);
    gdinfo->t1 = 0;
    printdata(list, NULL, (const double **) gZ, gdinfo, OPT_O, prn);
    gdinfo->t1 = t1;
    gretl_print_destroy(prn);
}
#endif

/* Gauss-Newton regression to calculate standard errors for the NLS
   parameters (see Davidson and MacKinnon).  This model is taken
   as the basis for the model struct returned by the nls function,
   which is why we make the artificial dataset, gZ, full length.
*/

static MODEL GNR (nlspec *spec, const double **Z, DATAINFO *pdinfo, 
		  PRN *prn)
{
    double *uhat = spec->fvec;
    double **gZ = NULL;
    DATAINFO *gdinfo;
    int *glist;
    MODEL gnr;
    gretlopt lsqopt;
    int i, j, t, v;
    int T = spec->nobs;
    int iters = spec->iters;
    int perfect = 0;
    int err = 0;

#if NLS_DEBUG
    fprintf(stderr, "GNR: T = %d\n", T);
    v = pdinfo->v;
    fprintf(stderr, "dinfo->v = %d\n", v);
    fprintf(stderr, "Z = %p\n", (void *) Z);
    fprintf(stderr, "Z[%d] = %p\n", v-1, (void *) Z[v-1]);
#endif

    if (gretl_iszero(0, spec->nobs - 1, uhat)) {
	pputs(prn, _("Perfect fit achieved\n"));
	perfect = 1;
	for (t=0; t<spec->nobs; t++) {
	    uhat[t] = 1.0; /* will be adjusted later */
	}
	spec->crit = 0.0;
    }

    /* number of variables = 1 (const) + 1 (depvar) + spec->ncoeff
       (derivatives) */
    gdinfo = create_auxiliary_dataset(&gZ, spec->ncoeff + 2, pdinfo->n);
    if (gdinfo == NULL) {
	gretl_model_init(&gnr);
	gnr.errcode = E_ALLOC;
	return gnr;
    }

    /* transcribe sample info */
    gdinfo->t1 = spec->t1;
    gdinfo->t2 = spec->t2;
    gdinfo->pd = pdinfo->pd;
    gdinfo->structure = pdinfo->structure;

    glist = gretl_list_new(spec->ncoeff + 1);

    if (glist == NULL) {
	destroy_dataset(gZ, gdinfo);
	gretl_model_init(&gnr);
	gnr.errcode = E_ALLOC;
	return gnr;
    }

    j = 0;

    /* dependent variable (NLS residual) */
    glist[1] = 1;
    strcpy(gdinfo->varname[1], "gnr_y");
    for (t=0; t<gdinfo->n; t++) {
	if (t < gdinfo->t1 || t > gdinfo->t2) {
	    gZ[1][t] = NADBL;
	} else {
	    gZ[1][t] = uhat[j++];
	}
    }

    for (i=0; i<spec->ncoeff; i++) {
	/* independent vars: derivatives wrt NLS params */
	v = i + 2;
	glist[v] = v;
	sprintf(gdinfo->varname[v], "gnr_x%d", i + 1);
    }

    if (analytic_mode(spec)) {
	for (i=0; i<spec->ncoeff; i++) {
	    v = i + 2;
	    for (t=0; t<gdinfo->t1; t++) {
		gZ[v][t] = NADBL;
	    }
	    for (t=gdinfo->t2; t<gdinfo->n; t++) {
		gZ[v][t] = NADBL;
	    }
	}
#if NLS_DEBUG
	fprintf(stderr, "GNR: calling get_nls_derivs\n");
	fprintf(stderr, "Z[%d] = %p\n", pdinfo->v-1, (void *) Z[pdinfo->v-1]);
#endif
	get_nls_derivs(spec->nparam, T, spec->t1, NULL, gZ + 2, spec);
    } else {
	for (i=0; i<spec->ncoeff; i++) {
	    v = i + 2;
	    j = T * i; /* calculate offset into jac */
	    for (t=0; t<gdinfo->n; t++) {
		if (t < gdinfo->t1 || t > gdinfo->t2) {
		    gZ[v][t] = NADBL;
		} else {
		    gZ[v][t] = spec->jac[j++];
		}
	    }
	}
    }

#if NLS_DEBUG > 1
    print_GNR_dataset(glist, gZ, gdinfo);
#endif

    lsqopt = OPT_A | OPT_B; /* OPT_B: don't calculate R^2 */
    if (spec->opt & OPT_R) {
	/* robust variance matrix, if wanted */
	lsqopt |= OPT_R;
    }

    gnr = lsq(glist, &gZ, gdinfo, OLS, lsqopt);

#if NLS_DEBUG
    gnr.name = gretl_strdup("GNR for NLS");
    printmodel(&gnr, gdinfo, OPT_NONE, prn);
    free(gnr.name);
    gnr.name = NULL;
#endif

    if (gnr.errcode) {
	pputs(prn, _("In Gauss-Newton Regression:\n"));
	errmsg(gnr.errcode, prn);
	err = 1;
    } 

    if (gnr.list[0] != glist[0]) {
	gnr.errcode = E_JACOBIAN;
    }

    if (gnr.errcode == 0) {
	gnr.ci = spec->ci;
	add_stats_to_model(&gnr, spec, Z);
	if (add_nls_std_errs_to_model(&gnr)) {
	    gnr.errcode = E_ALLOC;
	}
    }

    if (gnr.errcode == 0) {
	ls_criteria(&gnr);
	gnr.fstt = gnr.chisq = NADBL;
	add_coeffs_to_model(&gnr, spec->coeff);
	add_param_names_to_model(&gnr, spec, pdinfo);
	add_fit_resid_to_model(&gnr, spec, uhat, Z, pdinfo, perfect);
	gnr.list[1] = spec->dv;

	/* set relevant data on model to be shipped out */
	gretl_model_set_int(&gnr, "iters", iters);
	gretl_model_set_double(&gnr, "tol", spec->tol);
	transcribe_nls_function(&gnr, spec->nlfunc);
	if (spec->flags & NLS_AUTOREG) {
	    gretl_model_set_int(&gnr, "dynamic", 1);
	}
    }

    destroy_dataset(gZ, gdinfo);
    free(glist);

    return gnr;
}

/* allocate space to copy info into model struct */

static int nl_model_allocate (MODEL *pmod, nlspec *spec)
{
    int k = spec->ncoeff;
    int nvc = (k * k + k) / 2;

    pmod->coeff = malloc(k * sizeof *pmod->coeff);
    pmod->sderr = malloc(k * sizeof *pmod->sderr);
    pmod->vcv = malloc(nvc * sizeof *pmod->vcv);

    if (pmod->coeff == NULL || pmod->sderr == NULL || pmod->vcv == NULL) {
	pmod->errcode = E_ALLOC;
    } else {
	pmod->ncoeff = k;
    }

    return pmod->errcode;
}

/* Work up the results of estimation into the form of a gretl
   MODEL struct.  This is used for MLE and GMM; in the case of
   NLS the Gauss-Newton regression provides the basis for the
   final model structure.
*/

static int make_nl_model (MODEL *pmod, nlspec *spec, 
			  const DATAINFO *pdinfo,
			  gretlopt opt)
{
    nl_model_allocate(pmod, spec);
    if (pmod->errcode) {
	return pmod->errcode;
    }

    pmod->t1 = spec->t1;
    pmod->t2 = spec->t2;
    pmod->nobs = spec->nobs;
    
    /* hmm */
    pmod->dfn = pmod->ncoeff;
    pmod->dfd = pmod->nobs - pmod->ncoeff;

    pmod->ci = spec->ci;

    if (spec->ci == MLE) {
	pmod->lnL = spec->crit;
    } else {
	/* GMM */
	pmod->lnL = NADBL;
    }

    mle_criteria(pmod, 0);

    add_coeffs_to_model(pmod, spec->coeff);

    if (!(opt & OPT_G)) {
	/* not IVREG via GMM */
	pmod->errcode = add_param_names_to_model(pmod, spec, pdinfo);
    }

    if (!pmod->errcode) {
	if (pmod->ci == MLE) {
	    pmod->errcode = mle_add_vcv(pmod, spec);
	} else {
	    pmod->errcode = gmm_add_vcv(pmod, spec);
	}
    }

    if (!pmod->errcode && pmod->ci == GMM) {
	maybe_add_gmm_residual(pmod, spec, pdinfo);
    }

    if (!pmod->errcode) {
	gretl_model_set_int(pmod, "fncount", spec->fncount);
	gretl_model_set_int(pmod, "grcount", spec->grcount);
	gretl_model_set_double(pmod, "tol", spec->tol);
    }

    /* mask invalid stats */
    pmod->sigma = NADBL;
    pmod->rsq = pmod->adjrsq = NADBL;
    pmod->fstt = pmod->chisq = NADBL;
    pmod->rho = pmod->dw = NADBL;

    return pmod->errcode;
}

static int add_nls_coeffs (MODEL *pmod, nlspec *spec)
{
    pmod->ncoeff = spec->ncoeff;
    pmod->full_n = 0;

    pmod->errcode = gretl_model_allocate_storage(pmod);

    if (!pmod->errcode) {
	add_coeffs_to_model(pmod, spec->coeff);
    }

    return pmod->errcode;
}

/* free up resources associated with the nlspec struct */

static void clear_nlspec (nlspec *spec)
{
    int i;

    if (spec == NULL) {
	return;
    }

    if (spec->params != NULL) {
	for (i=0; i<spec->nparam; i++) {
	    free(spec->params[i].deriv);
	    spec->params[i].vec = NULL;
	    spec->params[i].dmat = NULL;
	}
	free(spec->params);
	spec->params = NULL;
    }
    spec->nparam = 0;

    free(spec->fvec);
    spec->fvec = NULL;

    free(spec->jac);
    spec->jac = NULL;

    free(spec->coeff);
    spec->coeff = NULL;
    spec->ncoeff = 0;
    spec->nvec = 0;

    if (spec->aux != NULL) {
	for (i=0; i<spec->naux; i++) {
	    free(spec->aux[i]);
	}
	free(spec->aux);
	spec->aux = NULL;
    }
    spec->naux = 0;

    if (spec->genrs != NULL) {
	for (i=0; i<spec->ngenrs; i++) {
	    destroy_genr(spec->genrs[i]);
	}
	free(spec->genrs);
	spec->genrs = NULL;
    }
    spec->ngenrs = 0;
    spec->generr = 0;

    free(spec->nlfunc);
    spec->nlfunc = NULL;

    free(spec->hessvec);
    spec->hessvec = NULL;

    spec->flags = 0;
    spec->opt = OPT_NONE;

    spec->dv = 0;
    spec->lhv = 0;
    spec->lvec = NULL;

    spec->iters = 0;
    spec->fncount = 0;
    spec->grcount = 0;

    spec->t1 = spec->t2 = 0;
    spec->real_t1 = spec->real_t2 = 0;
    spec->nobs = 0;

    spec->Z = NULL;
    spec->dinfo = NULL;
    spec->prn = NULL;

    if (spec->oc != NULL) {
	oc_set_destroy(spec->oc);
	spec->oc = NULL;
    }

    free(spec->missmask);
    spec->missmask = NULL;
}

/* 
   Next block: functions that interface with minpack.

   The details below may be obscure, but here's the basic idea: The
   minpack functions are passed an array ("fvec") that holds the
   calculated values of the function to be minimized, at a given value
   of the parameters, and also (in the case of analytic derivatives)
   an array holding the Jacobian ("jac").  Minpack is also passed a
   callback function that will recompute the values in these arrays,
   given a revised vector of parameter estimates.

   As minpack does its iterative thing, at each step it invokes the
   callback function, supplying its updated parameter estimates and
   saying via a flag variable ("iflag") whether it wants the function
   itself or the Jacobian re-evaluated.

   The libgretl strategy involves holding all the relevant values (nls
   residual, nls derivatives, and nls parameters) as variables in a
   gretl dataset (Z-array and datainfo-struct pair).  The callback
   function that we supply to minpack first transcribes the revised
   parameter estimates into the dataset, then invokes genr() to
   recalculate the residual and derivatives, then transcribes the
   results back into the fvec and jac arrays.
*/

/* callback for lm_calculate (below) to be used by minpack */

static int nls_calc (integer *m, integer *n, double *x, double *fvec, 
		     double *jac, integer *ldjac, integer *iflag,
		     void *p)
{
    nlspec *s = (nlspec *) p;

#if NLS_DEBUG
    fprintf(stderr, "nls_calc called by minpack with iflag = %d\n", 
	    (int) *iflag);
#endif

    /* write current coefficient values into dataset */
    update_coeff_values(x, s);

    if (*iflag == 1) {
	/* calculate function at x, results into fvec */
	if (nl_function_calc(fvec, p)) {
	    *iflag = -1;
	} 
    } else if (*iflag == 2) {
	/* calculate jacobian at x, results into jac */
	if (get_nls_derivs(*n, *m, 0, jac, NULL, p)) {
	    *iflag = -1; 
	}
    }

    return 0;
}

/* in case the user supplied analytical derivatives for the
   parameters, check them for sanity */

static int check_derivatives (integer m, integer n, double *x,
			      double *fvec, double *jac,
			      integer ldjac, PRN *prn,
			      nlspec *spec)
{
#if NLS_DEBUG > 1
    int T = spec->nobs * spec->ncoeff;
#endif
    integer mode, iflag;
    doublereal *xp = NULL;
    doublereal *err = NULL;
    doublereal *fvecp = NULL;
    int i, badcount = 0, zerocount = 0;

    xp = malloc(n * sizeof *xp);
    err = malloc(m * sizeof *err);
    fvecp = malloc(m * sizeof *fvecp);

    if (xp == NULL || err == NULL || fvecp == NULL) {
	free(err);
	free(xp);
	free(fvecp);
	return 1;
    }

#if NLS_DEBUG > 1
    fprintf(stderr, "\nchkder, starting: m=%d, n=%d, ldjac=%d\n",
	    (int) m, (int) n, (int) ldjac);
    for (i=0; i<spec->ncoeff; i++) {
	fprintf(stderr, "x[%d] = %.14g\n", i, x[i]);
    }    
    for (i=0; i<spec->nobs; i++) {
	fprintf(stderr, "fvec[%d] = %.14g\n", i, fvec[i]);
    }
#endif

    /* mode 1: x contains the point of evaluation of the function; on
       output xp is set to a neighboring point. */
    mode = 1;
    chkder_(&m, &n, x, fvec, jac, &ldjac, xp, fvecp, &mode, err);

    /* calculate gradient */
    iflag = 2;
    nls_calc(&m, &n, x, fvec, jac, &ldjac, &iflag, spec);
    if (iflag == -1) goto chkderiv_abort;

#if NLS_DEBUG > 1
    fprintf(stderr, "\nchkder, calculated gradient\n");
    for (i=0; i<T; i++) {
	fprintf(stderr, "jac[%d] = %.14g\n", i, jac[i]);
    }
#endif

    /* calculate function, at neighboring point xp */
    iflag = 1;
    nls_calc(&m, &n, xp, fvecp, jac, &ldjac, &iflag, spec);
    if (iflag == -1) goto chkderiv_abort; 

    /* mode 2: on input, fvec must contain the functions, the rows of
       fjac must contain the gradients evaluated at x, and fvecp must
       contain the functions evaluated at xp.  On output, err contains
       measures of correctness of the respective gradients.
    */
    mode = 2;
    chkder_(&m, &n, x, fvec, jac, &ldjac, xp, fvecp, &mode, err);

#if NLS_DEBUG > 1
    fprintf(stderr, "\nchkder, done mode 2:\n");
    for (i=0; i<m; i++) {
	fprintf(stderr, "%d: fvec = %.14g, fvecp = %.14g, err = %g\n", i, 
		fvec[i], fvecp[i], err[i]);
    }
#endif

    /* examine "err" vector */
    for (i=0; i<m; i++) {
	if (err[i] == 0.0) {
	    zerocount++;
	} else if (err[i] < 0.35) {
	    badcount++;
	}
    }

    if (zerocount > 0) {
	strcpy(gretl_errmsg, 
	       _("NLS: The supplied derivatives seem to be incorrect"));
	fprintf(stderr, _("%d out of %d tests gave zero\n"), zerocount, (int) m);
    } else if (badcount > 0) {
	pputs(prn, _("Warning: The supplied derivatives may be incorrect, or perhaps\n"
		     "the data are ill-conditioned for this function.\n"));
	pprintf(prn, _("%d out of %d gradients looked suspicious.\n\n"),
		badcount, (int) m);
    }

 chkderiv_abort:
    free(xp);
    free(err);
    free(fvecp);

    return (zerocount > m / 4);
}

/* drivers for BFGS code below, first for MLE */

static int mle_calculate (nlspec *s, PRN *prn)
{
    int err = 0;

    if (analytic_mode(s)) {
	integer m = s->nobs;
	integer n = s->ncoeff;
	integer ldjac = m; 

	if (!(s->opt & OPT_G)) {
	    err = check_derivatives(m, n, s->coeff, s->fvec, 
				    s->jac, ldjac, prn, s);
	}
    }

    if (!err) {
	BFGS_GRAD_FUNC gradfun = (analytic_mode(s))?
	    get_mle_gradient : NULL;
	int maxit = 500;

	err = BFGS_max(s->coeff, s->ncoeff, maxit, s->tol, 
		       &s->fncount, &s->grcount, 
		       get_mle_ll, C_LOGLIK, gradfun, s,
		       s->opt, s->prn);
	if (!err && (s->opt & (OPT_H | OPT_R))) {
	    /* doing Hessian or QML covariance matrix */
	    s->hessvec = numerical_hessian(s->coeff, s->ncoeff, 
					   get_mle_ll, s, &err);
	}
    }

    return err;    
}

/* Minpack driver for the case where analytical derivatives have been
   supplied */

static int lm_calculate (nlspec *spec, PRN *prn)
{
    integer info, lwa;
    integer m, n, ldjac;
    integer *ipvt;
    doublereal *wa;
    int err = 0;

    m = spec->nobs;    /* number of observations */
    n = spec->ncoeff;  /* number of coefficients */
    lwa = 5 * n + m;   /* work array size */
    ldjac = m;         /* leading dimension of jac array */

    wa = malloc(lwa * sizeof *wa);
    ipvt = malloc(n * sizeof *ipvt);

    if (wa == NULL || ipvt == NULL) {
	err = E_ALLOC;
	goto nls_cleanup;
    }

    err = check_derivatives(m, n, spec->coeff, spec->fvec, spec->jac, 
			    ldjac, prn, spec);
    if (err) {
	goto nls_cleanup; 
    }

    /* note: maxfev is automatically set to 100*(n + 1) */

    /* call minpack */
    lmder1_(nls_calc, &m, &n, spec->coeff, spec->fvec, spec->jac, &ldjac, 
	    &spec->tol, &info, ipvt, wa, &lwa, spec);

    switch ((int) info) {
    case -1: 
	err = 1;
	break;
    case 0:
	strcpy(gretl_errmsg, _("Invalid NLS specification"));
	err = 1;
	break;
    case 1:
    case 2:
    case 3:
    case 4:
	break;
    case 5:
    case 6:
    case 7:
	sprintf(gretl_errmsg, 
		_("NLS: failed to converge after %d iterations"),
		spec->iters);
	err = 1;
	break;
    default:
	break;
    }

 nls_cleanup:

    free(wa);
    free(ipvt);

    return err;    
}

/* callback for lm_approximate (below) to be used by minpack */

static int 
nls_calc_approx (integer *m, integer *n, double *x, double *fvec,
		 integer *iflag, void *p)
{
    int erru, errc, err;

    /* write current parameter values into dataset Z */
    err = erru = update_coeff_values(x, p);

    /* calculate function at x, results into fvec */  
    if (!err) {
	err = errc = nl_function_calc(fvec, p);
    }
    
    if (err) {
	/* flag error to minpack */
	fprintf(stderr, "nls_calc_approx: got error %d from %s\n", 
		err, (erru)? "update_coeff_values" : 
		"nl_function_calc");
	*iflag = -1;
    }

    return 0;
}

/* Minpack driver for the case where the Jacobian must be approximated
   numerically */

static int lm_approximate (nlspec *spec, PRN *prn)
{
    integer info, m, n, ldjac;
    integer maxfev, mode = 1, nprint = 0, nfev = 0;
    integer iflag = 0;
    integer *ipvt;
    doublereal gtol = 0.0;
    doublereal epsfcn = 0.0, factor = 100.;
    doublereal *dspace, *diag, *qtf;
    doublereal *wa1, *wa2, *wa3, *wa4;
    int err = 0;
    
    m = spec->nobs;              /* number of observations */
    n = spec->ncoeff;            /* number of parameters */
    ldjac = m;                   /* leading dimension of jac array */

    maxfev = 200 * (n + 1);      /* max iterations */

    dspace = malloc((5 * n + m) * sizeof *dspace);
    ipvt = malloc(n * sizeof *ipvt);

    if (dspace == NULL || ipvt == NULL) {
	err = E_ALLOC;
	goto nls_cleanup;
    }

    diag = dspace;
    qtf = diag + n;
    wa1 = qtf + n;
    wa2 = wa1 + n;
    wa3 = wa2 + n;
    wa4 = wa3 + n;

    /* call minpack */
    lmdif_(nls_calc_approx, &m, &n, spec->coeff, spec->fvec, 
	   &spec->tol, &spec->tol, &gtol, &maxfev, &epsfcn, diag, 
	   &mode, &factor, &nprint, &info, &nfev, spec->jac, &ldjac, 
	   ipvt, qtf, wa1, wa2, wa3, wa4, spec);

    spec->iters = nfev;

    switch ((int) info) {
    case -1: 
	err = 1;
	break;
    case 0:
	strcpy(gretl_errmsg, _("Invalid NLS specification"));
	err = 1;
	break;
    case 1:
    case 2:
    case 3:
    case 4:
	break;
    case 5:
    case 6:
    case 7:
    case 8:
	sprintf(gretl_errmsg, 
		_("NLS: failed to converge after %d iterations"),
		spec->iters);
	err = 1;
	break;
    default:
	break;
    }

    if (!err) {
	double ess = spec->crit;
	int iters = spec->iters;
	gretlopt opt = spec->opt;

	spec->opt = OPT_NONE;

	/* call minpack again */
	fdjac2_(nls_calc_approx, &m, &n, spec->coeff, spec->fvec, 
		spec->jac, &ldjac, &iflag, &epsfcn, wa4, spec);
	spec->crit = ess;
	spec->iters = iters;
	spec->opt = opt;
    }

 nls_cleanup:

    free(dspace);
    free(ipvt);

    return err;    
}

/* below: various public functions */

/**
 * nlspec_add_param_with_deriv:
 * @spec: pointer to nls specification.
 * @dstr: string specifying a derivative with respect to a
 *   parameter of the regression function.
 * @Z: data array.
 * @pdinfo: information on dataset.
 *
 * Adds an analytical derivative to @spec.  This pointer must
 * have previously been obtained by a call to nlspec_new().
 * The required format for @dstr is "%varname = %formula", where
 * %varname is the name of the (scalar) variable holding the parameter
 * in question, and %formula is an expression, of the sort that
 * is fed to gretl's %genr command, giving the derivative of the
 * regression function in @spec with respect to the parameter.
 * The variable holding the parameter must be already present in
 * the dataset.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int 
nlspec_add_param_with_deriv (nlspec *spec, const char *dstr,
			     const double **Z, const DATAINFO *pdinfo)
{
    const char *p = dstr;
    char *name = NULL;
    char *deriv = NULL;
    int err = 0;

    if (spec->ci == GMM) {
	strcpy(gretl_errmsg, _("Analytical derivatives cannot be used with GMM"));
	return E_DATA;
    }

    if (!strncmp(p, "deriv ", 6)) {
	/* make starting with "deriv" optional */
	p += 6;
    }

    err = equation_get_lhs_and_rhs(p, &name, &deriv);
    if (err) {
	fprintf(stderr, "parse error in deriv string: '%s'\n", dstr);
	return E_PARSE;
    }

    err = check_param_name(name);
    
    if (!err) {
	err = nlspec_push_param(spec, name, deriv);
	if (err) {
	    free(deriv);
	    deriv = NULL;
	}
    }

    free(name);

    if (!err) {
	set_analytic_mode(spec);
    }

    return err;
}

static int screen_bad_aux (const char *line, const DATAINFO *pdinfo)
{
    int n = gretl_namechar_spn(line);
    char word[FN_NAMELEN];
    int ci, err = E_DATA;

    *word = '\0';

    if (n > 0) {
	if (n > FN_NAMELEN - 1) {
	    n = FN_NAMELEN - 1;
	}
	strncat(word, line, n);
    } 

    ci = gretl_command_number(word);

    if (ci == GENR || ci == PRINT) {
	err = 0;
    } else if (plausible_genr_start(line, pdinfo)) {
	err = 0;
    } else if (get_user_function_by_name(word)) {
	err = 0;
    } else {
	sprintf(gretl_errmsg, _("command '%s' not valid in this context"), 
		word);
    }

    return err;
}

/**
 * nlspec_add_aux:
 * @spec: pointer to nls specification.
 * @s: string specifying an auxiliary command (primarily
 * for use in calculating function or derivatives).
 * @pdinfo: pointer to dataset information.
 *
 * Adds the specification of an auxiliary command to @spec, 
 * which pointer must have previously been obtained by a call 
 * to nlspec_new().
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int nlspec_add_aux (nlspec *spec, const char *s, const DATAINFO *pdinfo)
{
    int err;

#if NLS_DEBUG
    fprintf(stderr, "nlspec_add_aux: s = '%s'\n", s);
#endif

    err = screen_bad_aux(s, pdinfo);

    if (!err) {
	err = strings_array_add(&spec->aux, &spec->naux, s);
    }

    return err;
}

/**
 * nlspec_set_regression_function:
 * @spec: pointer to nls specification.
 * @fnstr: string specifying nonlinear regression function.
 * @pdinfo: information on dataset.
 *
 * Adds the regression function to @spec.  This pointer must
 * have previously been obtained by a call to nlspec_new().
 * The required format for @fnstr is "%varname = %formula", where
 * %varname is the name of the dependent variable and %formula
 * is an expression of the sort that is fed to gretl's %genr command.
 * The dependent variable must be already present in the
 * dataset.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int 
nlspec_set_regression_function (nlspec *spec, const char *fnstr, 
				const DATAINFO *pdinfo)
{    
    const char *p = fnstr;
    char *vname = NULL;
    char *rhs = NULL;
    int flen, err = 0;

    if (spec->nlfunc != NULL) {
	free(spec->nlfunc);
	spec->nlfunc = NULL;
    }

    spec->dv = 0;

    if (!strncmp(p, "nls ", 4) || 
	!strncmp(p, "mle ", 4) ||
	!strncmp(p, "gmm ", 4)) {
	p += 4;
    } else if (!strncmp(p, "gmm", 3)) {
	p += 3;
    }

    if (spec->ci == GMM && string_is_blank(p)) {
	/* GMM: we don't insist on a function on the first line */
	return 0;
    }

    if (equation_get_lhs_and_rhs(p, &vname, &rhs)) { 
	sprintf(gretl_errmsg, _("parse error in '%s'\n"), fnstr);
	err =  E_PARSE;
    } else if (spec->ci == NLS) {
	spec->dv = series_index(pdinfo, vname);
	if (spec->dv == pdinfo->v) {
	    sprintf(gretl_errmsg, _("Unknown variable '%s'"), vname);
	    err = E_UNKVAR;
	}
    } else {
	*spec->lhname = '\0';
	strncat(spec->lhname, vname, VNAMELEN - 1);
    } 

    if (!err) {
	if (spec->ci == MLE || spec->ci == GMM) {
	    spec->nlfunc = gretl_strdup(rhs);
	} else {
	    flen = strlen(vname) + strlen(rhs) + 6;
	    spec->nlfunc = malloc(flen);
	    if (spec->nlfunc != NULL) {
		/* the equation defining the NLS residual */
		sprintf(spec->nlfunc, "%s - (%s)", vname, rhs);
	    }
	}

	if (spec->nlfunc == NULL) {
	    err = E_ALLOC;
	} 
    }

    free(vname);
    free(rhs);

    return err;
}

/**
 * nlspec_set_t1_t2:
 * @spec: pointer to nls specification.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Sets the sample range for estimation of @spec.  This pointer must
 * have previously been obtained by a call to nlspec_new().
 */

void nlspec_set_t1_t2 (nlspec *spec, int t1, int t2)
{
    if (spec != NULL) {
	spec->t1 = t1;
	spec->t2 = t2;
	spec->nobs = t2 - t1 + 1;
    }
}

#define param_line(s) (!strncmp(s, "deriv ", 6) || \
                       !strncmp(s, "params ", 7) || \
                       !strncmp(s, "orthog ", 7) || \
                       !strncmp(s, "weights ", 8))

#define cmd_start(s) (!strncmp(s, "nls ", 4) || \
                      !strncmp(s, "mle ", 4) || \
                      !strncmp(s, "gmm ", 4) || \
                      !strcmp(s, "gmm")) 

/**
 * nl_parse_line:
 * @ci: %NLS, %MLE or %GMM.
 * @line: string containing information to be added to a 
 * nonlinear model specification.
 * @Z: data array.
 * @pdinfo: information on dataset.
 * @prn: gretl printing struct (for warning messages).
 *
 * This function is used to create the specification of a
 * nonlinear model, to be estimated via #nl_model (i.e.
 * via NLS, MLE or GMM).  It can be called several times to
 * build up the required information.  See the manual
 * entries for nls, mle and gmm for details.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int nl_parse_line (int ci, const char *line, const double **Z,
		   const DATAINFO *pdinfo, PRN *prn)
{
    nlspec *s = &private_spec;
    int err = 0;

    s->ci = ci;

#if NLS_DEBUG
    fprintf(stderr, "nls_parse_line: '%s'\n", line);
#endif

    if (cmd_start(line)) {
	if (s->nlfunc != NULL) {
	    clear_nlspec(s);
	}
	err = nlspec_set_regression_function(s, line, pdinfo);
	if (!err) {
	    nlspec_set_t1_t2(s, pdinfo->t1, pdinfo->t2);
	}	
    } else if (param_line(line)) {
	if (s->nlfunc == NULL && s->ci != GMM) {
	    strcpy(gretl_errmsg, _("No regression function has been specified"));
	    err = E_PARSE;
	} else {
	    if (!strncmp(line, "deriv", 5)) {
		if (numeric_mode(s) && s->params != NULL) {
		    strcpy(gretl_errmsg, _("You cannot supply both a \"params\" "
					   "line and analytical derivatives"));
		    err = E_PARSE;
		} else {
		    err = nlspec_add_param_with_deriv(s, line, Z, pdinfo);
		}
	    } else if (!strncmp(line, "params", 6)) {
		if (numeric_mode(s)) {
		    err = nlspec_add_params_from_line(s, line + 6, Z, pdinfo);
		} else {
		    pprintf(prn, _("Analytical derivatives supplied: "
				   "\"params\" line will be ignored"));
		    pputc(prn, '\n');
		}
	    } else if (!strncmp(line, "orthog", 6)) {
		err = nlspec_add_orthcond(s, line + 6, Z, pdinfo);
	    } else if (!strncmp(line, "weights", 7)) {
		err = nlspec_add_weights(s, line + 7);
	    }
	}
    } else {
	err = nlspec_add_aux(s, line, pdinfo);
    }

    if (err) {
	/* remember to clean up! */
	clear_nlspec(s);
    }

    return err;
}

static double default_nls_toler;

/**
 * get_default_nls_toler:
 *
 * Returns: the default value used in the convergence criterion
 * for estimation of models using nonlinear least squares.
 */

double get_default_nls_toler (void)
{
    if (default_nls_toler == 0.0) {
	default_nls_toler = pow(dpmpar_(&one), .75);
    }

    return default_nls_toler;
}

static int check_spec_requirements (nlspec *spec)
{
    if (spec->nparam < 1) {
	strcpy(gretl_errmsg, _("No parameters have been specified"));
	return 1;
    }

    if (spec->ci == GMM) {
	return check_gmm_requirements(spec);
    }

    return 0;
}

/* remedial treatment for an NLS model estimated using
   a trimmed dataset, to avoid missing observations
*/

static int nls_model_fix_sample (MODEL *pmod,
				 nlspec *spec, 
				 DATAINFO *pdinfo)
{
    double *uhat = NULL;
    double *yhat = NULL;

    if (pmod->uhat != NULL) {
	uhat = malloc(pdinfo->n * sizeof *uhat);
	if (uhat == NULL) {
	    return E_ALLOC;
	}
    }

    if (pmod->yhat != NULL) {
	yhat = malloc(pdinfo->n * sizeof *yhat);
	if (yhat == NULL) {
	    free(uhat);
	    return E_ALLOC;
	}
    }    

    if (uhat != NULL || yhat != NULL) {
	int t, s = 0;

	for (t=0; t<pdinfo->n; t++) {
	    if (t < spec->real_t1 || t > spec->real_t2 ||
		spec->missmask[t] == '1') {
		if (uhat != NULL) uhat[t] = NADBL;
		if (yhat != NULL) yhat[t] = NADBL;
	    } else {
		if (uhat != NULL) uhat[t] = pmod->uhat[s];
		if (yhat != NULL) yhat[t] = pmod->yhat[s];
		s++;
	    }
	}

	if (uhat != NULL) {
	    free(pmod->uhat);
	    pmod->uhat = uhat;
	}

	if (yhat != NULL) {
	    free(pmod->yhat);
	    pmod->yhat = yhat;
	}
    }

    pmod->full_n = pdinfo->n;
    pmod->t1 = spec->real_t1;
    pmod->t2 = spec->real_t2;

    pmod->missmask = spec->missmask;
    spec->missmask = NULL;

    return 0;
}

/* static function providing the real content for the two public
   wrapper functions below: does NLS, MLE or GMM */

static MODEL real_nl_model (nlspec *spec, double ***pZ, DATAINFO *pdinfo, 
			    gretlopt opt, PRN *prn)
{
    MODEL nlmod;
    int origv = pdinfo->v;
    int i, t, err = 0;

#if NLS_DEBUG
    fprintf(stderr, "real_nl_model: starting, pdinfo->v = %d\n", origv);
#endif

    gretl_model_init(&nlmod);
    gretl_model_smpl_init(&nlmod, pdinfo);

    if (spec == NULL) {
	/* we use the static spec composed via nl_parse_line() */
	spec = &private_spec;
    }

    if (spec->nlfunc == NULL && spec->ci != GMM) {
	strcpy(gretl_errmsg, _("No function has been specified"));
	nlmod.errcode = E_PARSE;
	goto bailout;
    } 

    spec->opt = opt;

    /* make pZ, pdinfo and prn available via spec */
    spec->Z = pZ;
    spec->dinfo = pdinfo;
    spec->prn = prn;

    if (opt & OPT_N) {
	/* ignore any analytical derivs */
	set_numeric_mode(spec);
    }

    if (numeric_mode(spec) && spec->nparam == 0 && spec->ci != GMM) {
	gretl_errmsg_sprintf(_("%s: no parameters were specified"),
			     gretl_command_word(spec->ci));
	nlmod.errcode = E_DATA;
	goto bailout;
    }

    if (check_spec_requirements(spec)) {
	nlmod.errcode = E_PARSE;
	goto bailout;
    } 

    if (spec->ci == GMM) {
	err = gmm_missval_check_etc(spec);
    } else {
	err = nl_missval_check(spec);
    }

    if (err) {
	nlmod.errcode = err;
	goto bailout;
    }

    /* allocate arrays to be passed to minpack */
    spec->fvec = malloc(spec->nobs * sizeof *spec->fvec);
    spec->jac = malloc(spec->nobs * spec->ncoeff * sizeof *spec->jac);

    if (spec->fvec == NULL || spec->jac == NULL) {
	nlmod.errcode = E_ALLOC;
	goto bailout;
    }

    if (spec->lvec != NULL) {
	for (t=0; t<spec->nobs; t++) {
	    spec->fvec[t] = spec->lvec->val[t];
	}
    } else {
	i = 0;
	for (t=spec->t1; t<=spec->t2; t++) {
	    spec->fvec[i++] = (*spec->Z)[spec->lhv][t];
	}
    }

    /* get tolerance from user setting or default */
    if (USES_BFGS(spec->ci)) {
	spec->tol = libset_get_double(BFGS_TOLER);
    } else {
	spec->tol = libset_get_double(NLS_TOLER);
    }

    pputs(prn, (numeric_mode(spec))?
	  _("Using numerical derivatives\n") :
	  _("Using analytical derivatives\n"));

    /* now start the actual calculations */

    if (spec->ci == MLE) {
	err = mle_calculate(spec, prn);
    } else if (spec->ci == GMM) {
	err = gmm_calculate(spec, prn);
    } else {
	/* NLS: invoke the appropriate minpack driver function */
	if (numeric_mode(spec)) {
	    err = lm_approximate(spec, prn);
	} else {
	    err = lm_calculate(spec, prn);
	    if (err) {
		fprintf(stderr, "lm_calculate returned %d\n", err);
	    }
	}
    }

    pprintf(prn, _("Tolerance = %g\n"), spec->tol);

    if (!err) {
	if (spec->ci == NLS) {
	    if (spec->opt & OPT_C) {
		/* coefficients only: don't bother with GNR */
		add_nls_coeffs(&nlmod, spec);
	    } else {
		/* Use Gauss-Newton Regression for covariance matrix,
		   standard errors */
		nlmod = GNR(spec, (const double **) *spec->Z, 
			    spec->dinfo, prn);
	    }
	} else {
	    /* MLE, GMM */
	    make_nl_model(&nlmod, spec, pdinfo, opt);
	}
    } else if (nlmod.errcode == 0) { 
	/* error code missing: supply one */
	if (spec->generr != 0) {
	    nlmod.errcode = spec->generr;
	} else {
	    nlmod.errcode = E_NOCONV;
	}
    }

    /* ensure that the canonical parameter values get back 
       into external scalars or vectors */
    update_coeff_values(spec->coeff, spec);

 bailout:

    if (spec->nvec > 0 && analytic_mode(spec)) {
	destroy_private_matrices();
    }

    destroy_private_scalars();

    if (spec->Z != pZ) {
	if (!nlmod.errcode) {
	    nls_model_fix_sample(&nlmod, spec, pdinfo);
	}
	destroy_dataset(*spec->Z, spec->dinfo);
	free(spec->Z);
    }

    clear_nlspec(spec);

#if NLS_DEBUG
    fprintf(stderr, "real_nl_model: finishing, dropping %d series\n",
	    pdinfo->v - origv);
#endif

    dataset_drop_last_variables(pdinfo->v - origv, pZ, pdinfo);

    if (!(opt & OPT_A)) {
	if (nlmod.errcode) {
	    if (opt & OPT_U) {
		nlmod.opt |= OPT_U; /* continue on error */
	    }
	} else {
	    set_model_id(&nlmod);
	}
    }

    return nlmod;
}

/**
 * nl_model:
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 * @opt: may include %OPT_V for verbose output, %OPT_R
 * for robust covariance matrix.
 * @prn: printing struct.
 *
 * Computes estimates of a model via nonlinear least squares,
 * maximum likelihood, or GMM.  The model must have been specified 
 * previously, via calls to the function #nl_parse_line.  Those
 * calls determine, among other things, which estimator will
 * be used.
 *
 * Returns: a model struct containing the parameter estimates
 * and associated statistics.
 */

MODEL nl_model (double ***pZ, DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    return real_nl_model(NULL, pZ, pdinfo, opt, prn);
}

/**
 * model_from_nlspec:
 * @spec: nonlinear model specification.
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 * @opt: may include %OPT_V for verbose output, %OPT_A to
 * treat as an auxiliary model, %OPT_C to produce coefficient
 * estimates only (don't bother with GNR to produce standard
 * errors).
 * @prn: printing struct.
 *
 * Computes estimates of the model specified in @spec.
 * The @spec must first be obtained using nlspec_new(), and
 * initialized using nlspec_set_regression_function().  If analytical
 * derivatives are to be used (which is optional but recommended)
 * these are set using nlspec_add_param_with_deriv().
 *
 * Returns: a model struct containing the parameter estimates
 * and associated statistics.
 */

MODEL model_from_nlspec (nlspec *spec, double ***pZ, DATAINFO *pdinfo, 
			 gretlopt opt, PRN *prn)
{
    return real_nl_model(spec, pZ, pdinfo, opt, prn);
}

/**
 * nlspec_new:
 * @ci: %NLS, %MLE or %GMM.
 * @pdinfo: information on dataset.
 *
 * Returns: a pointer to a newly allocated nonlinear model
 * specification, or %NULL on failure.
 */

nlspec *nlspec_new (int ci, const DATAINFO *pdinfo)
{
    nlspec *spec;

    spec = malloc(sizeof *spec);
    if (spec == NULL) {
	return NULL;
    }

    spec->nlfunc = NULL;

    spec->params = NULL;
    spec->nparam = 0;

    spec->aux = NULL;
    spec->naux = 0;
    
    spec->genrs = NULL;
    spec->ngenrs = 0;
    spec->generr = 0;

    spec->fvec = NULL;
    spec->jac = NULL;
    
    spec->coeff = NULL;
    spec->ncoeff = 0;
    spec->nvec = 0;

    spec->hessvec = NULL;

    spec->ci = ci;
    spec->flags = 0;
    spec->opt = OPT_NONE;

    spec->dv = 0;
    *spec->lhname = '\0';
    spec->lhv = 0;
    spec->lvec = NULL;

    spec->iters = 0;
    spec->fncount = 0;
    spec->grcount = 0;

    spec->t1 = spec->real_t1 = pdinfo->t1;
    spec->t2 = spec->real_t2 = pdinfo->t2;
    spec->nobs = spec->t2 - spec->t1 + 1;

    spec->Z = NULL;
    spec->dinfo = NULL;
    spec->prn = NULL;

    spec->oc = NULL;
    spec->missmask = NULL;

    return spec;
}

/**
 * nlspec_destroy:
 * @spec: pointer to nls specification.
 *
 * Frees all resources associated with @spec, and frees the
 * pointer itself.
 */

void nlspec_destroy (nlspec *spec)
{
    clear_nlspec(spec);
    free(spec);
}

/* below: apparatus for estimating a "tsls" model via GMM */

#define IVREG_RESIDNAME  "gmm___e"
#define IVREG_WEIGHTNAME "gmm___V"

static int ivreg_nlfunc_setup (nlspec *spec, MODEL *pmod,
			       DATAINFO *pdinfo)
{
    PRN *prn;
    int v = pmod->list[1];
    int k = pmod->ncoeff;
    const char *buf;
    int i, err = 0;

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
    if (err) {
	return err;
    }

    pprintf(prn, "gmm %s = %s-(", IVREG_RESIDNAME, pdinfo->varname[v]);
    for (i=0; i<k; i++) {
	if (i == 0 && pmod->ifc) {
	    pputs(prn, "b0");
	} else {
	    v = pmod->list[i+2];
	    pprintf(prn, "+b%d*%s", i, pdinfo->varname[v]);
	}
    }
    pputc(prn, ')');

    buf = gretl_print_get_buffer(prn);
    err = nlspec_set_regression_function(spec, buf, pdinfo);
    gretl_print_destroy(prn);

    return err;
}

static int ivreg_oc_setup (nlspec *spec, const int *ilist, 
			   MODEL *pmod, double ***pZ, 
			   DATAINFO *pdinfo, int *rv)
{
    int v = pdinfo->v;
    int i, err;

    /* add GMM residual */
    err = dataset_add_series(1, pZ, pdinfo);

    if (!err) {
	*rv = v;
	strcpy(pdinfo->varname[v], IVREG_RESIDNAME);
	for (i=0; i<pdinfo->n; i++) {
	    (*pZ)[v][i] = pmod->uhat[i];
	}
	err = nlspec_add_ivreg_oc(spec, v, ilist, (const double **) *pZ);
    }
    
    return err;
}

static int ivreg_weights_setup (nlspec *spec, const int *ilist,
				gretlopt opt)
{
    const char *mname = NULL;
    int k = ilist[0];
    gretl_matrix *V;
    int err;

    if (opt & OPT_W) {
	mname = get_optval_string(IVREG, OPT_W);
    }

    if (mname != NULL) {
	/* the user requested a specific weights matrix */
	V = get_matrix_by_name(mname);
	if (V == NULL) {
	    err = E_UNKVAR;
	} else if (V->rows != k || V->cols != k) {
	    err = E_NONCONF;
	} else {
	    err = nlspec_add_weights(spec, mname);
	}
    } else {
	/* use generic weights */
	V = gretl_identity_matrix_new(k);
	if (V == NULL) {
	    return E_ALLOC;
	}
	err = private_matrix_add(V, IVREG_WEIGHTNAME);
	if (err) {
	    gretl_matrix_free(V);
	} else {
	    err = nlspec_add_weights(spec, IVREG_WEIGHTNAME);
	}
    }

    return err;
}

static int ivreg_set_params (nlspec *spec, MODEL *pmod, 
			     double ***pZ, DATAINFO *pdinfo)
{
    int np = pmod->ncoeff;
    char **names;
    int i, err;

    names = strings_array_new_with_length(np, 8);
    if (names == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<np; i++) {
	sprintf(names[i], "b%d", i);
    }

    err = nlspec_add_param_list(spec, np, pmod->coeff, names, pZ, pdinfo);

    free_strings_array(names, np);

    return err;
}

static int finalize_ivreg_model (MODEL *pmod, MODEL *ols, 
				 const int *biglist,
				 const int *mlist,
				 int **ilist,
				 const double **Z,
				 int rv)
{
    int *endolist;
    int yno = mlist[1];
    int t, err = 0;

    pmod->ci = IVREG;
    pmod->opt = OPT_G;

    /* attach tsls-style list */
    pmod->list = gretl_list_copy(biglist);
    if (pmod->list == NULL) {
	return E_ALLOC;
    }

    /* get dep. var. info from OLS model */
    pmod->ybar = ols->ybar;
    pmod->sdy = ols->sdy;

    /* and steal uhat, yhat */
    pmod->uhat = ols->uhat;
    pmod->yhat = ols->yhat;
    ols->uhat = ols->yhat = NULL;

    /* insert residuals and fitted */
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->uhat[t] = Z[rv][t];
	pmod->yhat[t] = Z[yno][t] - Z[rv][t];
    }
    
    endolist = tsls_make_endolist(mlist, ilist, NULL, &err);

    if (!err) {
	gretl_model_set_list_as_data(pmod, "endolist", endolist);
    }

    return err;
}

MODEL ivreg_via_gmm (const int *list, double ***pZ,
		     DATAINFO *pdinfo, gretlopt opt)
{
    int orig_v = pdinfo->v;
    int *mlist = NULL, *ilist = NULL;
    nlspec *spec = NULL;
    MODEL olsmod, model;
    int rv = 0;
    int err = 0;

    gretl_model_init(&olsmod);
    gretl_model_init(&model);

    err = ivreg_process_lists(list, &mlist, &ilist);

    if (!err) {
	spec = nlspec_new(GMM, pdinfo);
	if (spec == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* OLS baseline and check */
	olsmod = lsq(mlist, pZ, pdinfo, OLS, OPT_A | OPT_M | OPT_Z);
	err = olsmod.errcode;
    }

    if (!err) {
	/* add the definition of the GMM residual */
	err = ivreg_nlfunc_setup(spec, &olsmod, pdinfo);
    }    

    if (!err) {
	/* add the orthogonality conditions */
	err = ivreg_oc_setup(spec, ilist, &olsmod, pZ, pdinfo, &rv);
    }

    if (!err) {
	/* add the weights matrix */
	err = ivreg_weights_setup(spec, ilist, opt);
    }

    set_auxiliary_scalars();

    if (!err) {
	/* add the scalar parameters */
	err = ivreg_set_params(spec, &olsmod, pZ, pdinfo);
    }

    if (!err) {
	/* call the GMM estimation routine */
	model = model_from_nlspec(spec, pZ, pdinfo, opt | OPT_G, NULL);
    }

    unset_auxiliary_scalars();

    if (err) {
	model.errcode = err;
    } else {
	/* turn the output model into an "ivreg" type */
	finalize_ivreg_model(&model, &olsmod, list, mlist, &ilist,
			     (const double **) *pZ, rv);
    }

    nlspec_destroy(spec);
    free(mlist);
    free(ilist);
    clear_model(&olsmod);

    dataset_drop_last_variables(pdinfo->v - orig_v, pZ, pdinfo);
    user_matrix_destroy_by_name(IVREG_WEIGHTNAME, NULL);

    return model;
}

/* apparatus for bootstrapping NLS forecast errors */

/* reconstitute the nlspec based on the information saved in @pmod */

static int set_nlspec_from_model (const MODEL *pmod, const double **Z,
				  const DATAINFO *pdinfo)
{
    char line[MAXLEN];
    const char *buf;
    int err = 0;

    buf = gretl_model_get_data(pmod, "nlinfo");
    if (buf == NULL) {
	return E_DATA;
    }

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf) && !err) {
	tailstrip(line);
	if (!strcmp(line, "end nls")) {
	    break;
	}
	err = nl_parse_line(NLS, line, Z, pdinfo, NULL);
    }

    bufgets_finalize(buf);

    if (err) {
	clear_nlspec(&private_spec);
    } 

    return err;
}

/* Given an NLS model @pmod, fill out the forecast error series
   @fcerr from @ft1 to @ft2 using the bootstrap, based on
   resampling the original residuals. See also forecast.c,
   from where this function is called.  Here we concentrate
   solely on parameter uncertainty.  At present we do this
   only for static forecasts.
*/

int nls_boot_calc (const MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		   int ft1, int ft2, double *fcerr) 
{
    nlspec *spec;
    gretl_matrix *fcmat = NULL;
    double *orig_y = NULL;
    double *resu = NULL;
    int origv = pdinfo->v;
    int yno = pmod->list[1];
    int iters = 100; /* just testing */
    int i, s, t, fT;
    int err = 0;

    /* build the 'private' spec, based on pmod */
    err = set_nlspec_from_model(pmod, (const double **) *pZ, pdinfo);
    if (err) {
	return err;
    }

    spec = &private_spec;

    spec->t1 = spec->real_t1 = pmod->t1;
    spec->t2 = spec->real_t2 = pmod->t2;
    spec->nobs = spec->t2 - spec->t1 + 1;

    /* back up the original dependent variable */
    orig_y = copyvec((*pZ)[yno], pdinfo->n);

    /* allocate storage for resampled residuals */
    resu = malloc(pdinfo->n * sizeof *resu);

    /* allocate forecast matrix */
    fT = ft2 - ft1 + 1;
    fcmat = gretl_matrix_alloc(iters, fT);

    if (orig_y == NULL || resu == NULL || fcmat == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* allocate arrays to be passed to minpack */
    spec->fvec = malloc(spec->nobs * sizeof *spec->fvec);
    spec->jac = malloc(spec->nobs * spec->ncoeff * sizeof *spec->jac);

    if (spec->fvec == NULL || spec->jac == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    spec->Z = pZ;
    spec->dinfo = pdinfo;

    for (i=0; i<iters && !err; i++) {

	/* resample from pmod->uhat into resu */
	resample_series(pmod->uhat, resu, pdinfo);

	/* construct y^* = \hat{y} + u^* */
	for (t=spec->t1; t<=spec->t2; t++) {
	    /* FIXME rescale the residuals */
	    (*pZ)[yno][t] = pmod->yhat[t] + resu[t];
	}

#if 0
	for (t=spec->t1; t<=spec->t2; t++) {
	    fprintf(stderr, "%d: y = %g, y* = %g\n", 
		    t, orig_y[t], (*pZ)[yno][t]);
	}
#endif	

	for (t=spec->t1, s=0; t<=spec->t2; t++, s++) {
	    spec->fvec[s] = resu[s];
	}

	/* re-estimate the model on the resampled data */
	if (numeric_mode(spec)) {
	    err = lm_approximate(spec, NULL);
	} else {
	    err = lm_calculate(spec, NULL);
	}

	if (!err) {
	    /* generate and record the forecast */
	    err = generate(pmod->depvar, pZ, pdinfo, OPT_P, NULL);
	    for (t=ft1, s=0; t<=ft2; t++, s++) {
		gretl_matrix_set(fcmat, i, s, (*pZ)[yno][t]);
	    }
	}
    }

    if (!err) {
	/* compute and record the standard deviations of the forecasts:
	   should we add the square root of the residual variance here? 
	*/
	gretl_matrix *sd;

#if 0
	gretl_matrix_print(fcmat, "Forecast matrix");
#endif
	sd = gretl_matrix_column_sd(fcmat, &err);
	if (!err) {
	    double cfac = sqrt((double) iters / (iters - 1));

	    for (t=ft1, s=0; t<=ft2; t++, s++) {
		fcerr[t] = sd->val[s] * cfac;
	    }
	    gretl_matrix_free(sd);
	}
    }

    if (spec->nvec > 0 && analytic_mode(spec)) {
	destroy_private_matrices();
    }

 bailout:

    clear_nlspec(spec);

    dataset_drop_last_variables(pdinfo->v - origv, pZ, pdinfo);

    /* restore original y */
    if (orig_y != NULL) {
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[yno][t] = orig_y[t];
	}
    }

    free(orig_y);
    free(resu);
    gretl_matrix_free(fcmat);

    return err;
}

