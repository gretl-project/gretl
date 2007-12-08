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

/* Nonlinear least squares for libgretl, using minpack; also
   Maximum Likelihood and GMM estimation using BFGS.  Much of
   the GMM-specific material is in the companion file, gmm.c
*/

#include "libgretl.h" 
#include "libset.h"
#include "usermat.h"
#include "gretl_func.h"
#include "nlspec.h"
#include "cmd_private.h"

#include "../../minpack/minpack.h"  

#define NLS_DEBUG 0
#define ML_DEBUG 0
#define RSTEPS 4

enum {
    NUMERIC_DERIVS,
    ANALYTIC_DERIVS
} nls_modes;

struct parm_ {
    char name[VNAMELEN];  /* name of parameter (scalar or vector) */
    char *deriv;          /* string representation of derivative of regression
			     function with respect to param (or NULL) */
    int vnum;             /* ID number of scalar variable in dataset */
    int dnum;             /* ID number of variable holding the derivative */
    int nc;               /* number of individual coefficients associated
			     with the parameter */
    gretl_matrix *vec;    /* pointer to vector parameter */
    gretl_matrix *dmat;   /* pointer to matrix derivative */
};

#define scalar_param(s,i) (s->params[i].vnum > 0)
#define matrix_deriv(s,i) (s->params[i].dmat != NULL)

/* file-scope global variables */

static nlspec private_spec;
static integer one = 1;

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
    np = (s->mode == ANALYTIC_DERIVS)? s->nparam : 0;
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

#if 0 
    /* namespace problems when two or more non-nested functions
       are called in an NLS/MLE/GMM block! */
    set_drop_function_vars(0);
#endif

    j = 0;

    for (i=0; i<ngen && !err; i++) {
	if (i < s->naux) {
	    /* auxiliary variables */
	    strcpy(formula, s->aux[i]);
	} else if (i == s->naux) {
	    /* residual or likelihood function */
	    if (*s->lhname != '\0') {
		sprintf(formula, "%s = %s", s->lhname, s->nlfunc);
	    } else {
		sprintf(formula, "$nl_y = %s", s->nlfunc);
	    }
	} else {
	    /* derivatives/gradients */
	    if (scalar_param(s, j)) {
		sprintf(formula, "$nl_x%d = %s", i, s->params[j].deriv);
	    } else {
		sprintf(formula, "matrix $nl_x%d = %s", i, 
			s->params[j].deriv);
	    }
	    j++;
	}
	
	genrs[i] = genr_compile(formula, s->Z, s->dinfo, &err);

	if (err) {
	    fprintf(stderr, "genr_compile: genrs[%d] = %p, err = %d\n", i, 
		    (void *) genrs[i], err);
	    fprintf(stderr, "formula: '%s'\n", formula);
	    break;
	}

	/* see if the formula actually works, and flush out NAs
	   while we're at it
	*/
	genr_set_na_check(genrs[i]);
	err = execute_genr(genrs[i], s->Z, s->dinfo, s->prn);
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
	    if (genr_is_print(genrs[i])) {
		continue;
	    } else {
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

    if (err) {
	set_drop_function_vars(1);
    }

    return err;
}

/* if i == 0 we're calculating the function; if i > 0 we're calculating
   a derivative */

static int nls_auto_genr (nlspec *s, int i)
{
    int j;

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
	s->generr = execute_genr(s->genrs[j], s->Z, s->dinfo, s->prn);
    }

    if (i == 0 && s->nlfunc == NULL) {
	/* we're done */
	return s->generr;
    }

    j = s->naux + i;
#if NLS_DEBUG
    fprintf(stderr, " j = naux + i = %d+%d = %d: executing genr[%d]\n", 
	    s->naux, i, j, j);
#endif
    s->generr = execute_genr(s->genrs[j], s->Z, s->dinfo, s->prn);

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

static int nlspec_push_param (nlspec *s, const char *name,
			      int v, const double **Z,
			      char *deriv)
{
    parm *params;
    int np = s->nparam;
    int err;
    
    params = realloc(s->params, (np + 1) * sizeof *params);
    if (params == NULL) {
	return E_ALLOC;
    }

    params[np].name[0] = '\0';
    strncat(params[np].name, name, VNAMELEN - 1);
    params[np].deriv = deriv;
    params[np].vnum = v;
    params[np].dnum = 0;
    params[np].nc = 0;
    params[np].vec = NULL;
    params[np].dmat = NULL;

#if NLS_DEBUG
    fprintf(stderr, "added param[%d] = '%s' (v = %d)\n", np, params[np].name, v);
#endif

    s->params = params;
    s->nparam = np + 1;

    if (v > 0) {
	s->params[np].nc = 1;
	err = push_scalar_coeff(s, Z[v][0]);
    } else {
	gretl_matrix *m = get_matrix_by_name(name);
	int k = gretl_vector_get_length(m);

#if NLS_DEBUG
	fprintf(stderr, "vector param: m = %p, k = %d\n", (void *) m, k);
#endif

	s->params[np].vec = m;
	s->params[np].nc = k;
	err = push_vec_coeffs(s, m, k);
	if (!err) {
	    s->nvec += 1;
	}
    }

    return err;
}

static int 
check_param_name (const char *name, const DATAINFO *pdinfo, int *pv)
{
    int v = varindex(pdinfo, name);

    if (v == 0) {
	return E_DATA;
    } else if (v >= pdinfo->v) {
	gretl_matrix *m = get_matrix_by_name(name);
	
	if (m != NULL) {
	    if (gretl_vector_get_length(m) > 0) {
		*pv = 0;
	    } else {
		return E_DATATYPE;
	    }
	} else {
	    return E_UNKVAR;
	}
    } else if (!var_is_scalar(pdinfo, v)) {
	return E_DATATYPE;
    } else {
	*pv = v;
    }

    return 0;
}

/* Scrutinize word and see if it's a new scalar that should be added
   to the NLS specification; if so, add it.  This function is used by
   the fallback get_params_from_nlfunc(), which comes into play only
   if no "deriv" or "params" lines have been provided.
*/

static int 
maybe_add_param_to_spec (nlspec *s, const char *word, 
			 const double **Z, const DATAINFO *pdinfo)
{
    int i, v, err;

#if NLS_DEBUG
    fprintf(stderr, "maybe_add_param: looking at '%s'\n", word);
#endif

    if (function_from_string(word) || const_lookup(word)) {
	return 0;
    }

    v = varindex(pdinfo, word);

    if (v < pdinfo->v) {
	if (var_is_series(pdinfo, v)) {
	    /* only scalars are wanted here */
	    return 0;
	}
    } else {
	sprintf(gretl_errmsg, _("Unknown variable '%s'"), word);
	return E_UNKVAR;
    }

    /* if this term is already present in the spec, skip it */
    for (i=0; i<s->nparam; i++) {
	if (!strcmp(word, s->params[i].name)) {
	    return 0;
	}
    }

#if NLS_DEBUG
    fprintf(stderr, "maybe_add_param: adding '%s'\n", word);
#endif

    err = nlspec_push_param(s, word, v, Z, NULL);

    return err;
}

/* Parse NLS function specification string to find names of variables
   that may figure as parameters of the regression function.  We need
   to do this only if we haven't been given derivatives or a "params"
   line.
*/

static int 
get_params_from_nlfunc (nlspec *s, const double **Z,
			const DATAINFO *pdinfo)
{
    const char *f = s->nlfunc;
    const char *p;
    char name[VNAMELEN];
    int n, np = 0;
    int err = 0;

    if (f == NULL) {
	return 0;
    }

#if NLS_DEBUG
    fprintf(stderr, "get_params: looking at '%s'\n", f);
#endif

    while (*f && !err) {
	p = f;
	if (isalpha(*f) && *(f + 1)) { 
	    /* find a variable name */
	    n = gretl_varchar_spn(f);
	    p += n;
	    if (n > VNAMELEN - 1) {
		/* variable name is too long */
		return 1;
	    }
	    *name = 0;
	    strncat(name, f, n);
	    if (np > 0) {
		/* right-hand side term */
		err = maybe_add_param_to_spec(s, name, Z, pdinfo);
	    }
	    np++;
	} else {
	    p++;
	}
	f = p;
    }

    return err;
}

/* For case where analytical derivatives are not given, the user
   may supply a line like:

   params b0 b1 b2 b3

   specifying the parameters to be estimated.  Here we parse
   such a list and add the parameter info to spec.  The terms
   in the list must be pre-existing scalar variables or
   matrices (vectors).
*/

static int 
nlspec_add_params_from_line (nlspec *s, const char *str,
			     const double **Z, const DATAINFO *pdinfo)
{
    int i, nf = count_fields(str);
    int v, err = 0;

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
	}

	if (!err) {
	    err = check_param_name(name, pdinfo, &v);
	}

	if (!err) {
	    err = nlspec_push_param(s, name, v, Z, NULL);
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
 * @list: list of variables by ID number.
 * @Z: data array.
 * @pdinfo: information on dataset.
 *
 * Adds to @spec a list of (scalar) parameters to be estimated, as
 * given in @list.  For an example of use see arma.c in the
 * gretl plugin directory.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int nlspec_add_param_list (nlspec *spec, const int *list,
			   const double **Z, const DATAINFO *pdinfo)
{
    int i, v, np = list[0];
    int err = 0;

    if (spec->params != NULL || np == 0) {
	return E_DATA;
    }

    for (i=0; i<np && !err; i++) {
	v = list[i+1];
	if (v > 0 && v < pdinfo->v && var_is_scalar(pdinfo, v)) {
	    err = nlspec_push_param(spec, pdinfo->varname[v], v, Z, NULL);
	} else {
	    err = E_DATA;
	}
    }

    if (err) {
	nlspec_destroy_arrays(spec);
    } 

    return err;
}

/* return the address in the dataset, or in a user matrix, that
   holds the "external" value of a given coefficient */

static double *coeff_address (nlspec *s, int i)
{
    double **Z = *(s->Z);
    int j, k, pos = 0;

    if (s->nvec == 0) {
	/* the mapping is simple */
	return &(Z[s->params[i].vnum][0]);
    }

    for (j=0; j<s->nparam; j++) {
	if (scalar_param(s, j)) {
	    if (i == pos) {
		return &(Z[s->params[j].vnum][0]);
	    }
	    pos++;
	} else {
	    gretl_matrix *m = get_matrix_by_name(s->params[j].name);

	    if (m != s->params[j].vec) {
		fprintf(stderr, "*** coeff_address: by name, '%s' is at %p; "
			"stored addr = %p\n", s->params[j].name,
			(void *) m, (void *) s->params[j].vec);
		s->params[j].vec = m;
	    }

	    k = s->params[j].nc;
	    if (i >= pos && i < pos + k) {
		return &(s->params[j].vec->val[i - pos]);
	    }
	    pos += k;
	}
    }

    fprintf(stderr, "Couldn't find location for coeff %d\n", i);
	
    return NULL;
}

int update_coeff_values (const double *x, nlspec *s)
{
    double *d;
    int i;

    /* write the values produced by the optimizer into the dataset */

    for (i=0; i<s->ncoeff; i++) {
	d = coeff_address(s, i);
	if (d == NULL) {
	    return 1;
	}
	*d = x[i];
#if NLS_DEBUG
	fprintf(stderr, "update params: revised coeff[%d] = %.14g\n", i, x[i]);
#endif
    }

    return 0;
}

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

    if (t2 - t1 + 1 < s->ncoeff) {
	return E_DF;
    }

    for (t=t1; t<=t2; t++) {
	if (na((*s->Z)[v][t])) {
	    fprintf(stderr, "  after setting t1=%d, t2=%d, "
		    "got NA for var %d at obs %d\n", t1, t2, v, t);
	    return E_MISSDATA;
	}
    }  

 nl_miss_exit:

    s->t1 = t1;
    s->t2 = t2;
    s->nobs = t2 - t1 + 1;

#if NLS_DEBUG
    fprintf(stderr, "  after: spec->t1 = %d, spec->t2 = %d, spec->nobs = %d\n\n",
	    s->t1, s->t2, s->nobs);
#endif

    return 0;
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

/* default numerical calculation of gradient in context of BFGS */

int BFGS_numeric_gradient (double *b, double *g, int n,
			   BFGS_CRIT_FUNC func, void *data)
{
    double bi0, f1, f2;
    gretlopt opt = OPT_NONE;
    int i;

    if (opt == OPT_R) {
	/* Richardson */
	double df[RSTEPS];
	double eps = 1.0e-4;
	double d = 0.0001;
	double v = 2.0;
	double h, p4m;
	int r = RSTEPS;
	int k, m;

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
		    return 1;
		}		    
		df[k] = (f2 - f1) / (2.0 * h); 
		h /= v;
	    }
	    b[i] = bi0;
	    p4m = 4.0;
	    for (m=0; m<r-1; m++) {
		for (k=0; k<r-m; k++) {
		    df[k] = (df[k+1] * p4m - df[k]) / (p4m - 1.0);
		}
		p4m *= 4.0;
	    }
	    g[i] = df[0];
	}
    } else {
	/* simple gradient calculation */
	const double h = 1.0e-8;

	for (i=0; i<n; i++) {
	    bi0 = b[i];
	    b[i] = bi0 - h;
	    f1 = func(b, data);
	    b[i] = bi0 + h;
	    f2 = func(b, data);
	    b[i] = bi0;
	    if (na(f1) || na(f2)) {
		return 1;
	    }
	    g[i] = (f2 - f1) / (2.0 * h);
	}
    }

    return 0;
}

/* this function is used in the context of the minpack callback, and
   also for checking derivatives in the MLE case
*/

static int nl_function_calc (double *f, void *p)
{
    nlspec *s = (nlspec *) p;
    const double *y;
    int t;

#if NLS_DEBUG > 1
    fprintf(stderr, "\n*** nl_function_calc called\n");
#endif

    /* calculate residual given current parameter estimates */
    if (nl_calculate_fvec(s)) {
	return 1;
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
	} else {
	    /* derivative is scalar or series */
	    int s, ser, v = spec->params[j].dnum;

	    if (v == 0) {
		/* FIXME? */
		v = spec->params[j].dnum = spec->dinfo->v - 1;
	    }
	    ser = var_is_series(spec->dinfo, v);

#if NLS_DEBUG
	    fprintf(stderr, "param[%d]: dnum = %d, series = %d\n", j, v, ser);
#endif

	    /* transcribe from dataset to array g */
	    for (t=0; t<T; t++) {
		s = (ser)? (t + spec->t1) : 0;
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
    int t1, t2;
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
	fprintf(stderr, "param[%d]: done nls_calculate_deriv\n", j);
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
	} else {
	    /* derivative may be series or scalar */
	    int v = spec->params[j].dnum;

#if ML_DEBUG > 1
	    fprintf(stderr, "param[%d], dnum = %d\n", j, v);
#endif

	    if (var_is_series(spec->dinfo, v)) {
		t1 = spec->t1;
		t2 = spec->t2;
	    } else {
		t1 = t2 = 0;
	    }

	    g[i] = 0.0;

	    for (t=t1; t<=t2; t++) {
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
   model struct. */

static void add_stats_to_model (MODEL *pmod, nlspec *spec,
				const double **Z)
{
    int dv = spec->dv;
    double d, tss;
    int t;

    pmod->ess = spec->crit;
    pmod->sigma = sqrt(pmod->ess / (pmod->nobs - spec->ncoeff));
    
    pmod->ybar = gretl_mean(pmod->t1, pmod->t2, Z[dv]);
    pmod->sdy = gretl_stddev(pmod->t1, pmod->t2, Z[dv]);

    tss = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	d = Z[dv][t] - pmod->ybar;
	tss += d * d;
    }  

    if (tss == 0.0) {
	pmod->rsq = NADBL;
    } else {
	pmod->rsq = 1.0 - pmod->ess / tss;
    }

    pmod->adjrsq = NADBL;
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

    /* expand Hinv (this would be unnecessary if we had
       a special multiplication routine for vech's) 
    */
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

    if (spec->mode == NUMERIC_DERIVS) {
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
	    } else {
		v = spec->params[i].dnum;
		if (var_is_scalar(spec->dinfo, v)) {
		    x = (*spec->Z)[v][0];
		}
		for (t=0; t<T; t++) {
		    if (var_is_series(spec->dinfo, v)) {
			x = (*spec->Z)[v][t + spec->t1];
		    }
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
	for (i=0; i<k; i++) {
	    for (j=0; j<=i; j++) {
		x = gretl_matrix_get(V, i, j);
		pmod->vcv[ijton(i, j, k)] = x;
	    }
	}
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
	vcvopt = VCV_HESSIAN;
    } else {
	/* either OPG or QML */
	err = mle_build_vcv(pmod, spec, &vcvopt);
    }

    if (!err) {
	for (i=0; i<k; i++) {
	    x = pmod->vcv[ijton(i, i, k)];
	    pmod->sderr[i] = sqrt(x);
	}
	gretl_model_set_int(pmod, "ml_vcv", vcvopt);
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
	pmod->depvar = gretl_strdup(pdinfo->varname[spec->dv]);
	if (pmod->depvar == NULL) {
	    err = E_ALLOC;
	}
    } else if (spec->ci == MLE || (spec->ci == GMM && spec->nlfunc != NULL)) {
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

    return 0;
}

static void 
add_fit_resid_to_model (MODEL *pmod, nlspec *spec, double *uhat, 
			const double **Z, int perfect)
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
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR);
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

static MODEL GNR (double *uhat, double *jac, nlspec *spec,
		  double ***pZ, const DATAINFO *pdinfo, 
		  PRN *prn)
{
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

    if (gretl_iszero(0, spec->nobs - 1, uhat)) {
	pputs(prn, _("Perfect fit achieved\n"));
	perfect = 1;
	for (t=0; t<spec->nobs; t++) {
	    uhat[t] = 1.0;
	}
	spec->crit = 0.0;
    }

    /* number of variables = 1 (const) + 1 (depvar) + spec->ncoeff
       (derivatives) */
    gdinfo = create_new_dataset(&gZ, spec->ncoeff + 2, pdinfo->n, 0);
    if (gdinfo == NULL) {
	gretl_model_init(&gnr);
	gnr.errcode = E_ALLOC;
	return gnr;
    }

    /* transcribe sample info */
    gdinfo->t1 = spec->t1;
    gdinfo->t2 = spec->t2;

#if 0
    fprintf(stderr, "pdinfo->n = %d, gdinfo->t1 = %d, gdinfo->t2 = %d\n",
	    pdinfo->n, gdinfo->t1, gdinfo->t2);
#endif
    
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

    if (spec->mode == ANALYTIC_DERIVS) {
	for (i=0; i<spec->ncoeff; i++) {
	    v = i + 2;
	    for (t=0; t<gdinfo->t1; t++) {
		gZ[v][t] = NADBL;
	    }
	    for (t=gdinfo->t2; t<gdinfo->n; t++) {
		gZ[v][t] = NADBL;
	    }
	}
	get_nls_derivs(spec->nparam, T, spec->t1, NULL, gZ + 2, spec);
    } else {
	for (i=0; i<spec->ncoeff; i++) {
	    v = i + 2;
	    j = T * i; /* calculate offset into jac */
	    for (t=0; t<gdinfo->n; t++) {
		if (t < gdinfo->t1 || t > gdinfo->t2) {
		    gZ[v][t] = NADBL;
		} else {
		    gZ[v][t] = jac[j++];
		}
	    }
	}
    }

#if NLS_DEBUG > 1
    print_GNR_dataset(glist, gZ, gdinfo);
#endif

    lsqopt = OPT_A;
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
	strcpy(gretl_errmsg, _("Failed to calculate Jacobian"));
	gnr.errcode = E_DATA;
    }

    if (gnr.errcode == 0) {
	gnr.ci = spec->ci;
	add_stats_to_model(&gnr, spec, (const double **) *pZ);
	if (add_nls_std_errs_to_model(&gnr)) {
	    gnr.errcode = E_ALLOC;
	}
    }

    if (gnr.errcode == 0) {
	ls_criteria(&gnr);
	add_coeffs_to_model(&gnr, spec->coeff);
	add_param_names_to_model(&gnr, spec, pdinfo);
	add_fit_resid_to_model(&gnr, spec, uhat, (const double **) *pZ, 
			       perfect);
	gnr.list[1] = spec->dv;

	/* set relevant data on model to be shipped out */
	gretl_model_set_int(&gnr, "iters", iters);
	gretl_model_set_double(&gnr, "tol", spec->tol);
	transcribe_nls_function(&gnr, spec->nlfunc);
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

/* work up the results of estimation into the form of a gretl
   MODEL struct */

static int make_nl_model (MODEL *pmod, nlspec *spec, 
			  const DATAINFO *pdinfo)
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
	mle_criteria(pmod, 0);
    }

    add_coeffs_to_model(pmod, spec->coeff);

    pmod->errcode = add_param_names_to_model(pmod, spec, pdinfo);

    if (!pmod->errcode) {
	if (pmod->ci == MLE) {
	    pmod->errcode = mle_add_vcv(pmod, spec);
	} else {
	    pmod->errcode = gmm_add_vcv(pmod, spec);
	}
    }

    if (!pmod->errcode) {
	gretl_model_set_int(pmod, "fncount", spec->fncount);
	gretl_model_set_int(pmod, "grcount", spec->grcount);
	gretl_model_set_double(pmod, "tol", spec->tol);
    }

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

    spec->mode = NUMERIC_DERIVS;
    spec->opt = OPT_NONE;

    spec->dv = 0;
    spec->lhv = 0;
    spec->lvec = NULL;

    spec->iters = 0;
    spec->fncount = 0;
    spec->grcount = 0;

    spec->t1 = spec->t2 = 0;
    spec->nobs = 0;

    spec->Z = NULL;
    spec->dinfo = NULL;
    spec->prn = NULL;

    if (spec->oc != NULL) {
	oc_set_destroy(spec->oc);
	spec->oc = NULL;
    }
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

/* drivers for BFGS code below */

static int mle_calculate (nlspec *s, double *fvec, double *jac, PRN *prn)
{
    int err = 0;

    if (s->mode == ANALYTIC_DERIVS) {
	integer m = s->nobs;
	integer n = s->ncoeff;
	integer ldjac = m; 

	err = check_derivatives(m, n, s->coeff, fvec, jac, ldjac, 
				prn, s);
    }

    if (!err) {
	BFGS_GRAD_FUNC gradfun = (s->mode == ANALYTIC_DERIVS)?
	    get_mle_gradient : NULL;
	int maxit = libset_get_int(BFGS_MAXITER);

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

/* driver for minpack Levenberg-Marquandt code for use when analytical
   derivatives have been supplied */

static int lm_calculate (nlspec *spec, double *fvec, double *jac, 
			 PRN *prn)
{
    integer info, lwa;
    integer m, n, ldjac;
    integer *ipvt;
    doublereal *wa;
    int err = 0;

    m = spec->nobs;              /* number of observations */
    n = spec->ncoeff;            /* number of coefficients */
    lwa = 5 * n + m;             /* work array size */
    ldjac = m;                   /* leading dimension of jac array */

    wa = malloc(lwa * sizeof *wa);
    ipvt = malloc(n * sizeof *ipvt);

    if (wa == NULL || ipvt == NULL) {
	err = E_ALLOC;
	goto nls_cleanup;
    }

    err = check_derivatives(m, n, spec->coeff, fvec, jac, ldjac, 
			    prn, spec);
    if (err) {
	goto nls_cleanup; 
    }

    /* note: maxfev is automatically set to 100*(n + 1) */

    /* call minpack */
    lmder1_(nls_calc, &m, &n, spec->coeff, fvec, jac, &ldjac, &spec->tol, 
	    &info, ipvt, wa, &lwa, spec);

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
    case 4: /* is this right? */
	pprintf(prn, _("Convergence achieved after %d iterations\n"),
		spec->iters);
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
    /* write current parameter values into dataset Z */
    update_coeff_values(x, p);

    /* calculate function at x, results into fvec */    
    if (nl_function_calc(fvec, p)) {
	*iflag = -1;
    }

    return 0;
}

/* driver for minpack Levenberg-Marquandt code for use when the
   Jacobian must be approximated numerically */

static int 
lm_approximate (nlspec *spec, double *fvec, double *jac, PRN *prn)
{
    integer info, m, n, ldjac;
    integer maxfev, mode = 1, nprint = 0, nfev = 0;
    integer iflag = 0;
    integer *ipvt;
    doublereal gtol = 0.0;
    doublereal epsfcn = 0.0, factor = 100.;
    doublereal *diag, *qtf;
    doublereal *wa1, *wa2, *wa3, *wa4;
    int err = 0;
    
    m = spec->nobs;              /* number of observations */
    n = spec->ncoeff;            /* number of parameters */
    ldjac = m;                   /* leading dimension of jac array */

    maxfev = 200 * (n + 1);      /* max iterations */

    diag = malloc(n * sizeof *diag);
    qtf = malloc(n * sizeof *qtf);
    wa1 = malloc(n * sizeof *wa1);
    wa2 = malloc(n * sizeof *wa2);
    wa3 = malloc(n * sizeof *wa3);
    wa4 = malloc(m * sizeof *wa4);
    ipvt = malloc(n * sizeof *ipvt);

    if (diag == NULL || qtf == NULL ||
	wa1 == NULL || wa2 == NULL || wa3 == NULL || wa4 == NULL ||
	ipvt == NULL) {
	err = E_ALLOC;
	goto nls_cleanup;
    }

    /* call minpack */
    lmdif_(nls_calc_approx, &m, &n, spec->coeff, fvec, 
	   &spec->tol, &spec->tol, &gtol, &maxfev, &epsfcn, diag, 
	   &mode, &factor, &nprint, &info, &nfev, jac, &ldjac, 
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
	pprintf(prn, _("Convergence achieved after %d iterations\n"),
		spec->iters);
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
	fdjac2_(nls_calc_approx, &m, &n, spec->coeff, fvec, jac, 
		&ldjac, &iflag, &epsfcn, wa4, spec);
	spec->crit = ess;
	spec->iters = iters;
	spec->opt = opt;
    }

 nls_cleanup:

    free(diag);
    free(qtf);
    free(wa1);
    free(wa2);
    free(wa3);
    free(wa4);
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
    int v, err = 0;

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

    err = check_param_name(name, pdinfo, &v);
    
    if (!err) {
	err = nlspec_push_param(spec, name, v, Z, deriv);
	if (err) {
	    free(deriv);
	    deriv = NULL;
	}
    }

    free(name);

    if (!err) {
	spec->mode = ANALYTIC_DERIVS;
    }

#if NLS_DEBUG
    if (!err) {
	fprintf(stderr, "add_param_with_deriv: '%s'\n"
		" set vnum = %d, initial value = %.14g\n", dstr, 
		v, Z[v][0]);
    }
#endif

    return err;
}

static void get_aux_command_word (char *word, const char *line)
{
    int n = gretl_varchar_spn(line);

    *word = '\0';

    if (n > 0) {
	if (n > FN_NAMELEN - 1) {
	    n = FN_NAMELEN - 1;
	}
	strncat(word, line, n);
    } 
}

static int screen_bad_aux (const char *line, const DATAINFO *pdinfo)
{
    char word[FN_NAMELEN];
    int ci, err = E_DATA;

    get_aux_command_word(word, line);
    ci = gretl_command_number(word);

    if (ci == GENR) {
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
    char **aux;
    char *this;
    int nx = spec->naux + 1;
    int err = 0;

#if NLS_DEBUG
    fprintf(stderr, "nlspec_add_aux: s = '%s'\n", s);
#endif

    err = screen_bad_aux(s, pdinfo);
    if (err) {
	return err;
    }

    this = gretl_strdup(s);
    if (this == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	aux = realloc(spec->aux, nx * sizeof *spec->aux);
	if (aux == NULL) {
	    free(this);
	    err = E_ALLOC;
	}
    }

    if (!err) {
	spec->aux = aux;
	spec->aux[nx - 1] = this;
	spec->naux += 1;
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
	spec->dv = varindex(pdinfo, vname);
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
 * nls_parse_line:
 * @ci: either %NLS or %MLE (docs not finished on this)
 * @line: specification of regression function or derivative
 *        of this function with respect to a parameter.
 * @Z: data array.
 * @pdinfo: information on dataset.
 * @prn: gretl printing struct (for warning messages).
 *
 * This function is used to create the specification of a
 * nonlinear regression function, to be estimated via #nls.
 * It should first be called with a @line containing a
 * string specification of the regression function.  Optionally,
 * it can then be called one or more times to specify 
 * analytical derivatives for the parameters of the regression
 * function.  
 *
 * The format of @line should be that used in the #genr function.
 * When specifying the regression function, the formula may
 * optionally be preceded by the string %nls.  When specifying
 * a derivative, @line must start with the string %deriv.  See
 * the gretl manual for details.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int nls_parse_line (int ci, const char *line, const double **Z,
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
		if (s->mode != ANALYTIC_DERIVS && s->params != NULL) {
		    strcpy(gretl_errmsg, _("You cannot supply both a \"params\" "
					   "line and analytical derivatives"));
		    err = E_PARSE;
		} else {
		    err = nlspec_add_param_with_deriv(s, line, Z, pdinfo);
		}
	    } else if (!strncmp(line, "params", 6)) {
		if (s->mode != ANALYTIC_DERIVS) {
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

/* static function providing the real content for the two public
   wrapper functions below */

static MODEL real_nls (nlspec *spec, double ***pZ, DATAINFO *pdinfo, 
		       gretlopt opt, PRN *prn)
{
    MODEL nlsmod;
    double *fvec = NULL;
    double *jac = NULL;
    int origv = pdinfo->v;
    int i, t, err = 0;

#if NLS_DEBUG
    fprintf(stderr, "real_nls: starting\n");
#endif

    gretl_model_init(&nlsmod);
    gretl_model_smpl_init(&nlsmod, pdinfo);

    if (spec == NULL) {
	/* we use the static spec composed via nls_parse_line() */
	spec = &private_spec;
    }

    if (spec->nlfunc == NULL && spec->ci != GMM) {
	strcpy(gretl_errmsg, _("No function has been specified"));
	nlsmod.errcode = E_PARSE;
	goto bailout;
    } 

    spec->opt = opt;

    /* make pZ, pdinfo and prn available */
    spec->Z = pZ;
    spec->dinfo = pdinfo;
    spec->prn = prn;

    if (spec->mode == NUMERIC_DERIVS && spec->nparam == 0 &&
	spec->ci != GMM) {
	err = get_params_from_nlfunc(spec, (const double **) *pZ, pdinfo);
	if (err) {
	    if (err == 1) {
		nlsmod.errcode = E_PARSE;
	    } else {
		nlsmod.errcode = err;
	    }
	    goto bailout;
	}
    }

    if (check_spec_requirements(spec)) {
	nlsmod.errcode = E_PARSE;
	goto bailout;
    } 

    if (spec->ci == GMM) {
	err = gmm_missval_check_etc(spec);
    } else {
	err = nl_missval_check(spec);
    }

    if (err) {
	nlsmod.errcode = err;
	goto bailout;
    }

    /* allocate arrays to be passed to minpack */
    fvec = malloc(spec->nobs * sizeof *fvec);
    jac = malloc(spec->nobs * spec->ncoeff * sizeof *jac);

    if (fvec == NULL || jac == NULL) {
	nlsmod.errcode = E_ALLOC;
	goto bailout;
    }

    if (spec->lvec != NULL) {
	for (t=0; t<spec->nobs; t++) {
	    fvec[t] = spec->lvec->val[t];
	}
    } else {
	i = 0;
	for (t=spec->t1; t<=spec->t2; t++) {
	    fvec[i++] = (*spec->Z)[spec->lhv][t];
	}
    }

    /* get tolerance from user setting or default */
    if (USES_BFGS(spec->ci)) {
	spec->tol = libset_get_double(BFGS_TOLER);
    } else {
	spec->tol = libset_get_double(NLS_TOLER);
    }

    pputs(prn, (spec->mode == NUMERIC_DERIVS)?
	  _("Using numerical derivatives\n") :
	  _("Using analytical derivatives\n"));

    if (spec->ci == MLE) {
	err = mle_calculate(spec, fvec, jac, prn);
    } else if (spec->ci == GMM) {
	err = gmm_calculate(spec, fvec, jac, prn);
    } else {
	/* invoke appropriate minpack driver function */
	if (spec->mode == NUMERIC_DERIVS) {
	    err = lm_approximate(spec, fvec, jac, prn);
	} else {
	    err = lm_calculate(spec, fvec, jac, prn);
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
		add_nls_coeffs(&nlsmod, spec);
	    } else {
		/* Use Gauss-Newton Regression for covariance matrix,
		   standard errors */
		nlsmod = GNR(fvec, jac, spec, pZ, pdinfo, prn);
	    }
	} else {
	    make_nl_model(&nlsmod, spec, pdinfo);
	}
    } else {
	if (nlsmod.errcode == 0) { 
	    if (spec->generr != 0) {
		nlsmod.errcode = spec->generr;
	    } else {
		nlsmod.errcode = E_NOCONV;
	    }
	}
    }

 bailout:

    free(fvec);
    free(jac);

    if (spec->nvec > 0 && spec->mode == ANALYTIC_DERIVS) {
	destroy_private_matrices();
    }

    clear_nlspec(spec);

    dataset_drop_last_variables(pdinfo->v - origv, pZ, pdinfo);

    if (nlsmod.errcode == 0 && !(opt & OPT_A)) {
	set_model_id(&nlsmod);
    }

    set_drop_function_vars(1);

    return nlsmod;
}

/**
 * nls:
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 * @opt: may include %OPT_V for verbose output, %OPT_R
 * for robust covariance matrix.
 * @prn: printing struct.
 *
 * Computes estimates of a model via nonlinear least squares.
 * The model must have been specified previously, via calls to
 * the function #nls_parse_line.  
 *
 * Returns: a model struct containing the parameter estimates
 * and associated statistics.
 */

MODEL nls (double ***pZ, DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    return real_nls(NULL, pZ, pdinfo, opt, prn);
}

/**
 * model_from_nlspec:
 * @spec: nls specification.
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 * @opt: may include %OPT_V for verbose output, %OPT_A to
 * treat as an auxiliary model, %OPT_C to produce coefficient
 * estimates only (don't bother with GNR to produce standard
 * errors).
 * @prn: printing struct.
 *
 * Computes estimates of the model specified in @spec, via nonlinear 
 * least squares. The @spec must first be obtained using nlspec_new(), and
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
    return real_nls(spec, pZ, pdinfo, opt, prn);
}

/**
 * nlspec_new:
 * @ci: %NLS, %MLE or %GMM.
 * @pdinfo: information on dataset.
 *
 * Returns: a pointer to a newly allocated nls specification,
 * or %NULL on failure.
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
    
    spec->coeff = NULL;
    spec->ncoeff = 0;
    spec->nvec = 0;

    spec->hessvec = NULL;

    spec->ci = ci;
    spec->mode = NUMERIC_DERIVS;
    spec->opt = OPT_NONE;

    spec->dv = 0;
    *spec->lhname = '\0';
    spec->lhv = 0;
    spec->lvec = NULL;

    spec->iters = 0;
    spec->fncount = 0;
    spec->grcount = 0;

    spec->t1 = pdinfo->t1;
    spec->t2 = pdinfo->t2;
    spec->nobs = spec->t2 - spec->t1 + 1;

    spec->Z = NULL;
    spec->dinfo = NULL;
    spec->prn = NULL;

    spec->oc = NULL;

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

static void free_triangular_array (double **m, int n)
{
    int i;

    if (m != NULL) {
	for (i=0; i<n; i++) {
	    free(m[i]);
	}
	free(m);
    }
}

static double **triangular_array_new (int n)
{
    double **m;
    int i;

    m = malloc(n * sizeof *m);

    if (m != NULL) {
	for (i=0; i<n; i++) {
	    m[i] = NULL;
	}
	for (i=0; i<n; i++) {
	    m[i] = malloc((i + 1) * sizeof **m);
	    if (m[i] == NULL) {
		free_triangular_array(m, n);
		m = NULL;
		break;
	    }
	}
    }

    return m;
}

#define BFGS_DEBUG 0

/* apparatus for constructing numerical approximation to
   the Hessian */

static void hess_h_init (double *h, double *h0, int n)
{
    int i;

    for (i=0; i<n; i++) {
	h[i] = h0[i];
    }
}

static void hess_h_reduce (double *h, double v, int n)
{
    int i;

    for (i=0; i<n; i++) {
	h[i] /= v;
    }
}

static void hess_b_adjust_i (double *c, double *b, double *h, int n, 
			     int i, double sgn)
{
    int k;

    for (k=0; k<n; k++) {
	c[k] = b[k] + (k == i) * sgn * h[i];
    }
}

static void hess_b_adjust_ij (double *c, double *b, double *h, int n, 
			      int i, int j, double sgn)
{
    int k;

    for (k=0; k<n; k++) {
	c[k] = b[k] + (k == i) * sgn * h[i] +
	    (k == j) * sgn * h[j];
    }
}

/* The algorithm below implements the method of Richardson
   Extrapolation.  It is derived from code in the gnu R package
   "numDeriv" by Paul Gilbert, which was in turn derived from C code
   by Xinqiao Liu.  Turned back into C and modified for gretl by
   Allin Cottrell, June 2006.
*/

double *numerical_hessian (double *b, int n, BFGS_CRIT_FUNC func, void *data,
			   int *err)
{
    double Dx[RSTEPS];
    double Hx[RSTEPS];
    double *c = NULL;
    double *D = NULL;
    double *h0 = NULL;
    double *h = NULL;
    double *Hd = NULL;

    gretl_matrix *V = NULL;
    double *vcv = NULL;

    /* numerical parameters */
    int r = RSTEPS;      /* number of Richardson steps */
    double eps = 1.0e-4;
    double d = 0.0001;
    double v = 2.0;      /* reduction factor for h */

    double f0, f1, f2;
    double p4m;

    int vn = (n * (n + 1)) / 2;
    int dn = vn + n;
    int i, j, k, m, u;

    c  = malloc(n * sizeof *c);
    h0 = malloc(n * sizeof *h0);
    h  = malloc(n * sizeof *h);
    Hd = malloc(n * sizeof *Hd);
    D  = malloc(dn * sizeof *D);

    if (c == NULL || h0 == NULL || h == NULL || 
	Hd == NULL || D == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    /* vech form of variance matrix */
    V = gretl_column_vector_alloc(vn);
    if (V == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }	

    for (i=0; i<n; i++) {
	h0[i] = (fabs(b[i]) < 0.01)? eps : d * b[i];
    }

    f0 = func(b, data);

    /* first derivatives and Hessian diagonal */

    for (i=0; i<n; i++) {
	hess_h_init(h, h0, n);
	for (k=0; k<r; k++) {
	    hess_b_adjust_i(c, b, h, n, i, 1);
	    f1 = func(c, data);
	    if (na(f1)) {
		*err = E_NAN;
		goto bailout;
	    }
	    hess_b_adjust_i(c, b, h, n, i, -1);
	    f2 = func(c, data);
	    if (na(f2)) {
		*err = E_NAN;
		goto bailout;
	    }
	    /* F'(i) */
	    Dx[k] = (f1 - f2) / (2.0 * h[i]); 
	    /* F''(i) */
	    Hx[k] = (f1 - 2.0*f0 + f2) / (h[i] * h[i]);
	    hess_h_reduce(h, v, n);
	}
	p4m = 4;
	for (m=0; m<r-1; m++) {
	    for (k=0; k<r-m; k++) {
		Dx[k] = (Dx[k+1] * p4m - Dx[k]) / (p4m - 1);
		Hx[k] = (Hx[k+1] * p4m - Hx[k]) / (p4m - 1);
	    }
	    p4m *= 4;
	}
	D[i] = Dx[0];
	Hd[i] = Hx[0];
    }

    /* second derivatives: lower half of Hessian only */

    u = n;
    for (i=0; i<n; i++) {
	for (j=0; j<=i; j++) {
	    if (i == j) {
		D[u] = Hd[i];
	    } else {
		hess_h_init(h, h0, n);
		for (k=0; k<r; k++) {
		    hess_b_adjust_ij(c, b, h, n, i, j, 1);
		    f1 = func(c, data);
		    if (na(f1)) {
			*err = E_NAN;
			goto bailout;
		    }
		    hess_b_adjust_ij(c, b, h, n, i, j, -1);
		    f2 = func(c, data);
		    if (na(f2)) {
			*err = E_NAN;
			goto bailout;
		    }
		    /* cross-partial */
		    Dx[k] = (f1 - 2.0*f0 + f2 - Hd[i]*h[i]*h[i]
			     - Hd[j]*h[j]*h[j]) / (2.0*h[i]*h[j]);
		    hess_h_reduce(h, v, n);
		}
		p4m = 4.0;
		for (m=0; m<r-1; m++) {
		    for (k=0; k<r-m; k++) {
			Dx[k] = (Dx[k+1] * p4m - Dx[k]) / (p4m - 1);
		    }
		    p4m *= 4.0;
		}
		D[u] = Dx[0];
	    }
	    u++;
	}
    }

    /* transcribe the negative of the Hessian */
    u = n;
    for (i=0; i<n; i++) {
	for (j=0; j<=i; j++) {
	    k = ijton(i, j, n);
	    V->val[k] = -D[u++];
	}
    }

    *err = gretl_invert_packed_symmetric_matrix(V);
    if (!*err) {
	vcv = gretl_matrix_steal_data(V);
    } else {
	fprintf(stderr, "numerical hessian: failed to invert V\n");
	gretl_packed_matrix_print(V, "V");
    }

    gretl_matrix_free(V);

 bailout:

    if (*err == E_NAN) {
	fprintf(stderr, "Got E_NAN in numerical_hessian()\n");
    }

    free(c);
    free(D);
    free(h0);
    free(h);
    free(Hd);

    return vcv;
}

#define stepfrac	0.2
#define acctol		1.0e-7 /* alt: 0.0001 or 1.0e-7 (?) */
#define reltest		10.0

#define coeff_unchanged(a,b) (reltest + a == reltest + b)

static void reverse_gradient (double *g, int n)
{
    int i;

    for (i=0; i<n; i++) {
	g[i] = -g[i];
    }
}

static int broken_gradient (double *g, int n)
{
    int i;

    for (i=0; i<n; i++) {
	if (isnan(g[i])) {
	    return 1;
	}
    }

    return 0;
}

/* 
   If "set inivals" has been used, supersede whatever initial
   values were there by those given by the user (the customer is
   always right).  In addition, respect user settings for the
   maximum number of iterations and the convergence tolerance.
*/

static void BFGS_get_user_values (double *b, int n, int *maxit,
				  double *reltol, gretlopt opt,
				  PRN *prn)
{
    const gretl_matrix *uinit;
    int uilen, umaxit;
    double utol;
    int i;

    uinit = get_init_vals();
    uilen = gretl_vector_get_length(uinit);

    if (uilen > 0) {
	/* the user has given something */
	if (uilen < n) {
	    fprintf(stderr, "Only %d initial values given, but %d "
		    "are necessary\n", uilen, n);
	} else {
	    for (i=0; i<n; i++) {
		b[i] = uinit->val[i];
	    }
	    if (opt & OPT_V) {
		pputs(prn, _("\n\n*** User-specified starting values:\n"));
		for (i=0; i<n; i++) {
		    pprintf(prn, " %12.6f", b[i]);
		    if (i%6 == 5) {
			pputc(prn, '\n');
		    }
		}
		pputs(prn, "\n\n");
	    }
	    free_init_vals();
	}
    }

    umaxit = libset_get_int(BFGS_MAXITER);
    if (umaxit > 0) {
	*maxit = umaxit;
    }
    
    utol = libset_get_double(BFGS_TOLER);
    if (utol != get_default_nls_toler()) {
	/* it has actually been set */
	*reltol = utol;
	if (!(opt & OPT_Q)) {
	    fprintf(stderr, "user-specified BFGS tolerance = %g\n", utol);
	}
    }	
}

/**
 * BFGS_max:
 * @b: array of adjustable coefficients.
 * @n: number elements in array @b.
 * @maxit: the maximum number of iterations to allow.
 * @reltol: relative tolerance for terminating iteration.
 * @fncount: location to receive count of function evaluations.
 * @grcount: location to receive count of gradient evaluations.
 * @cfunc: pointer to function used to calculate maximand.
 * @crittype: code for type of the maximand/minimand: should
 * be %C_LOGLIK, %C_GMM or %C_OTHER.  Used only in printing
 * iteration info.
 * @gradfunc: pointer to function used to calculate the 
 * gradient, or %NULL for default numerical calculation.
 * @data: pointer that will be passed as the last
 * parameter to the callback functions @cfunc and @gradfunc.
 * @opt: may contain %OPT_V for verbose operation.
 * @prn: printing struct (or %NULL).  Only used if @opt
 * includes %OPT_V.
 *
 * Obtains the set of values for @b which jointly maximize the
 * criterion value as calculated by @cfunc.  Uses the BFGS
 * variable-metric method.  Based on Pascal code in J. C. Nash,
 * "Compact Numerical Methods for Computers," 2nd edition, converted
 * by p2c then re-crafted by B. D. Ripley for gnu R.  Revised for 
 * gretl by Allin Cottrell.
 * 
 * Returns: 0 on successful completion, non-zero error code
 * on error.
 */

int BFGS_orig (double *b, int n, int maxit, double reltol,
	       int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	       int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
	       gretlopt opt, PRN *prn)
{
    int crit_ok, done;
    double *g = NULL, *t = NULL, *X = NULL, *c = NULL, **H = NULL;
    int verbose = (opt & OPT_V);
    int ndelta, fcount, gcount;
    double d, fmax, f, f0, sumgrad;
    int i, j, ilast, iter;
    double s, steplen = 0.0;
    double D1, D2;
    int err = 0;

    BFGS_get_user_values(b, n, &maxit, &reltol, opt, prn);

    if (gradfunc == NULL) {
	gradfunc = BFGS_numeric_gradient;
    }

    g = malloc(n * sizeof *g);
    t = malloc(n * sizeof *t);
    X = malloc(n * sizeof *X);
    c = malloc(n * sizeof *c);
    H = triangular_array_new(n);

    if (g == NULL || t == NULL || X == NULL || c == NULL || H == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    f = cfunc(b, data);

    if (na(f)) {
	fprintf(stderr, "initial value of f is not finite\n");
	err = E_DATA;
	goto bailout;
    }

    f0 = fmax = f;
    iter = ilast = fcount = gcount = 1;
    gradfunc(b, g, n, cfunc, data);
    reverse_gradient(g, n);

    do {
	if (verbose) {
	    reverse_gradient(g, n);
	    print_iter_info(iter, f, crittype, n, b, g, steplen, prn);
	    reverse_gradient(g, n);
	}

	if (ilast == gcount) {
	    /* (re-)start: initialize curvature matrix */
	    for (i=0; i<n; i++) {
		for (j=0; j<i; j++) {
		    H[i][j] = 0.0;
		}
		H[i][i] = 1.0;
	    }
	}

	for (i=0; i<n; i++) {
	    /* copy coefficients to X, gradient to c */
	    X[i] = b[i];
	    c[i] = g[i];
	}

	sumgrad = 0.0;
	for (i=0; i<n; i++) {
	    s = 0.0;
	    for (j=0; j<=i; j++) {
		s -= H[i][j] * g[j];
	    }
	    for (j=i+1; j<n; j++) {
		s -= H[j][i] * g[j];
	    }
	    t[i] = s;
	    sumgrad += s * g[i];
	}

	if (sumgrad < 0.0) {	
	    steplen = 1.0;
	    crit_ok = 0;
	    do {
		/* loop so long as (a) we haven't achieved an
		   acceptable value of the criterion and (b) there is
		   still some prospect of doing so */
		ndelta = n;
		for (i=0; i<n; i++) {
		    b[i] = X[i] + steplen * t[i];
		    if (coeff_unchanged(b[i], X[i])) {
			ndelta--;
		    }
		}
		if (ndelta > 0) {
		    f = cfunc(b, data);
		    d = sumgrad * steplen * acctol;
		    fcount++;
		    crit_ok = !na(f) && (f >= fmax + d);
		    if (!crit_ok) {
			/* calculated criterion no good: try smaller step */
			steplen *= stepfrac;
		    }
		}
	    } while (ndelta != 0 && !crit_ok);

	    done = fabs(fmax - f) <= reltol * (fabs(fmax) + reltol);

#if BFGS_DEBUG
	    fprintf(stderr, "LHS=%g, RHS=%g; done = %d\n",
		    fabs(fmax - f), reltol * (fabs(fmax) + reltol),
		    done);
#endif

	    /* prepare to stop if relative change is small enough */
	    if (done) {
		ndelta = 0;
		fmax = f;
	    }

	    if (ndelta > 0) {
		/* making progress */
		fmax = f;
		gradfunc(b, g, n, cfunc, data);
		reverse_gradient(g, n);
		gcount++;
		iter++;
		D1 = 0.0;
		for (i=0; i<n; i++) {
		    t[i] = steplen * t[i];
		    c[i] = g[i] - c[i];
		    D1 += t[i] * c[i];
		}
		if (D1 > 0.0) {
		    D2 = 0.0;
		    for (i=0; i<n; i++) {
			s = 0.0;
			for (j=0; j<=i; j++) {
			    s += H[i][j] * c[j];
			}
			for (j=i+1; j<n; j++) {
			    s += H[j][i] * c[j];
			}
			X[i] = s;
			D2 += s * c[i];
		    }
		    D2 = 1.0 + D2 / D1;
		    for (i=0; i<n; i++) {
			for (j=0; j<=i; j++) {
			    H[i][j] += (D2 * t[i]*t[j] - X[i]*t[j] - t[i]*X[j]) / D1;
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
	} else {
	    /* heading in the wrong direction */
	    if (ilast == gcount) {
		/* we just reset: don't reset again; set ndelta = 0 so
		   that we exit the main loop
		*/
		ndelta = 0;
		if (gcount == 1) {
		    err = (broken_gradient(g, n))? E_NAN : E_DATA;
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

#if BFGS_DEBUG
    fprintf(stderr, "terminated: ndelta=%d, ilast=%d, gcount=%d\n",
	    ndelta, ilast, gcount);
#endif

    if (iter >= maxit) {
	fprintf(stderr, _("stopped after %d iterations\n"), iter);
	err = E_NOCONV;
    } else if (fmax < f0) {
	/* FIXME this should never happen */
	fprintf(stderr, "failed to match initial value of objective function, %g\n", f0);
	err = E_NOCONV;
    }

    *fncount = fcount;
    *grcount = gcount;

    if (verbose) {
	reverse_gradient(g, n);
	print_iter_info(-1, f, crittype, n, b, g, steplen, prn);
	pputc(prn, '\n');
    }

 bailout:

    free(g);
    free(t);
    free(X);
    free(c);
    free_triangular_array(H, n);

    return err;
}

int LBFGS_max (double *b, int n, int maxit, double reltol,
	       int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	       int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
	       gretlopt opt, PRN *prn)
{
    double *g = NULL;
    double *l = NULL;
    double *u = NULL;
    double *wa = NULL;
    int *nbd = NULL;
    int *iwa = NULL;

    int i, m, wadim;
    char task[60];
    char csave[60];
    double f, pgtol;
    double dsave[29];
    int isave[44];
    int lsave[4];
    int iter, ibak = 0;
    int err = 0;

    *fncount = *grcount = 0;    

    BFGS_get_user_values(b, n, &maxit, &reltol, opt, prn);

    /*
      m: the number of corrections used in the limited memory matrix.
      It is not altered by the routine.  Values of m < 3 are not
      recommended, and large values of m can result in excessive
      computing time. The range 3 <= m <= 20 is recommended.
    */
    m = 10; /* was initially set to 5 */

    wadim = (2*m+4)*n + 12*m*m + 12*m;

    g = malloc(n * sizeof *g);
    l = malloc(n * sizeof *l);
    u = malloc(n * sizeof *u);
    wa = malloc(wadim * sizeof *wa);
    nbd = malloc(n * sizeof *nbd);
    iwa = malloc(3*n * sizeof *iwa);

    if (g == NULL || l == NULL || u == NULL ||
	wa == NULL || nbd == NULL || iwa == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (gradfunc == NULL) {
	gradfunc = BFGS_numeric_gradient;
    }

    /* Gradient convergence criterion (not used -- we use reltol instead) */
    pgtol = 0;

    /* Bounds on the parameters: for now we just set them all to be
       less than some ridiculously large number */
    for (i=0; i<n; i++) {
	nbd[i] = 3; /* case 3: upper bound only */
	u[i] = NADBL / 100;
    }	

    /* Start the iteration by initializing 'task' */
    strcpy(task, "START");

    while (1) {
	/* Call the L-BFGS-B code */
	setulb_(&n, &m, b, l, u, nbd, &f, g, &reltol, &pgtol, wa, iwa, 
		task, csave, lsave, isave, dsave);

	iter = isave[29] + 1;

	if (!strncmp(task, "FG", 2)) {

	    /* Compute function value, f */
	    f = cfunc(b, data);
	    if (!na(f)) {
		f = -f;
	    } else if (*fncount == 0) {
		fprintf(stderr, "initial value of f is not finite\n");
		err = E_DATA;
		break;
	    }
	    *fncount += 1;

	    /* Compute gradient, g */
	    gradfunc(b, g, n, cfunc, data);
	    reverse_gradient(g, n);
	    *grcount += 1;
	    
	} else if (!strncmp(task, "NEW_X", 5)) {
	    /* The optimizer has produced a new set of parameter values */
	    if (isave[33] >= maxit) {
		strcpy(task, "STOP: TOTAL NO. of f AND g "
		       "EVALUATIONS EXCEEDS LIMIT");
		err = E_NOCONV;
		break;
	    } 
	} else {
	    fprintf(stderr, "%s\n", task);
	    break;
	}

	if (opt & OPT_V) {
	    if (iter != ibak) {
		double steplen = (iter == 1)? NADBL : dsave[13];

		reverse_gradient(g, n);
		print_iter_info(iter, -f, crittype, n, b, g, steplen, prn);
		reverse_gradient(g, n);
	    }
	    ibak = iter;
	}
    }

    if (!err && crittype == C_GMM) {
	/* finalize GMM computations */
	f = cfunc(b, data);
    }

    if (opt & OPT_V) {
	reverse_gradient(g, n);
	print_iter_info(-1, -f, crittype, n, b, g, dsave[13], prn);
	pputc(prn, '\n');
    }

 bailout:

    free(g);
    free(l);
    free(u);
    free(wa);
    free(nbd);
    free(iwa);

    return err;
}

/**
 * BFGS_max:
 * @b: array of adjustable coefficients.
 * @n: number elements in array @b.
 * @maxit: the maximum number of iterations to allow.
 * @reltol: relative tolerance for terminating iteration.
 * @fncount: location to receive count of function evaluations.
 * @grcount: location to receive count of gradient evaluations.
 * @cfunc: pointer to function used to calculate maximand.
 * @crittype: code for type of the maximand/minimand: should
 * be %C_LOGLIK, %C_GMM or %C_OTHER.  Used only in printing
 * iteration info.
 * @gradfunc: pointer to function used to calculate the 
 * gradient, or %NULL for default numerical calculation.
 * @data: pointer that will be passed as the last
 * parameter to the callback functions @cfunc and @gradfunc.
 * @opt: may contain %OPT_V for verbose operation.
 * @prn: printing struct (or %NULL).  Only used if @opt
 * includes %OPT_V.
 *
 * Obtains the set of values for @b which jointly maximize the
 * criterion value as calculated by @cfunc.  Uses the BFGS
 * variable-metric method.  Based on Pascal code in J. C. Nash,
 * "Compact Numerical Methods for Computers," 2nd edition, converted
 * by p2c then re-crafted by B. D. Ripley for gnu R.  Revised for 
 * gretl by Allin Cottrell.
 * 
 * Returns: 0 on successful completion, non-zero error code
 * on error.
 */

int BFGS_max (double *b, int n, int maxit, double reltol,
	      int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	      int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
	      gretlopt opt, PRN *prn)
{
    if (libset_get_bool(USE_LBFGS)) {
	return LBFGS_max(b, n, maxit, reltol,
			 fncount, grcount, cfunc, 
			 crittype, gradfunc, data, 
			 opt, prn);
    } else {
	return BFGS_orig(b, n, maxit, reltol,
			 fncount, grcount, cfunc, 
			 crittype, gradfunc, data, 
			 opt, prn);
    }
}

/* user-level access to BFGS */

typedef struct umax_ umax;

struct umax_ {
    gretl_matrix *b;
    int ncoeff;
    GENERATOR *g;
    double x_out;
    gretl_matrix *m_out;
    double ***Z;
    DATAINFO *dinfo;
    PRN *prn;
};

static void umax_init (umax *u)
{
    u->b = NULL;
    u->ncoeff = 0;
    u->g = NULL;
    u->x_out = NADBL;
    u->m_out = NULL;
    u->Z = NULL;
    u->dinfo = NULL;
    u->prn = NULL;
}

static void umax_clear (umax *u)
{
    /* drop any "$" variables created */
    dataset_drop_listed_variables(NULL, u->Z, u->dinfo, NULL);
    destroy_genr(u->g);
}

static double user_get_criterion (const double *b, void *p)
{
    umax *u = (umax *) p;
    double x = NADBL;
    int i, v;
    int err;

    for (i=0; i<u->ncoeff; i++) {
	u->b->val[i] = b[i];
    }

    err = execute_genr(u->g, u->Z, u->dinfo, u->prn); 

    if (err) {
	return NADBL;
    }

    v = genr_get_output_varnum(u->g);
    if (v > 0) {
	x = (*u->Z)[v][0];
    } else {
	gretl_matrix *m = genr_get_output_matrix(u->g);

	if (m != NULL && m->rows == 1 && m->cols == 1) {
	    x = m->val[0];
	}
    }

    u->x_out = x;
    
    return x;
}

static int user_gen_setup (umax *u,
			   const char *fncall,
			   double ***pZ, 
			   DATAINFO *pdinfo)
{
    char formula[MAXLINE];
    GENERATOR *g;
    int err = 0;

    sprintf(formula, "$umax=%s", fncall);

    g = genr_compile(formula, pZ, pdinfo, &err);

    if (!err) {
	/* see if the formula actually works */
	err = execute_genr(g, pZ, pdinfo, u->prn);
    }

    if (!err) {
	u->g = g;
	u->Z = pZ;
	u->dinfo = pdinfo;
	u->m_out = genr_get_output_matrix(g);
    } else {
	destroy_genr(g);
    }

    return err;
}

double user_BFGS (gretl_matrix *b, const char *fncall,
		  double ***pZ, DATAINFO *pdinfo,
		  PRN *prn, int *err)
{
    umax u;
    double ret = NADBL;
    gretlopt opt = OPT_NONE;
    int maxit;
    int fcount = 0, gcount = 0;
    double tol;

    umax_init(&u);

    u.ncoeff = gretl_vector_get_length(b);
    if (u.ncoeff == 0) {
	*err = E_DATA;
	goto bailout;
    }

    u.b = b;

    *err = user_gen_setup(&u, fncall, pZ, pdinfo);
    if (*err) {
	return NADBL;
    }

    maxit = libset_get_int(BFGS_MAXITER);
    tol = libset_get_double(BFGS_TOLER);

    if (libset_get_bool(MAX_VERBOSE)) {
	opt = OPT_V;
	u.prn = prn;
    }

    *err = BFGS_max(b->val, u.ncoeff, maxit, tol, &fcount, &gcount,
		    user_get_criterion, C_OTHER, NULL, &u, opt, prn);

    if (fcount > 0) {
	pprintf(prn, _("Function evaluations: %d\n"), fcount);
	pprintf(prn, _("Evaluations of gradient: %d\n"), gcount);
    }

    if (!*err) {
	ret = u.x_out;
    }

 bailout:

    umax_clear(&u);

    return ret;
}

static int user_calc_fvec (integer *m, integer *n, double *x, double *fvec,
			   integer *iflag, void *p)
{
    umax *u = (umax *) p;
    gretl_matrix *v;
    int i, err;

    for (i=0; i<*n; i++) {
	u->b->val[i] = x[i];
    }

    err = execute_genr(u->g, u->Z, u->dinfo, u->prn); 
    if (err) {
	fprintf(stderr, "execute_genr: err = %d\n", err); 
    }

    if (err) {
	*iflag = -1;
	return 0;
    }

    v = genr_get_output_matrix(u->g);

    if (v == NULL || gretl_vector_get_length(v) != *m) {
	fprintf(stderr, "user_calc_fvec: got bad matrix\n"); 
	*iflag = -1;
    } else {
	for (i=0; i<*m; i++) {
	    fvec[i] = v->val[i];
	}
    }
    
    return 0;
}

static int fdjac_allocate (integer m, integer n,
			   gretl_matrix **J, 
			   double **w, double **f)
{
    *J = gretl_matrix_alloc(m, n);
    if (*J == NULL) {
	return E_ALLOC;
    }

    *w = malloc(m * sizeof **w);
    *f = malloc(m * sizeof **f);

    if (*w == NULL || *f == NULL) {
	return E_ALLOC;
    }

    return 0;
}

gretl_matrix *fdjac (gretl_matrix *theta, const char *fncall,
		     double ***pZ, DATAINFO *pdinfo,
		     int *err)
{
    umax u;
    gretl_matrix *J = NULL;
    integer m, n;
    integer iflag = 0;
    double *wa = NULL;
    double *fvec = NULL;
    double epsfcn = 0.0;
    int i;

    *err = 0;

    umax_init(&u);

    n = gretl_vector_get_length(theta);
    if (n == 0) {
	*err = E_DATA;
	return NULL;
    }

    u.b = theta;
    u.ncoeff = n;

    *err = user_gen_setup(&u, fncall, pZ, pdinfo);
    if (*err) {
	fprintf(stderr, "ldjac: error %d from user_gen_setup\n", *err);
	goto bailout;
    }

    if (u.m_out == NULL) {
	*err = E_TYPES; /* FIXME */
	goto bailout;
    }

    m = gretl_vector_get_length(u.m_out);
    if (m == 0) {
	*err = E_DATA;
	goto bailout;
    }
    
    *err = fdjac_allocate(m, n, &J, &wa, &fvec);
    if (*err) {
	goto bailout;
    }

    for (i=0; i<m; i++) {
	fvec[i] = u.m_out->val[i];
    }

    fdjac2_(user_calc_fvec, &m, &n, theta->val, fvec, J->val, 
	    &m, &iflag, &epsfcn, wa, &u);

 bailout:

    free(wa);
    free(fvec);

    if (*err) {
	gretl_matrix_free(J);
	J = NULL;
    }

    umax_clear(&u);

    return J;
}
