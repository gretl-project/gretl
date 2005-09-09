/*
 *  Copyright (c) 2003-2005 by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* Nonlinear least squares for libgretl, using minpack */

#include "libgretl.h" 
#include "libset.h"

#include "f2c.h"
#include "../../minpack/minpack.h"  

#undef NLS_DEBUG

enum {
    NUMERIC_DERIVS,
    ANALYTIC_DERIVS
} nls_modes;

typedef struct _nls_param nls_param;

struct _nls_param {
    char name[VNAMELEN]; /* name of parameter (scalar variable) */
    char *deriv;         /* string representation of derivative of regression
			    function with respect to param (or NULL) */
    int varnum;          /* ID number of the scalar variable in dataset */
    int dernum;          /* ID number of the variable holding the derivative */
};

struct _nls_spec {
    int ci;             /* NLS or MLE */
    int mode;           /* derivatives: numeric or analytic */
    int depvar;         /* ID number of dependent variable */
    int uhatnum;        /* ID number of variable holding residuals */
    char *nlfunc;       /* string representation of nonlinear function,
			   expressed in terms of the residuals */
    int nparam;         /* number of parameters to be estimated */
    int iters;          /* number of iterations performed */
    int t1;             /* starting observation */
    int t2;             /* ending observation */
    double ess;         /* error sum of squares */
    double ll;          /* log likelihood */
    double tol;         /* tolerance for stopping iteration */
    nls_param *params;  /* array of information on function parameters
			   (see the _nls_param struct above) */
    doublereal *coeff;  /* coefficient estimates */
};

/* file-scope global variables: we need to access these variables in
   the context of the callback functions that we register with the
   minpack routines, and there is no provision in minpack for passing
   additional pointers into the callbacks.
*/

static double ***nZ;
static DATAINFO *ndinfo;
static PRN *nprn;

static nls_spec private_spec;
static nls_spec *pspec;
static integer one = 1;
static int genr_err;

/* Use the genr() function to create/update the values of the residual
   from the regression function (if i = 0), or the derivatives of the
   regression function with respect to the parameters.  The formulae
   for generating these values are stored in string form in the nlspec
   struct.
*/

static int nls_auto_genr (int i, double ***pZ, DATAINFO *pdinfo)
{
    char formula[MAXLINE];

    if (i == 0) {
	/* note: $nl_y is the residual */
	sprintf(formula, "$nl_y = %s", pspec->nlfunc);
    } else {
	/* derivatives: artificial independent variables */
	sprintf(formula, "$nl_x%d = %s", i, pspec->params[i-1].deriv);
    }

    if (pZ != NULL && pdinfo != NULL) {
	genr_err = generate(formula, pZ, pdinfo, NULL, OPT_P);
    } else {
	/* using "global" nZ and ndinfo pointers here */
	genr_err = generate(formula, nZ, ndinfo, NULL, OPT_P);
    }

#if NLS_DEBUG
    if (genr_err) {
	errmsg(genr_err, nprn);
    } else {
	fprintf(stderr, "nls_auto_genr: i=%d, formula='%s'\n", i, formula);
    }
#endif

    return genr_err;
}

/* three wrappers for the above to enhance comprehensibility below */

static int nls_calculate_fvec (void)
{
    return nls_auto_genr(0, NULL, NULL);
}

static int nls_test_calculate_fvec (double ***pZ, DATAINFO *pdinfo)
{
    return nls_auto_genr(0, pZ, pdinfo);
}

static int nls_calculate_deriv (int i)
{
    return nls_auto_genr(i + 1, NULL, NULL);
}

static int nlspec_allocate_param (nls_spec *spec)
{
    nls_param *params;
    double *coeff;
    int nt = spec->nparam + 1;

    params = realloc(spec->params, nt * sizeof *spec->params);
    if (params == NULL) {
	return 1;
    }

    spec->params = params;

    coeff = realloc(spec->coeff, nt * sizeof *spec->coeff);
    if (coeff == NULL) {
	free(params);
	return 1;
    }

    spec->coeff = coeff;

    spec->nparam += 1;

    return 0;
}

/* allocate space for an additional regression parameter in the
   nls_spec struct and add its info */

static int 
real_add_param_to_spec (const char *vname, int vnum, double initval,
			nls_spec *spec)
{
    int i;

#if NLS_DEBUG
    fprintf(stderr, "real_add_param: adding '%s'\n", vname);
#endif

    if (nlspec_allocate_param(spec)) {
	return E_ALLOC;
    }

    i = spec->nparam - 1;
    
    spec->params[i].varnum = vnum;
    spec->params[i].dernum = 0;
    strcpy(spec->params[i].name, vname);
    spec->params[i].deriv = NULL;

    spec->coeff[i] = initval;

    return 0;
}

/* scrutinize word and see if it's a new scalar that should
   be added to the NLS specification; if so, add it */

static int 
maybe_add_param_to_spec (nls_spec *spec, const char *word, 
			 const double **Z, const DATAINFO *pdinfo)
{
    int i, v;

#if NLS_DEBUG
    fprintf(stderr, "maybe_add_param: looking at '%s'\n", word);
#endif

    /* if word represents a math function or constant, skip it */
    if (genr_function_from_string(word) || !strcmp(word, "pi")) {
	return 0;
    }

    /* try looking up word as the name of a variable */
    v = varindex(pdinfo, word);
    if (v >= pdinfo->v) {
	sprintf(gretl_errmsg, _("Unknown variable '%s'"), word);
	return E_UNKVAR;
    }

    /* if term is not a scalar, skip it: only scalars can figure
       as regression parameters */
    if (pdinfo->vector[v]) {
	return 0;
    }

    /* if this term is already present in the specification, skip it */
    for (i=0; i<spec->nparam; i++) {
	if (strcmp(word, spec->params[i].name) == 0) {
	    return 0;
	}
    }

    /* else: add this param to the NLS specification */
    return real_add_param_to_spec(word, v, Z[v][0], spec);
}

/* Parse NLS function specification string to find names of variables
   that may figure as parameters of the regression function.  We need
   to do this only if analytical derivatives have not been supplied.
*/

static int 
get_params_from_nlfunc (nls_spec *spec, const double **Z,
			const DATAINFO *pdinfo)
{
    const char *s = spec->nlfunc;
    const char *p;
    char vname[VNAMELEN];
    int n, np = 0;
    int err = 0;

#if NLS_DEBUG
    fprintf(stderr, "get_params: looking at '%s'\n", s);
#endif

    while (*s && !err) {
	p = s;
	if (isalpha(*s) && *(s + 1)) { 
	    /* find a variable name */
	    p++;
	    n = 1;
	    while (isalnum(*p) || *p == '_') {
		p++; n++;
	    }
	    if (n > 8) {
		/* variable name is too long */
		return 1;
	    }
	    *vname = 0;
	    strncat(vname, s, n);
	    if (np > 0) {
		/* right-hand side term */
		err = maybe_add_param_to_spec(spec, vname, Z, pdinfo);
	    }
	    np++;
	} else {
	    p++;
	}
	s = p;
    }

    return err;
}

/* Adjust starting and ending points of sample if need be, to avoid
   missing values; abort if there are missing values within the
   (possibly reduced) sample range.  For this purpose we generate
   the nls residual variable.
*/

static int 
nls_missval_check (double ***pZ, DATAINFO *pdinfo, nls_spec *spec)
{
    int t, v, miss = 0;
    int t1 = spec->t1, t2 = spec->t2;
    int err = 0;

    /* generate the nls residual variable */
    err = nls_test_calculate_fvec(pZ, pdinfo);
    if (err) {
	return err;
    }
	
    /* ID number of last variable generated, above */
    v = pdinfo->v - 1;

    for (t=spec->t1; t<=spec->t2; t++) {
	if (na((*pZ)[v][t])) {
	    t1++;
	} else {
	    break;
	}
    }

    for (t=spec->t2; t>=spec->t1; t--) {
	if (na((*pZ)[v][t])) {
	    t2--;
	} else {
	    break;
	}
    }

    if (t2 - t1 + 1 < spec->nparam) {
	return E_DF;
    }

    for (t=t1; t<=t2; t++) {
	if (na((*pZ)[v][t])) {
	    miss = 1;
	    break;
	}
    }  

    if (miss) {
	strcpy(gretl_errmsg, _("There were missing data values"));
	return 1;
    }

    spec->t1 = t1;
    spec->t2 = t2;

    return 0;
}

/* this function is used in the context of the minpack callback */

static int get_nls_fvec (double *fvec)
{
    int j, t, v;

    /* calculate residual given current parameter estimates */
    if (nls_calculate_fvec()) {
	return 1;
    }

    if (pspec->uhatnum == 0) {
	/* look up ID number of the fvec variable if we don't know
	   it already */
	v = varindex(ndinfo, "$nl_y");
	if (v == ndinfo->v) {
	    return 1;
	}
	pspec->uhatnum = v;
    } else {
	v = pspec->uhatnum;
    }

    pspec->ess = 0.0;
    pspec->ll = 0.0;
    j = 0;

    /* transcribe from dataset to fvec array */
    for (t=pspec->t1; t<=pspec->t2; t++) {
	fvec[j] = (*nZ)[v][t];
	if (pspec->ci == MLE) {
	    pspec->ll -= fvec[j];
	} else {
	    pspec->ess += fvec[j] * fvec[j];
	}
	j++;
    }

    pspec->iters += 1;

    if (pspec->ci == MLE) {
	pprintf(nprn, "iteration %2d: ll = %.8g\n", pspec->iters, pspec->ll);
    }

#if NLS_DEBUG
    if (pspec->ci == NLS) {
	fprintf(stderr, "iteration %2d: SSR = %.8g\n", pspec->iters, pspec->ess);
    }
#endif

    return 0;
}

/* this function is used in the context of the minpack callback */

static int get_nls_deriv (int i, double *deriv)
{
    int j, t, v, vec;

    /* calculate value of deriv with respect to param */
    if (nls_calculate_deriv(i)) {
	return 1;
    }

    if (pspec->params[i].dernum == 0) {
	char varname[VNAMELEN];

	/* look up ID number of the derivative if not known */
	sprintf(varname, "$nl_x%d", i + 1);
	v = varindex(ndinfo, varname);
	if (v == ndinfo->v) {
	    return 1;
	}
	pspec->params[i].dernum = v;
    } else {
	v = pspec->params[i].dernum;
    }

    /* derivative may be vector or scalar */
    vec = ndinfo->vector[v];

    j = 0;
    /* transcribe from dataset to deriv array */
    for (t=pspec->t1; t<=pspec->t2; t++) {
	if (vec) {
	    deriv[j] = - (*nZ)[v][t];
	} else {
	    deriv[j] = - (*nZ)[v][0];
	}
	j++;
    }

    return 0;
}

/* compute auxiliary statistics and add them to the NLS 
   model struct */

static void add_stats_to_model (MODEL *pmod, nls_spec *spec,
				const double **Z)
{
    int dv = spec->depvar;
    double d, tss;
    int t;

    pmod->ess = spec->ess;
    pmod->sigma = sqrt(spec->ess / (pmod->nobs - spec->nparam));
    
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

/* add coefficient covariance matrix and standard errors */

static int add_std_errs_to_model (MODEL *pmod)
{
    int i, k;

    if (pmod->vcv == NULL && makevcv(pmod)) {
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
add_param_names_to_model (MODEL *pmod, nls_spec *spec, const DATAINFO *pdinfo)
{
    int i;
    int np = pmod->ncoeff + 1;

    pmod->params = malloc(np * sizeof *pmod->params);
    if (pmod->params == NULL) {
	return 1;
    }

    pmod->nparams = np;

    pmod->params[0] = malloc(VNAMELEN);
    if (pmod->params[0] == NULL) {
	free(pmod->params);
	return 1;
    }

    strcpy(pmod->params[0], pdinfo->varname[spec->depvar]);

    for (i=1; i<=pmod->ncoeff; i++) {
	pmod->params[i] = malloc(VNAMELEN);
	if (pmod->params[i] == NULL) {
	    int j;

	    for (j=0; j<i; j++) {
		free(pmod->params[j]);
	    }
	    free(pmod->params);
	    pmod->params = NULL;
	    pmod->nparams = 0;
	    return 1;
	}
	strcpy(pmod->params[i], spec->params[i-1].name);
    }

    return 0;
}

static void 
add_fit_resid_to_model (MODEL *pmod, nls_spec *spec, double *uhat, 
			const double **Z, int perfect)
{
    int t, j = 0;

    if (perfect) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    pmod->uhat[t] = 0.0;
	    pmod->yhat[t] = Z[spec->depvar][t];
	}
    } else {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    pmod->uhat[t] = uhat[j];
	    pmod->yhat[t] = Z[spec->depvar][t] - uhat[j];
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
	gretl_model_set_data(pmod, "nl_regfunc", formula, 
			     strlen(formula) + 1);
    } else {
	err = E_ALLOC;
    }

    return err;
}

/* Gauss-Newton regression to calculate standard errors for the
   NLS parameters (see Davidson and MacKinnon) */

static MODEL GNR (double *uhat, double *jac, nls_spec *spec,
		  const double **Z, const DATAINFO *pdinfo, 
		  PRN *prn)
{
    double **gZ = NULL;
    DATAINFO *gdinfo;
    int *glist;
    MODEL gnr;
    int i, j, t;
    int iters = spec->iters;
    int perfect = 0;
    int err = 0;

    if (gretl_iszero(spec->t1, spec->t2, uhat)) {
	pputs(prn, _("Perfect fit achieved\n"));
	perfect = 1;
	for (t=spec->t1; t<=spec->t2; t++) {
	    uhat[t] = 1.0;
	}
	spec->ess = 0.0;
    }

    gdinfo = create_new_dataset(&gZ, spec->nparam + 1, pdinfo->n, 0);
    if (gdinfo == NULL) {
	gretl_model_init(&gnr);
	gnr.errcode = E_ALLOC;
	return gnr;
    }

    /* transcribe sample info */
    gdinfo->t1 = spec->t1;
    gdinfo->t2 = spec->t2;
    
    glist = gretl_list_new(spec->nparam + 1);

    if (glist == NULL) {
	destroy_dataset(gZ, gdinfo);
	gretl_model_init(&gnr);
	gnr.errcode = E_ALLOC;
	return gnr;
    }
    
    for (i=0; i<=spec->nparam; i++) {
	glist[i+1] = i;
	if (i == 0) {
	    /* dependent variable (NLS residual) */
	    j = 0;
	    for (t=0; t<gdinfo->n; t++) {
		if (t < gdinfo->t1 || t > gdinfo->t2) {
		    gZ[i][t] = NADBL;
		} else {
		    gZ[i][t] = uhat[j++];
		}
	    }
	} else {
	    /* independent vars: derivatives wrt params */
	    if (spec->mode == ANALYTIC_DERIVS) {
		get_nls_deriv(i-1, gZ[i]);
	    } else {
		j = gdinfo->n * (i - 1);
		for (t=0; t<gdinfo->n; t++) {
		    if (t < gdinfo->t1 || t > gdinfo->t2) {
			gZ[i][t] = NADBL;
		    } else {
			gZ[i][t] = jac[j++];
		    }
		}
	    }
	}
    }

    gnr = lsq(glist, &gZ, gdinfo, OLS, OPT_A, 0.0);

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
	add_stats_to_model(&gnr, spec, Z);
	if (add_std_errs_to_model(&gnr)) {
	    gnr.errcode = E_ALLOC;
	}
    }

    if (gnr.errcode == 0) {
	ls_aic_bic(&gnr);
	gnr.ci = NLS;
	add_coeffs_to_model(&gnr, spec->coeff);
	add_param_names_to_model(&gnr, spec, pdinfo);
	add_fit_resid_to_model(&gnr, spec, uhat, Z, perfect);
	gnr.list[1] = spec->depvar;

	/* set relevant data on model to be shipped out */
	gretl_model_set_int(&gnr, "iters", iters);
	gretl_model_set_double(&gnr, "tol", spec->tol);
	transcribe_nls_function(&gnr, pspec->nlfunc);
    }

    destroy_dataset(gZ, gdinfo);
    free(glist);

    return gnr;
}

/* free up resources associated with the nlspec struct */

static void clear_nls_spec (nls_spec *spec)
{
    int i;

    if (spec == NULL) {
	return;
    }

    if (spec->params != NULL) {
	for (i=0; i<spec->nparam; i++) {
	    free(spec->params[i].deriv);
	}
	free(spec->params);
	spec->params = NULL;
    }

    free(spec->nlfunc);
    spec->nlfunc = NULL;

    free(spec->coeff);
    spec->coeff = NULL;

    spec->ci = NLS;
    spec->mode = NUMERIC_DERIVS;
    spec->nparam = 0;

    spec->depvar = 0;
    spec->uhatnum = 0;

    spec->iters = 0;
    spec->t1 = spec->t2 = 0;
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

static void update_nls_param_values (const double *x)
{
    int i, v;

    /* write the values produced by minpack into the dataset */
    for (i=0; i<pspec->nparam; i++) {
	v = pspec->params[i].varnum;
	(*nZ)[v][0] = x[i];
	fprintf(stderr, "revised param[%d] = %g\n", i, x[i]);
    }
}

/* callback for lm_calculate (below) to be used by minpack */

static int nls_calc (integer *m, integer *n, double *x, double *fvec, 
		     double *jac, integer *ldjac, integer *iflag)
{
    int T = *m;
    int i;

#if NLS_DEBUG
    fprintf(stderr, "nls_calc called by minpack with iflag = %d\n", 
	    (int) *iflag);
#endif

    /* write current parameter values into dataset Z */
    update_nls_param_values(x);

    if (*iflag == 1) {
	/* calculate function at x, results into fvec */
	if (get_nls_fvec(fvec)) {
	    *iflag = -1;
	}
    } else if (*iflag == 2) {
	/* calculate jacobian at x, results into jac */
	for (i=0; i<*n; i++) {
	    if (get_nls_deriv(i, &jac[i*T])) {
		*iflag = -1; 
	    }
	}	
    }

    return 0;
}

/* in case the user supplied analytical derivatives for the
   NLS regression parameters, check them for sanity */

static int check_nls_derivs (integer m, integer n, double *x,
			     double *fvec, double *jac,
			     integer ldjac, PRN *prn)
{
    integer mode = 1;
    integer iflag;
    doublereal *xp = NULL;
    doublereal *err = NULL;
    doublereal *fvecp = NULL;
    int i, badcount = 0, zerocount = 0;

    xp = malloc(m * sizeof *xp);
    err = malloc(m * sizeof *err);
    fvecp = malloc(m * sizeof *fvecp);

    if (xp == NULL || err == NULL || fvecp == NULL) {
	free(err);
	free(xp);
	free(fvecp);
	return 1;
    }

    iflag = 1;
    nls_calc(&m, &n, x, fvec, jac, &ldjac, &iflag);
    if (iflag == -1) goto chkderiv_abort;

    chkder_(&m, &n, x, fvec, jac, &ldjac, xp, fvecp, &mode, err);

    iflag = 2;
    nls_calc(&m, &n, x, fvec, jac, &ldjac, &iflag);
    if (iflag == -1) goto chkderiv_abort;

    iflag = 1;
    nls_calc(&m, &n, xp, fvecp, jac, &ldjac, &iflag);
    if (iflag == -1) goto chkderiv_abort; 

    mode = 2;
    chkder_(&m, &n, x, fvec, jac, &ldjac, xp, fvecp, &mode, err);

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
	fprintf(stderr, "%d out of %d tests gave zero\n", zerocount, (int) m);
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

    return (zerocount > m/4);
}

/* driver for minpack levenberg-marquandt code for use when analytical
   derivatives have been supplied */

static int lm_calculate (nls_spec *spec, double *fvec, double *jac, PRN *prn)
{
    integer info, lwa;
    integer m, n, ldjac;
    integer *ipvt;
    doublereal *wa;
    int err = 0;

    m = spec->t2 - spec->t1 + 1; /* number of observations */
    n = spec->nparam;            /* number of parameters */
    lwa = 5 * n + m;             /* work array size */
    ldjac = m;                   /* leading dimension of jac array */

    wa = malloc(lwa * sizeof *wa);
    ipvt = malloc(n * sizeof *ipvt);

    if (wa == NULL || ipvt == NULL) {
	err = E_ALLOC;
	goto nls_cleanup;
    }

    err = check_nls_derivs(m, n, spec->coeff, fvec, jac, ldjac, prn);
    if (err) {
	goto nls_cleanup; 
    }

    /* call minpack */
    if (spec->ci == MLE) {
	lmmle1_(nls_calc, &m, &n, spec->coeff, fvec, jac, &ldjac, &spec->tol, 
		&info, ipvt, wa, &lwa);
    } else {
	lmder1_(nls_calc, &m, &n, spec->coeff, fvec, jac, &ldjac, &spec->tol, 
		&info, ipvt, wa, &lwa);
    }

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
		 integer *iflag)
{
    /* write current parameter values into dataset Z */
    update_nls_param_values(x);

    /* calculate function at x, results into fvec */    
    if (get_nls_fvec(fvec)) {
	*iflag = -1;
    }

    return 0;
}

/* driver for minpack levenberg-marquandt code for use when the
   Jacobian must be approximated numerically */

static int 
lm_approximate (nls_spec *spec, double *fvec, double *jac, PRN *prn)
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
    
    m = spec->t2 - spec->t1 + 1; /* number of observations */
    n = spec->nparam;            /* number of parameters */
    ldjac = m;                   /* leading dimension of jac array */

    maxfev = 200 * (n + 1);

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
	   &spec->tol, &spec->tol, &gtol, &maxfev, &epsfcn, diag, &mode, &factor,
	   &nprint, &info, &nfev, jac, &ldjac, 
	   ipvt, qtf, wa1, wa2, wa3, wa4);

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
	double ess = spec->ess;
	int iters = spec->iters;

	/* call minpack again */
	fdjac2_(nls_calc_approx, &m, &n, spec->coeff, fvec, jac, 
		&ldjac, &iflag, &epsfcn, wa4);
	spec->ess = ess;
	spec->iters = iters;
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

/* below: public functions */

/**
 * nls_spec_add_param_with_deriv:
 * @spec: pointer to nls specification.
 * @dstr: string specifying a derivative with respect to a
 *   parameter of the regression function.
 * @Z: data array.
 * @pdinfo: information on dataset.
 *
 * Adds an analytical derivative to @spec.  This pointer must
 * have previously been obtained by a call to #nls_spec_new.
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
nls_spec_add_param_with_deriv (nls_spec *spec, const char *dstr,
			       const double **Z, const DATAINFO *pdinfo)
{
    nls_param *param = NULL;
    const char *p = dstr;
    char *vname = NULL;
    int i, v, err = 0;

    if (nlspec_allocate_param(spec)) {
	return E_ALLOC;
    }

    i = spec->nparam - 1;

    param = &spec->params[i];

    if (!strncmp(p, "deriv ", 6)) {
	/* make starting with "deriv" optional */
	p += 6;
    }

    err = equation_get_lhs_and_rhs(p, &vname, &param->deriv);
    if (err) {
	fprintf(stderr, "parse error in deriv string: '%s'\n", dstr);
	return E_PARSE;
    }

    *param->name = '\0';
    strncat(param->name, vname, 8);
    free(vname);

    v = varindex(pdinfo, param->name);
    if (v < pdinfo->v) {
	param->varnum = v;
	param->dernum = 0;
	spec->coeff[i] = Z[v][0];
    } else {
	free(param->deriv);
	param->deriv = NULL;
	sprintf(gretl_errmsg, _("Unknown variable '%s'"), param->name);
	err = E_UNKVAR;
    }

    if (!err) {
	spec->mode = ANALYTIC_DERIVS;
    }

#if NLS_DEBUG
    if (param->deriv != NULL) {
	fprintf(stderr, "add_param_with_deriv: '%s'\n"
		" set varnum = %d, initial value = %g\n", dstr, 
		param->varnum, spec->coeff[i]);
    }
#endif

    return err;
}

/**
 * nls_spec_set_regression_function:
 * @spec: pointer to nls specification.
 * @fnstr: string specifying nonlinear regression function.
 * @pdinfo: information on dataset.
 *
 * Adds the regression function to @spec.  This pointer must
 * have previously been obtained by a call to #nls_spec_new.
 * The required format for @fnstr is "%varname = %formula", where
 * %varname is the name of the dependent variable and %formula
 * is an expression of the sort that is fed to gretl's %genr command.
 * The dependent variable must be already present in the
 * dataset.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int 
nls_spec_set_regression_function (nls_spec *spec, const char *fnstr, 
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

    if (!strncmp(p, "nls ", 4) || !strncmp(p, "mle ", 4)) {
	/* starting with "nls" / "mle" is optional here */
	p += 4;
    }    

    if (equation_get_lhs_and_rhs(p, &vname, &rhs)) { 
	sprintf(gretl_errmsg, _("parse error in '%s'\n"), fnstr);
	err = E_PARSE;
    }

    if (!err) {
	/* FIXME mle */
	spec->depvar = varindex(pdinfo, vname);

	if (spec->depvar == pdinfo->v) {
	    sprintf(gretl_errmsg, _("Unknown variable '%s'"), vname);
	    err = E_UNKVAR;
	}
    }

    if (spec->ci == MLE) {
	flen = strlen(rhs) + 4;
    } else {
	flen = strlen(vname) + strlen(rhs) + 6;
    }

    spec->nlfunc = malloc(flen);
    if (spec->nlfunc == NULL) {
	err = E_ALLOC;
    } else {
	if (spec->ci == MLE) {
	    sprintf(spec->nlfunc, "-(%s)", rhs);
	} else {
	    sprintf(spec->nlfunc, "%s - (%s)", vname, rhs);
	}
    }

    free(vname);
    free(rhs);

    return err;
}

/**
 * nls_spec_set_t1_t2:
 * @spec: pointer to nls specification.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Sets the sample range for estimation of @spec.  This pointer must
 * have previously been obtained by a call to #nls_spec_new.
 */

void nls_spec_set_t1_t2 (nls_spec *spec, int t1, int t2)
{
    if (spec != NULL) {
	spec->t1 = t1;
	spec->t2 = t2;
    }
}

/**
 * nls_parse_line:
 * @ci: either %NLS or %MLE (doco not finished on this)
 * @line: specification of regression function or derivative
 *        of this function with respect to a parameter.
 * @Z: data array.
 * @pdinfo: information on dataset.
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
		    const DATAINFO *pdinfo)
{
    int err = 0;

    pspec = &private_spec;
    pspec->ci = ci;

    if (strncmp(line, "deriv", 5) == 0) {
	if (pspec->nlfunc == NULL) {
	    strcpy(gretl_errmsg, _("No regression function has been specified"));
	    err = E_PARSE;
	} else {
	    err = nls_spec_add_param_with_deriv(pspec, line, Z, pdinfo);
	}
    } else {
	/* do we already have an nls specification under way? */
	if (pspec->nlfunc != NULL) {
	    clear_nls_spec(pspec);
	}
	err = nls_spec_set_regression_function(pspec, line, pdinfo);
	if (!err) {
	    nls_spec_set_t1_t2(pspec, pdinfo->t1, pdinfo->t2);
	}
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

/* static function providing the real content for the two public
   wrapper functions below */

static MODEL real_nls (nls_spec *spec, double ***pZ, DATAINFO *pdinfo, 
		       PRN *prn)
{
    MODEL nlsmod;
    double *fvec = NULL;
    double *jac = NULL;
    int origv = pdinfo->v;
    int err = 0;

    genr_err = 0;
    gretl_model_init(&nlsmod);
    gretl_model_smpl_init(&nlsmod, pdinfo);

    if (pspec->nlfunc == NULL) {
	strcpy(gretl_errmsg, _("No regression function has been specified"));
	nlsmod.errcode = E_PARSE;
	goto bailout;
    } 

    /* publish pdinfo and prn */
    ndinfo = pdinfo;
    nprn = prn;

    if (spec != NULL) {
	/* the caller supplied an nls specification directly */
	pspec = spec;
    } else {
	/* we use the static spec composed via nls_parse_line() */
	pspec = &private_spec;
    }

    if (pspec->mode == NUMERIC_DERIVS) {
	err = get_params_from_nlfunc(pspec, (const double **) *pZ, pdinfo);
	if (err) {
	    if (err == 1) {
		nlsmod.errcode = E_PARSE;
	    } else {
		nlsmod.errcode = err;
	    }
	    goto bailout;
	}
    }

    if (pspec->nparam == 0) {
	strcpy(gretl_errmsg, _("No regression function has been specified"));
	nlsmod.errcode = E_PARSE;
	goto bailout;
    } 

    err = nls_missval_check(pZ, pdinfo, pspec);
    if (err) {
	nlsmod.errcode = err;
	goto bailout;
    }

    /* allocate (full-length) arrays to be passed to minpack */
    fvec = malloc(pdinfo->n * sizeof *fvec);
    jac = malloc(pdinfo->n * pspec->nparam * sizeof *jac);

    if (fvec == NULL || jac == NULL) {
	nlsmod.errcode = E_ALLOC;
	goto bailout;
    }

    /* get tolerance from user setting or default */
    pspec->tol = get_nls_toler();

    /* export Z pointer for minpack's benefit */
    nZ = pZ;

    /* invoke minpack driver function */
    if (pspec->mode == NUMERIC_DERIVS) {
	pputs(prn, _("Using numerical derivatives\n"));
	err = lm_approximate(pspec, fvec, jac, prn);
    } else {
	pputs(prn, _("Using analytical derivatives\n"));
	err = lm_calculate(pspec, fvec, jac, prn);
    }

    /* re-attach data array pointer: may have moved! */
    *pZ = *nZ;

    pprintf(prn, _("Tolerance = %g\n"), pspec->tol);

    if (!err) {
	/* We have parameter estimates; now use a Gauss-Newton
	   regression for covariance matrix, standard errors. */
	nlsmod = GNR(fvec, jac, pspec, (const double **) *pZ, 
		     pdinfo, prn);
    } else {
	if (nlsmod.errcode == 0) { 
	    if (genr_err != 0) {
		nlsmod.errcode = genr_err;
	    } else {
		nlsmod.errcode = E_NOCONV;
	    }
	}
    }

 bailout:

    free(fvec);
    free(jac);

    if (spec == NULL) {
	clear_nls_spec(pspec);
    }

    dataset_drop_last_variables(pdinfo->v - origv, pZ, pdinfo);

    if (nlsmod.errcode == 0) {
	set_model_id(&nlsmod);
    }

    return nlsmod;
}

/**
 * nls:
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 * @prn: printing struct.
 *
 * Computes estimates of a model via nonlinear least squares.
 * The model must have been specified previously, via calls to
 * the function #nls_parse_line.  
 *
 * Returns: a model struct containing the parameter estimates
 * and associated statistics.
 */

MODEL nls (double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    return real_nls(NULL, pZ, pdinfo, prn);
}

/**
 * model_from_nls_spec:
 * @spec: nls specification.
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 * @prn: printing struct.
 *
 * Computes estimates of the model specified in @spec, via nonlinear 
 * least squares. The @spec must first be obtained using #nls_spec_new, and
 * initialized using #nls_spec_set_regression_function.  If analytical
 * derivatives are to be used (which is optional but recommended)
 * these are set using #nls_spec_add_param_with_deriv.
 *
 * Returns: a model struct containing the parameter estimates
 * and associated statistics.
 */

MODEL model_from_nls_spec (nls_spec *spec, double ***pZ, DATAINFO *pdinfo, 
			   PRN *prn)
{
    return real_nls(spec, pZ, pdinfo, prn);
}

/**
 * nls_spec_new:
 * @ci: either %NLS or %MLE.
 * @pdinfo: information on dataset.
 *
 * Returns: a pointer to a newly allocated nls specification,
 * or %NULL on failure.
 */

nls_spec *nls_spec_new (int ci, const DATAINFO *pdinfo)
{
    nls_spec *spec;

    spec = malloc(sizeof *spec);
    if (spec == NULL) {
	return NULL;
    }

    spec->nlfunc = NULL;

    spec->params = NULL;
    spec->nparam = 0;
    
    spec->coeff = NULL;

    spec->ci = ci;
    spec->mode = NUMERIC_DERIVS;

    spec->nparam = 0;
    spec->depvar = 0;
    spec->iters = 0;

    spec->t1 = pdinfo->t1;
    spec->t2 = pdinfo->t2;

    return spec;
}

/**
 * nls_spec_destroy:
 * @spec: pointer to nls specification.
 *
 * Frees all resources associated with @spec, and frees the
 * pointer itself.
 */

void nls_spec_destroy (nls_spec *spec)
{
    clear_nls_spec(spec);
    free(spec);
}







