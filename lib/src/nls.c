/*
 *  Copyright (c) 2003 by Allin Cottrell
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "libgretl.h" 
#include "gretl_private.h"
#include "f2c.h"
#include "../../minpack/minpack.h"  

enum {
    NUMERIC_DERIVS,
    ANALYTIC_DERIVS
} nls_modes;

typedef struct _nls_spec nls_spec;
typedef struct _nls_term nls_term;

struct _nls_term {
    char name[VNAMELEN]; /* name of parameter */
    char *deriv;         /* string representation of derivative */
    int varnum;          /* ID number of the var holding the derivative */
};

struct _nls_spec {
    int mode;           /* derivatives: numeric or analytic */
    int depvar;         /* ID number of dependent variable */
    char *nlfunc;       /* string representation of regression function */
    int nparam;         /* number of parameters to be estimated */
    int iters;          /* number of iterations performed */
    int t1;             /* starting observation */
    int t2;             /* ending observation */
    double ess;         /* error sum of squares */
    double tol;         /* tolerance for stopping iteration */
    nls_term *terms;    /* array of info on terms in the function */
    doublereal *coeff;  /* coefficient estimates */
};

/* file-scope global variables */
static double ***pZ;
static DATAINFO *pdinfo;
static PRN *prn;
static nls_spec nlspec;
integer one = 1;
double toler;

static void print_iter_ess (void)
{
    nlspec.iters += 1;    
#if 0
    pprintf(prn, "iteration %2d: SSR = %g\n", nlspec.iters, nlspec.ess);
#endif
}

static int nls_auto_gen (int i)
{
    char formula[MAXLEN];
    int err;

    if (i == 0) {
	sprintf(formula, "$nls_y = %s", nlspec.nlfunc);
    } else {
	sprintf(formula, "$nls_x%d = %s", i, nlspec.terms[i-1].deriv);
    }

    err = generate(pZ, pdinfo, formula, 0, NULL, 0);

    if (err) {
	errmsg(err, prn);
    }

    return err;
}

static int add_term_from_nlfunc (const char *vname)
{
    nls_term *terms;
    double *coeff;
    int i, v, nt = nlspec.nparam + 1; 

    /* if term is math function or constant, skip */
    if (get_function(vname) || !strcmp(vname, "pi")) return 0;

    v = varindex(pdinfo, vname);
    if (v >= pdinfo->v) {
	sprintf(gretl_errmsg, _("Unknown variable '%s'"), vname);
	return E_UNKVAR;
    }

    /* if term is not a scalar, skip it */
    if (pdinfo->vector[v]) return 0;

    /* if term is already present, skip it */
    for (i=0; i<nlspec.nparam; i++) {
	if (strcmp(vname, nlspec.terms[i].name) == 0) 
	    return 0;
    }

    terms = realloc(nlspec.terms, nt * sizeof *nlspec.terms);
    if (terms == NULL) return E_ALLOC;

    coeff = realloc(nlspec.coeff, nt * sizeof *nlspec.coeff);
    if (coeff == NULL) {
	free(terms);
	return E_ALLOC;
    }

    nlspec.terms = terms;
    nlspec.coeff = coeff;

    terms[nt-1].varnum = v;
    strcpy(terms[nt-1].name, vname);
    terms[nt-1].deriv = NULL;
    nlspec.coeff[nt-1] = (*pZ)[v][0];
    nlspec.nparam += 1;

    return 0;
}

static int get_params_from_nlfunc (void)
{
    const char *s = nlspec.nlfunc;
    const char *p;
    char vname[VNAMELEN];
    int n, np = 0;
    int err = 0;

    while (*s && !err) {
	p = s;
	if (isalpha(*s) && *(s + 1)) {
	    p++;
	    n = 1;
	    while (isalnum(*p) || *p == '_') {
		p++; n++;
	    }
	    if (n > 8) return 1;
	    *vname = 0;
	    strncat(vname, s, n);
	    if (np > 0) {
		err = add_term_from_nlfunc(vname);
	    }
	    np++;
	} else p++;
	s = p;
    }

    return err;
}

static int check_for_missing_vals (void)
{
    int t, v, miss = 0;
    int t1 = nlspec.t1, t2 = nlspec.t2;
    int err;

    err = nls_auto_gen(0);
    if (err) return err;

    v = pdinfo->v - 1;

    for (t=nlspec.t1; t<=nlspec.t2; t++) {
	if (na((*pZ)[v][t])) t1++;
	else break;
    }

    for (t=nlspec.t2; t>=nlspec.t1; t--) {
	if (na((*pZ)[v][t])) t2--;
	else break;
    }

    if (t2 - t1 + 1 < nlspec.nparam) return E_DF;

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

    nlspec.t1 = t1;
    nlspec.t2 = t2;

    return 0;
}

static void update_params (const double *x)
{
    int i, v;

    for (i=0; i<nlspec.nparam; i++) {
	v = nlspec.terms[i].varnum;
	(*pZ)[v][0] = x[i];
    }
}

static int get_resid (double *fvec)
{
    int j, t, v;

    if (nls_auto_gen(0)) return 1;

    v = varindex(pdinfo, "$nls_y");
    if (v == pdinfo->v) return 1;

    nlspec.ess = 0.0;
    j = 0;
    for (t=nlspec.t1; t<=nlspec.t2; t++) {
	fvec[j] = (*pZ)[v][t];
	nlspec.ess += fvec[j] * fvec[j];
	j++;
    }
    print_iter_ess();

    return 0;
}

static int get_deriv (int i, double *deriv)
{
    int j, t, v, vec;
    char varname[VNAMELEN];

    if (nls_auto_gen(i + 1)) return 1;

    sprintf(varname, "$nls_x%d", i + 1);
    v = varindex(pdinfo, varname);
    if (v == pdinfo->v) return 1;

    vec = pdinfo->vector[v];

    j = 0;
    for (t=nlspec.t1; t<=nlspec.t2; t++) {
	if (vec) {
	    deriv[j] = - (*pZ)[v][t];
	} else {
	    deriv[j] = - (*pZ)[v][0];
	}
	j++;
    }

    return 0;
}

int nls_calc (integer *m, integer *n, double *x, double *fvec, 
	      double *fjac, integer *ldfjac, integer *iflag)
{
    int i;
    int T = *m;

    /* get current x[] values into dataset Z */
    update_params(x);

    if (*iflag == 1) {
	/* calculate functions at x and return results in fvec */
	if (get_resid(fvec)) *iflag = -1;
    }
    else if (*iflag == 2) {
	/* calculate jacobian at x and return results in fjac */
	for (i=0; i<*n; i++) {
	    if (get_deriv(i, &fjac[i*T])) *iflag = -1; 
	}	
    }
    return 0;
}

int nls_calc_approx (integer *m, integer *n, double *x, double *fvec,
		     integer *iflag)
{
    update_params(x);
    if (get_resid(fvec)) *iflag = -1;

    return 0;
}

static void add_stats_to_model (MODEL *pmod)
{
    double d, tss;
    int t, t1, t2;

    t1 = pmod->t1;
    t2 = pmod->t2;

    pmod->ess = nlspec.ess;
    pmod->sigma = sqrt(nlspec.ess/(pmod->nobs - nlspec.nparam));
    
    pmod->ybar = gretl_mean(t1, t2, (*pZ)[nlspec.depvar]);
    pmod->sdy = gretl_stddev(t1, t2, (*pZ)[nlspec.depvar]);

    tss = 0.0;
    for (t=t1; t<=t2; t++) {
	d = (*pZ)[nlspec.depvar][t] - pmod->ybar;
	tss += d * d;
    }    

    pmod->rsq = 1.0 - nlspec.ess / tss;
    pmod->adjrsq = NADBL;
}

static int add_std_errs_to_model (MODEL *pmod)
{
    int i, k;

    if (pmod->vcv == NULL && makevcv(pmod)) return E_ALLOC;

    for (i=0; i<pmod->ncoeff; i++) {
	k = ijton(i+1, i+1, pmod->ncoeff);
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

static void add_coeffs_to_model (MODEL *pmod, double *coeff)
{
    int i;

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = coeff[i];
    }
}

static int add_param_names_to_model (MODEL *pmod)
{
    int i;
    int np = pmod->ncoeff + 1;

    pmod->params = malloc(np * sizeof *pmod->params);
    if (pmod->params == NULL) return 1;

    pmod->nparams = np;

    pmod->params[0] = malloc(VNAMELEN);
    if (pmod->params[0] == NULL) {
	free(pmod->params);
	return 1;
    }
    strcpy(pmod->params[0], pdinfo->varname[nlspec.depvar]);

    for (i=1; i<=pmod->ncoeff; i++) {
	pmod->params[i] = malloc(VNAMELEN);
	if (pmod->params[i] == NULL) {
	    int j;

	    for (j=0; j<i; j++) free(pmod->params[j]);
	    free(pmod->params);
	    pmod->params = NULL;
	    pmod->nparams = 0;
	    return 1;
	}
	strcpy(pmod->params[i], nlspec.terms[i-1].name);
    }

    return 0;
}

static void add_fit_resid_to_model (MODEL *pmod, double *fvec)
{
    int t, j;

    j = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->uhat[t] = fvec[j];
	pmod->yhat[t] = (*pZ)[nlspec.depvar][t] - fvec[j];
	j++;
    }
}

static MODEL GNR (double *fvec, double *fjac)
{
    double **nZ = NULL;
    DATAINFO *ninfo;
    int *nlist;
    MODEL gnr;
    int i, j, t;
    int t1 = nlspec.t1, t2 = nlspec.t2;
    int T = t2 - t1 + 1;
    int iters = nlspec.iters;
    int err = 0;

    nlspec.t1 = 0;
    nlspec.t2 = pdinfo->n - 1;

    ninfo = create_new_dataset(&nZ, nlspec.nparam + 1, pdinfo->n, 0);
    if (ninfo == NULL) {
	gretl_model_init(&gnr, NULL);
	gnr.errcode = E_ALLOC;
	return gnr;
    }

    nlist = malloc((nlspec.nparam + 2) * sizeof *nlist);
    if (nlist == NULL) {
	free_Z(nZ, ninfo);
	free_datainfo(ninfo);
	gretl_model_init(&gnr, NULL);
	gnr.errcode = E_ALLOC;
	return gnr;
    }
    
    nlist[0] = nlspec.nparam + 1;

    for (i=0; i<=nlspec.nparam; i++) {
	nlist[i+1] = i;
	if (i == 0) {
	    /* dependent variable (NLS residual) */
	    j = 0;
	    for (t=0; t<pdinfo->n; t++) {
		if (t < t1 || t > t2) nZ[i][t] = NADBL;
		else nZ[i][t] = fvec[j++];
	    }
	} else {
	    if (nlspec.mode == ANALYTIC_DERIVS) {
		get_deriv(i-1, nZ[i]);
	    } else {
		j = T * (i - 1);
		for (t=0; t<pdinfo->n; t++) {
		    if (t < t1 || t > t2) nZ[i][t] = NADBL;
		    else nZ[i][t] = fjac[j++];
		}
	    }
	}
    }

    gnr = lsq(nlist, &nZ, ninfo, OLS, OPT_A, 0.0);

    if (gnr.errcode) {
	pputs(prn, _("In Gauss-Newton Regression:\n"));
	errmsg(gnr.errcode, prn);
	err = 1;
    } 

    if (gnr.list[0] != nlist[0]) {
	strcpy(gretl_errmsg, _("Failed to calculate Jacobian"));
	gnr.errcode = E_DATA;
    }

    if (gnr.errcode == 0) {
	add_stats_to_model(&gnr);
	if (add_std_errs_to_model(&gnr)) {
	    gnr.errcode = E_ALLOC;
	}
    }

    if (gnr.errcode == 0) {
	gretl_aic_etc(&gnr);
	gnr.ci = NLS;
	add_coeffs_to_model(&gnr, nlspec.coeff);
	add_param_names_to_model(&gnr);
	add_fit_resid_to_model(&gnr, fvec);
	gnr.list[1] = nlspec.depvar;
	gretl_model_set_int(&gnr, "iters", iters);
	gretl_model_set_double(&gnr, "tol", nlspec.tol);
    }

    nlspec.t1 = t1;
    nlspec.t2 = t2;

    free_Z(nZ, ninfo); 
    free_datainfo(ninfo);
    free(nlist);

    return gnr;
}

static void clear_nls_spec (void)
{
    int i;

    if (nlspec.terms != NULL) {
	for (i=0; i<nlspec.nparam; i++) {
	    free(nlspec.terms[i].deriv);
	}
	free(nlspec.terms);
	nlspec.terms = NULL;
    }

    free(nlspec.nlfunc);
    nlspec.nlfunc = NULL;

    free(nlspec.coeff);
    nlspec.coeff = NULL;

    nlspec.mode = NUMERIC_DERIVS;
    nlspec.nparam = 0;
    nlspec.depvar = 0;
    nlspec.iters = 0;
    nlspec.t1 = nlspec.t2 = 0;
}

static int nls_spec_start (const char *nlfunc, const DATAINFO *dinfo)
{
    char depvarname[VNAMELEN];
    const char *p;
    int v;

    /* do we already have an nls specification under way? */
    if (nlspec.nlfunc != NULL) {
	clear_nls_spec();
    }

    if (strncmp(nlfunc, "nls ", 4) == 0) { 
	p = nlfunc + 4;
    } else {
	p = nlfunc;
    }

    if (sscanf(p, "%8s = %*s", depvarname) != 1) {
	sprintf(gretl_errmsg, _("parse error in '%s'\n"), p);
	return E_PARSE;
    }

    v = varindex(dinfo, depvarname);
    if (v == dinfo->v) {
	sprintf(gretl_errmsg, _("Unknown variable '%s'"), depvarname);
	return E_UNKVAR;
    }

    nlspec.nlfunc = malloc(strlen(p) + 4);
    if (nlspec.nlfunc == NULL) return E_ALLOC;

    p = strchr(p, '=') + 1;
    while (isspace(*p)) p++;
    sprintf(nlspec.nlfunc, "%s - (%s)", depvarname, p);

    nlspec.depvar = v;
    nlspec.nparam = 0;
    nlspec.iters = 0;
    nlspec.terms = NULL;

    nlspec.t1 = dinfo->t1;
    nlspec.t2 = dinfo->t2;

    return 0;
}

static int parse_deriv_line (const char *line, int i, nls_term *term,
			     const double **Z, const DATAINFO *dinfo)
{
    int v;
    int err = 0;
    const char *p;

    term->deriv = malloc(strlen(line) - 10);
    if (term->deriv == NULL) return E_ALLOC;

    if (sscanf(line, "deriv %8s = %s", term->name, term->deriv) != 2) {
	free(term->deriv);
	term->deriv = NULL;
	fprintf(stderr, "parse error in line: '%s'\n", line);
	return E_PARSE;
    }

    p = strchr(line, '=') + 1;
    while (isspace(*p)) p++;
    strcpy(term->deriv, p);

    v = varindex(dinfo, term->name);
    if (v < dinfo->v) {
	term->varnum = v;
	nlspec.coeff[i] = Z[v][0];
    } else {
	free(term->deriv);
	term->deriv = NULL;
	sprintf(gretl_errmsg, _("Unknown variable '%s'"), term->name);
	err = E_UNKVAR;
    }	

    return err;
}

static int nls_spec_add_term (const char *line, const double **Z,
			      const DATAINFO *datainfo)
{
    nls_term *terms;
    double *coeff;
    int nt = nlspec.nparam + 1; 
    int err = 0;

    if (nlspec.nlfunc == NULL) {
	strcpy(gretl_errmsg, _("No regression function has been specified"));
	return E_PARSE;
    }
  
    terms = realloc(nlspec.terms, nt * sizeof *nlspec.terms);
    if (terms == NULL) return E_ALLOC;
    nlspec.terms = terms;

    coeff = realloc(nlspec.coeff, nt * sizeof *nlspec.coeff);
    if (coeff == NULL) {
	free(nlspec.terms);
	nlspec.terms = NULL;
	return E_ALLOC;
    }
    nlspec.coeff = coeff;

    err = parse_deriv_line(line, nt-1, &terms[nt-1], Z, datainfo);
    if (!err) {
	nlspec.nparam += 1;
	nlspec.mode = ANALYTIC_DERIVS;
    }
	
    return err;
}

int nls_parse_line (const char *line, const double **Z,
		    const DATAINFO *dinfo)
{
    if (strncmp(line, "deriv", 5) == 0) {
	return nls_spec_add_term(line, Z, dinfo);
    }
    else {
	return nls_spec_start(line, dinfo);
    }
}

static int check_derivs (integer m, integer n, double *x,
			 double *fvec, double *fjac,
			 integer ldfjac)
{
    integer mode = 1;
    integer iflag;
    doublereal *xp;
    doublereal *err;
    doublereal *fvecp;
    int i;
    int badcount = 0, zerocount = 0;

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
    nls_calc(&m, &n, x, fvec, fjac, &ldfjac, &iflag);
    if (iflag == -1) goto chkderiv_abort;
    chkder_(&m, &n, x, fvec, fjac, &ldfjac, xp, fvecp, &mode, err);

    iflag = 2;
    nls_calc(&m, &n, x, fvec, fjac, &ldfjac, &iflag);
    if (iflag == -1) goto chkderiv_abort;
    iflag = 1;
    nls_calc(&m, &n, xp, fvecp, fjac, &ldfjac, &iflag);
    if (iflag == -1) goto chkderiv_abort; 
    mode = 2;
    chkder_(&m, &n, x, fvec, fjac, &ldfjac, xp, fvecp, &mode, err);

    /* examine "err" vector */
    for (i=0; i<m; i++) {
	if (err[i] == 0.0) zerocount++;
	else if (err[i] < 0.35) badcount++;
    }

    if (zerocount > 0) {
	strcpy(gretl_errmsg, 
	       _("NLS: The supplied derivatives seem to be incorrect"));
	fprintf(stderr, "%d out of %d tests gave zero\n", zerocount, (int) m);
    }    
    else if (badcount > 0) {
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

/* Below: version of levenberg-marquandt code for use when analytical
   derivatives have been supplied */

static int lm_calculate (double *fvec, double *fjac)
{
    integer info, lwa;
    integer m, n, ldfjac;
    integer *ipvt;
    doublereal *wa;
    int err = 0;

    m = nlspec.t2 - nlspec.t1 + 1; /* number of observations */
    n = nlspec.nparam;             /* number of parameters */
    lwa = 5 * n + m;               /* work array size */
    ldfjac = m;                    /* leading dimension of fjac array */

    wa = malloc(lwa * sizeof *wa);
    ipvt = malloc(n * sizeof *ipvt);

    if (wa == NULL || ipvt == NULL) {
	err = E_ALLOC;
	goto nls_cleanup;
    }

    err = check_derivs(m, n, nlspec.coeff, fvec, fjac, ldfjac);
    if (err) goto nls_cleanup; 

    lmder1_(nls_calc, &m, &n, nlspec.coeff, fvec, fjac, &ldfjac, &nlspec.tol, 
	    &info, ipvt, wa, &lwa);

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
	pprintf(prn, _("Convergence achieved after %d iterations\n"),
		nlspec.iters);
	break;
    case 4:
    case 5:
    case 6:
    case 7:
	sprintf(gretl_errmsg, 
		_("NLS: failed to converge after %d iterations"),
		nlspec.iters);
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

/* Below: version of levenberg-marquandt code for use when the Jacobian
   must be approximated numerically */

static int lm_approximate (double *fvec, double *fjac)
{
    integer info, m, n, ldfjac;
    integer maxfev, mode = 1, nprint = 0, nfev = 0;
    integer iflag = 0;
    integer *ipvt;
    doublereal gtol = 0.0;
    doublereal epsfcn = 0.0, factor = 100.;
    doublereal *diag, *qtf;
    doublereal *wa1, *wa2, *wa3, *wa4;
    int err = 0;
    
    m = nlspec.t2 - nlspec.t1 + 1; /* number of observations */
    n = nlspec.nparam;             /* number of parameters */
    ldfjac = m;                    /* leading dimension of fjac array */

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

    lmdif_(nls_calc_approx, &m, &n, nlspec.coeff, fvec, 
	   &nlspec.tol, &nlspec.tol, &gtol, &maxfev, &epsfcn, diag, &mode, &factor,
	   &nprint, &info, &nfev, fjac, &ldfjac, 
	   ipvt, qtf, wa1, wa2, wa3, wa4);

    nlspec.iters = nfev;

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
		nlspec.iters);
	break;
    case 5:
    case 6:
    case 7:
    case 8:
	sprintf(gretl_errmsg, 
		_("NLS: failed to converge after %d iterations"),
		nlspec.iters);
	err = 1;
	break;
    default:
	break;
    }

    if (!err) {
	double ess = nlspec.ess;
	int iters = nlspec.iters;

	fdjac2_(nls_calc_approx, &m, &n, nlspec.coeff, fvec, fjac, 
		&ldfjac, &iflag, &epsfcn, wa4);
	nlspec.ess = ess;
	nlspec.iters = iters;
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

MODEL nls (double ***mainZ, DATAINFO *maininfo, PRN *mainprn)
{
    MODEL nlsmod;
    double *fvec, *fjac;
    int origv = maininfo->v;
    int err = 0;

    gretl_model_init(&nlsmod, maininfo);

    if (nlspec.nlfunc == NULL) {
	strcpy(gretl_errmsg, _("No regression function has been specified"));
	nlsmod.errcode = E_PARSE;
	return nlsmod;
    }   

    /* "export" these as file-scope globals */
    pZ = mainZ;
    pdinfo = maininfo;
    prn = mainprn; 

    if (nlspec.mode == NUMERIC_DERIVS) {
	err = get_params_from_nlfunc();
	if (err) {
	    clear_nls_spec();
	    if (err == 1) {
		nlsmod.errcode = E_PARSE;
	    } else {
		nlsmod.errcode = err;
	    }
	    return nlsmod;
	}
    }

    if (nlspec.nparam == 0) {
	strcpy(gretl_errmsg, _("No regression function has been specified"));
	clear_nls_spec();
	nlsmod.errcode = E_PARSE;
	return nlsmod;
    }   

    fvec = malloc(pdinfo->n * sizeof *fvec);
    fjac = malloc(pdinfo->n * nlspec.nparam * sizeof *fvec);
    if (fvec == NULL || fjac == NULL) {
	free(fvec);
	free(fjac);
	clear_nls_spec();
	nlsmod.errcode = E_ALLOC;
	return nlsmod;
    }

    nlsmod.errcode = err = check_for_missing_vals();

    if (!err) {
	if (toler > 0) nlspec.tol = toler;
#if 0
	else nlspec.tol = sqrt(dpmpar_(&one));
#else
	else nlspec.tol = pow(dpmpar_(&one), .75);
#endif
	if (nlspec.mode == NUMERIC_DERIVS) {
	    pputs(prn, _("Using numerical derivatives\n"));
	    err = lm_approximate(fvec, fjac);
	} else {
	    pputs(prn, _("Using analytical derivatives\n"));
	    err = lm_calculate(fvec, fjac);
	}
	pprintf(prn, _("Tolerance = %g\n"), nlspec.tol);
    }

    if (!err) {
	nlsmod = GNR(fvec, fjac);
    } else {
	if (nlsmod.errcode == 0) { 
	    nlsmod.errcode = E_UNSPEC;
	}
    }

    free(fvec);
    free(fjac);
    clear_nls_spec();

    dataset_drop_vars(pdinfo->v - origv, pZ, pdinfo);

    *mainZ = *pZ;
  
    return nlsmod;
}

void set_nls_toler (double x)
{
    toler = x;
}


