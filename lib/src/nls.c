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
#include "internal.h"
#include "../../plugin/f2c.h"
#include "../../minpack/minpack.h"  

typedef struct _nls_spec nls_spec;
typedef struct _nls_term nls_term;

struct _nls_term {
    char name[9];       /* name of parameter */
    char *deriv;        /* string representation of derivative */
    int varnum;         /* ID number of the var holding the derivative */
};

struct _nls_spec {
    int depvar;         /* ID number of dependent variable */
    char *nlfunc;       /* string representation of regression function */
    int nparam;         /* number of parameters to be estimated */
    int iters;          /* number of iterations performed */
    int t1;             /* starting observation */
    int t2;             /* ending observation */
    double ess;         /* error sum of squares */
    nls_term *terms;    /* array of info on terms in the function */
    doublereal *coeff;  /* coefficient estimates */
};

/* file-scope global variables */
static double ***pZ;
static DATAINFO *pdinfo;
static PRN *prn;
static nls_spec nlspec;

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
    char varname[9];

    if (nls_auto_gen(i + 1)) return 1;

    sprintf(varname, "$nls_x%d", i + 1);
    v = varindex(pdinfo, varname);
    if (v == pdinfo->v) return 1;

    vec = pdinfo->vector[v];

    j = 0;
    for (t=nlspec.t1; t<=nlspec.t2; t++) {
	if (vec) deriv[j] = - (*pZ)[v][t];
	else deriv[j] = - (*pZ)[v][0];
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

static void add_stats_to_model (MODEL *pmod)
{
    double d, tss;
    int t, t1, t2;

    t1 = pmod->t1;
    t2 = pmod->t2;

    pmod->ess = nlspec.ess;
    pmod->sigma = sqrt(nlspec.ess/(pmod->nobs - nlspec.nparam));
    
    pmod->ybar = _esl_mean(t1, t2, (*pZ)[nlspec.depvar]);
    pmod->sdy = _esl_stddev(t1, t2, (*pZ)[nlspec.depvar]);

    tss = 0.0;
    for (t=t1; t<=t2; t++) {
	d = (*pZ)[nlspec.depvar][t] - pmod->ybar;
	tss += d * d;
    }    
    pmod->rsq = 1.0 - nlspec.ess / tss;
    pmod->adjrsq = NADBL;
}

static void add_coeffs_to_model (MODEL *pmod, double *coeff)
{
    int i;

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i+1] = coeff[i];
    }
}

static int add_param_names_to_model (MODEL *pmod)
{
    int i;

    pmod->params = malloc((1 + pmod->ncoeff) * sizeof *pmod->params);
    if (pmod->params == NULL) return 1;

    pmod->params[0] = malloc(9);
    if (pmod->params[0] == NULL) {
	free(pmod->params);
	return 1;
    }
    strcpy(pmod->params[0], pdinfo->varname[nlspec.depvar]);

    for (i=1; i<=pmod->ncoeff; i++) {
	pmod->params[i] = malloc(9);
	if (pmod->params[i] == NULL) {
	    int j;

	    for (j=0; j<i; j++) free(pmod->params[j]);
	    free(pmod->params);
	    pmod->params = NULL;
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

static MODEL GNR (double *fvec)
{
    double **nZ = NULL;
    DATAINFO *ninfo;
    int *nlist;
    MODEL gnr;
    int i, j, t;
    int t1 = nlspec.t1, t2 = nlspec.t2;
    int err = 0;

    nlspec.t1 = 0;
    nlspec.t2 = pdinfo->n - 1;

    ninfo = create_new_dataset(&nZ, nlspec.nparam + 1, pdinfo->n, 0);
    if (ninfo == NULL) {
	_init_model(&gnr, NULL);
	gnr.errcode = E_ALLOC;
	return gnr;
    }

    nlist = malloc((nlspec.nparam + 2) * sizeof *nlist);
    if (nlist == NULL) {
	free_Z(nZ, ninfo);
	free_datainfo(ninfo);
	_init_model(&gnr, NULL);
	gnr.errcode = E_ALLOC;
	return gnr;
    }
    
    nlist[0] = nlspec.nparam + 1;

    for (i=0; i<=nlspec.nparam; i++) {
	nlist[i+1] = i;
	if (i == 0) {
	    j = 0;
	    for (t=0; t<pdinfo->n; t++) {
		if (t < t1 || t > t2) nZ[i][t] = NADBL;
		else nZ[i][t] = fvec[j++];
	    }
	} else {
	    get_deriv(i-1, nZ[i]);
	}
    }

    gnr = lsq(nlist, &nZ, ninfo, OLS, 0, 0.0);
    if (gnr.errcode) {
	pputs(prn, _("In Gauss-Newton Regression:\n"));
	errmsg(gnr.errcode, prn);
	err = 1;
    } else {
	add_coeffs_to_model(&gnr, nlspec.coeff);
	add_param_names_to_model(&gnr);
	add_fit_resid_to_model(&gnr, fvec);
	gnr.list[1] = nlspec.depvar;
    }

    nlspec.t1 = t1;
    nlspec.t2 = t2;

    free_Z(nZ, ninfo); 
    free_datainfo(ninfo);
    free(nlist);

    return gnr;
}

static int nls_spec_start (const char *nlfunc, const DATAINFO *dinfo)
{
    char depvarname[9];
    const char *p;
    int v;

    /* do we already have an nls specification under way? */
    if (nlspec.nparam > 0) return E_PARSE;

    if (strncmp(nlfunc, "nls ", 4) == 0) 
	p = nlfunc + 4;
    else 
	p = nlfunc;

    if (sscanf(p, "%8s = %*s", depvarname) != 1) return E_PARSE;

    v = varindex(dinfo, depvarname);
    if (v == dinfo->v) return E_UNKVAR;

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

    coeff = realloc(nlspec.coeff, nt * sizeof *nlspec.coeff);
    if (coeff == NULL) {
	free(terms);
	return E_ALLOC;
    }

    nlspec.coeff = coeff;

    err = parse_deriv_line(line, nt-1, &terms[nt-1], Z, datainfo);
    if (!err) {
	nlspec.terms = terms;
	nlspec.nparam += 1;
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

static void clear_nls_spec (void)
{
    int i;

    for (i=0; i<nlspec.nparam; i++) {
	free(nlspec.terms[i].deriv);
    }

    free(nlspec.terms);
    nlspec.terms = NULL;

    free(nlspec.nlfunc);
    nlspec.nlfunc = NULL;

    free(nlspec.coeff);
    nlspec.coeff = NULL;

    nlspec.nparam = 0;
    nlspec.depvar = 0;
    nlspec.iters = 0;
    nlspec.t1 = nlspec.t2 = 0;
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

    /* examine err vector */
    for (i=0; i<m; i++) {
	if (err[i] == 0.0) zerocount++;
	else if (err[i] < 0.35) badcount++;
    }

    if (zerocount > 0) {
	strcpy(gretl_errmsg, 
	       _("NLS: The supplied derivatives seem to be incorrect"));
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

    return zerocount;
}

static int lm_calculate (double *fvec, double toler)
{
    integer info, lwa;
    integer m, n, ldfjac;
    integer *ipvt;
    doublereal tol;
    doublereal *fjac, *wa;
    int err = 0;

    m = nlspec.t2 - nlspec.t1 + 1; /* number of observations */
    n = nlspec.nparam;             /* number of parameters */
    lwa = 5 * n + m;               /* work array size */
    ldfjac = m;                    /* leading dimension of fjac array */
    tol = toler;                   /* tolerance for convergence */

    wa = malloc(lwa * sizeof *wa);
    ipvt = malloc(n * sizeof *ipvt);
    fjac = malloc(m * n * sizeof *fjac);

    if (wa == NULL || ipvt == NULL || fjac == NULL) {
	err = E_ALLOC;
	goto nls_cleanup;
    }

    err = check_derivs(m, n, nlspec.coeff, fvec, fjac, ldfjac);
    if (err) goto nls_cleanup; 

    /* run levenberg-marquandt nonlinear least squares from minpack */
    lmder1_(nls_calc, &m, &n, nlspec.coeff, fvec, fjac, &ldfjac, &tol, 
	    &info, ipvt, wa, &lwa);

    switch ((int) info) {
    case -1: 
	err = 1;
	break;
    case 0:
	pputs(prn, _("Invalid NLS specification"));
	pputs(prn, "\n");
	err = 1;
	break;
    case 1:
    case 2:
    case 3:
	pprintf(prn, _("NLS: convergence achieved after %d iterations\n"),
		nlspec.iters);
	break;
    case 4:
    case 5:
    case 6:
    case 7:
	pputs(prn, _("NLS failed to converge\n"));
	err = 1;
	break;
    default:
	break;
    }

 nls_cleanup:
    free(wa);
    free(ipvt);
    free(fjac);

    return err;    
}

#if 0
static void print_nls_spec (void)
{
    int i;

    printf("NLS specification:\n");
    printf("residual: %s\n", nlspec.nlfunc);
    for (i=0; i<nlspec.nparam; i++) {
	printf("deriv[%d] (var number %d): %s = %s\n", i+1, 
	       (nlspec.terms[i]).varnum,
	       (nlspec.terms[i]).name,
	       (nlspec.terms[i]).deriv);
    }
}
#endif

MODEL nls (double ***mainZ, DATAINFO *maininfo, PRN *mainprn)
{
    MODEL nlsmod;
    double *fvec;
    double toler = 0.0001; /* make this configurable */
    int origv = maininfo->v;
    int err = 0;

    _init_model(&nlsmod, maininfo);

    if (nlspec.nlfunc == NULL) {
	sprintf(gretl_errmsg, _("No regression function has been specified"));
	nlsmod.errcode = E_PARSE;
	return nlsmod;
    }    

    if (nlspec.nparam == 0) {
	sprintf(gretl_errmsg, _("No derivatives have been specified"));
	nlsmod.errcode = E_PARSE;
	return nlsmod;
    }

    fvec = malloc(maininfo->n * sizeof *fvec);
    if (fvec == NULL) {
	nlsmod.errcode = E_ALLOC;
	return nlsmod;
    }

    /* "export" these as file-scope globals */
    pZ = mainZ;
    pdinfo = maininfo;
    prn = mainprn;

#if 0
    print_nls_spec();
#endif

    nlsmod.errcode = err = check_for_missing_vals();

    if (!err) err = lm_calculate(fvec, toler);
    if (!err) nlsmod = GNR(fvec);
    if (!err) {
	add_stats_to_model(&nlsmod);  
	_aicetc(&nlsmod);
	nlsmod.ci = NLS;
    } 

    free(fvec);
    clear_nls_spec();

    dataset_drop_vars(pdinfo->v - origv, pZ, pdinfo);

    *mainZ = *pZ;
  
    return nlsmod;
}


