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
    char name[9];
    char *deriv;
    int varnum;
};

struct _nls_spec {
    int depvar;
    char *nlfunc;
    int nparam;
    int iters;
    double ess;
    nls_term *terms;
    doublereal *coeff;
};

/* In fortran arrays, column entries are contiguous.
   Columns of data matrix X hold variables, rows hold observations.
   So in a fortran array, entries for a given variable are
   contiguous.
*/

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

/* functions below need to be generalized to deal with an arbitrary
   user-specified regression function */

static void update_params (const double *x)
{
    int i, v;

    for (i=0; i<nlspec.nparam; i++) {
	v = nlspec.terms[i].varnum;
	(*pZ)[v][0] = x[i];
    }
}

static int get_resid (double *fvec, int T)
{
    int t, v;

    if (nls_auto_gen(0)) return 1;

    v = varindex(pdinfo, "$nls_y");
    if (v == pdinfo->v) return 1;

    nlspec.ess = 0.0;
    for (t=0; t<T; t++) {
	fvec[t] = (*pZ)[v][t];
	nlspec.ess += fvec[t] * fvec[t];
    }
    print_iter_ess();

    return 0;
}

static int get_deriv (int i, double *deriv, int T)
{
    int t, v, vec;
    char varname[9];

    if (nls_auto_gen(i + 1)) return 1;

    sprintf(varname, "$nls_x%d", i + 1);
    v = varindex(pdinfo, varname);
    if (v == pdinfo->v) return 1;

    vec = pdinfo->vector[v];
    
    for (t=0; t<T; t++) {
	if (vec) deriv[t] = - (*pZ)[v][t];
	else deriv[t] = - (*pZ)[v][0];
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
	if (get_resid(fvec, T)) *iflag = -1;
    }
    else if (*iflag == 2) {
	/* calculate jacobian at x and return results in fjac */
	for (i=0; i<*n; i++) {
	    if (get_deriv(i, &fjac[i*T], T)) *iflag = -1; 
	}	
    }
    return 0;
}

static void add_stats_to_model (MODEL *pmod, int nparam)
{
    double d, tss;
    int t, t1, t2;

    t1 = pmod->t1;
    t2 = pmod->t2;

    pmod->ess = nlspec.ess;
    pmod->sigma = sqrt(nlspec.ess/(pmod->nobs - nparam));
    
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

static MODEL GNR (int nobs, int nparam, double *nle, double *coeff)
{
    double **nZ = NULL;
    DATAINFO *ninfo;
    int *nlist;
    MODEL gnr;
    int i, t;
    int err = 0;

    ninfo = create_new_dataset(&nZ, nparam + 1, nobs, 0);
    nlist = malloc((nparam + 2) * sizeof *nlist);
    
    nlist[0] = nparam + 1;

    for (i=0; i<=nparam; i++) {
	nlist[i+1] = i;
	if (i == 0) {
	    for (t=0; t<nobs; t++) {
		nZ[i][t] = nle[t];
	    }
	} else {	
	    get_deriv(i-1, nZ[i], nobs);
	}
    }

    gnr = lsq(nlist, &nZ, ninfo, OLS, 0, 0.0);
    if (gnr.errcode) {
	pputs(prn, "In Gauss-Newton Regression:\n");
	errmsg(gnr.errcode, prn);
	err = 1;
    } else {
	add_coeffs_to_model(&gnr, coeff);
	add_param_names_to_model(&gnr);
    }

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

    if (nlspec.nparam > 0) return 1;

    nlspec.depvar = -1;

    if (strncmp(nlfunc, "nls ", 4) == 0) 
	p = nlfunc + 4;
    else 
	p = nlfunc;

    if (sscanf(p, "%8s = %*s", depvarname) != 1) return 1;

    v = varindex(dinfo, depvarname);
    if (v == dinfo->v) return 1;

    nlspec.nlfunc = malloc(strlen(p) + 4);
    if (nlspec.nlfunc == NULL) return 1;

    p = strchr(p, '=') + 1;
    while (isspace(*p)) p++;
    sprintf(nlspec.nlfunc, "%s - (%s)", depvarname, p);

    nlspec.depvar = v;
    nlspec.nparam = 0;
    nlspec.iters = 0;
    nlspec.terms = NULL;

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
    nlspec.depvar = -1;
    nlspec.iters = 0;
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
    int i, badcount = 0;

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
	if (err[i] < 0.35) badcount++;
    }

    if (badcount > 0) {
	pputs(prn, "Warning: The supplied derivatives may be incorrect,\n"
	      "or perhaps the data matrix is ill-conditioned.\n");
	pprintf(prn, "%d out of %d gradients looked suspicious.\n\n",
		badcount, (int) m);
    }

 chkderiv_abort:
    free(xp);
    free(err);
    free(fvecp);

    return 0;
}

static int lm_calculate (int nobs, int nparam, double *x, double *fvec,
			 double toler)
{
    integer info, lwa;
    integer m, n, ldfjac;
    integer *ipvt;
    doublereal tol;
    doublereal *fjac, *wa;
    int err = 0;

    m = nobs;           /* number of observations */
    n = nparam;         /* number of parameters */
    lwa = 5 * n + m;    /* work array size */
    ldfjac = m;         /* leading dimension of fjac array */
    tol = toler;        /* tolerance for convergence */

    wa = malloc(lwa * sizeof *wa);
    ipvt = malloc(n * sizeof *ipvt);
    fjac = malloc(m * n * sizeof *fjac);

    if (wa == NULL || ipvt == NULL || fjac == NULL) {
	err = 1;
	goto nls_cleanup;
    }

    check_derivs(m, n, x, fvec, fjac, ldfjac);

    /* run levenberg-marquandt from minpack */
    lmder1_(nls_calc, &m, &n, x, fvec, fjac, &ldfjac, &tol, 
	    &info, ipvt, wa, &lwa);

    switch ((int) info) {
    case -1: 
	pputs(prn, "genr failed\n");
	err = 1;
	break;
    case 0:
	pputs(prn, "Invalid input for NLS\n");
	err = 1;
	break;
    case 1:
    case 2:
    case 3:
	pprintf(prn, "NLS: convergence achieved after %d iterations\n",
		nlspec.iters);
	break;
    case 4:
    case 5:
    case 6:
    case 7:
	pputs(prn, "NLS failed to converge\n");
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

    /* "export" these as local globals ;-) */
    pZ = mainZ;
    pdinfo = maininfo;
    prn = mainprn;

#if 0
    print_nls_spec();
#endif

    err = lm_calculate(pdinfo->n, nlspec.nparam, nlspec.coeff, fvec, toler); 

    if (err) {
	nlsmod.errcode = E_UNSPEC;
    } else {
	nlsmod = GNR(pdinfo->n, nlspec.nparam, fvec, nlspec.coeff);
    }

    if (!err) {
	add_stats_to_model(&nlsmod, nlspec.nparam);  
	_aicetc(&nlsmod);
	nlsmod.ci = NLS;
    } 

    free(fvec);
    clear_nls_spec();

    dataset_drop_vars(pdinfo->v - origv, pZ, pdinfo);

    *mainZ = *pZ;
  
    return nlsmod;
}


