/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* panel data plugin for gretl */

#include "libgretl.h"
#include "f2c.h"
#include "clapack_double.h"
#include "gretl_model.h"

#undef PDEBUG

typedef struct hausman_t_ hausman_t;

struct hausman_t_ {
    int ns;
    double sigma_e;
    double H;
    double *bdiff;
    double *sigma;
};

/* .................................................................. */

static void haus_init (hausman_t *haus)
{
    haus->ns = 0;
    haus->sigma_e = NADBL;
    haus->H = NADBL;
    haus->bdiff = NULL;
    haus->sigma = NULL;
}

static void haus_free (hausman_t *haus)
{
    free(haus->bdiff);
    free(haus->sigma);
}

static int haus_alloc (hausman_t *haus, int ns)
{
    haus->ns = ns;

    haus->bdiff = malloc(ns * sizeof *haus->bdiff);
    if (haus->bdiff == NULL) {
	return E_ALLOC;
    }

    haus->sigma = malloc(((ns * ns + ns) / 2) * sizeof *haus->sigma);
    if (haus->sigma == NULL) {
	free(haus->bdiff);
	haus->bdiff = NULL;
	return E_ALLOC; 
    }

    return 0;
}   

/* .................................................................. */

static void print_panel_coeff (MODEL *pmod, MODEL *panelmod,
			       DATAINFO *pdinfo, int i, 
			       PRN *prn)
{
    char numstr[18];

    sprintf(numstr, "(%.5g)", panelmod->sderr[i]);
    pprintf(prn, "%*s: %14.5g %15s\n", VNAMELEN,
	    pdinfo->varname[pmod->list[i+2]],
	    panelmod->coeff[i], numstr);
}

/* .................................................................. */

static double group_means_variance (MODEL *pmod, 
				    double **Z, DATAINFO *pdinfo,
				    double ***groupZ, DATAINFO **ginfo,
				    int nunits, int effn, int *unit_obs,
				    int T, int effT) 
{
    int i, j, k, t, s, start, *list;
    double xx;
    MODEL meanmod;

#ifdef PDEBUG
    fprintf(stderr, "group_means_variance: creating dataset\n"
	    " nvar=%d, nobs=%d, *groupZ=%p\n", pmod->list[0],
	    effn, (void *) *groupZ);
#endif

    *ginfo = create_new_dataset(groupZ, pmod->list[0], effn, 0);
    if (*ginfo == NULL) {
	return NADBL;
    }

    list = malloc((pmod->list[0] + 1) * sizeof *list);
    if (list == NULL) {
	clear_datainfo(*ginfo, CLEAR_FULL);
	free(*ginfo);
	return NADBL;
    }

#ifdef PDEBUG
    fprintf(stderr, "gmv: *groupZ=%p\n", (void *) *groupZ);
#endif

    list[0] = pmod->list[0];
    k = 1;
    for (j=1; j<=list[0]; j++) { 
	/* the variables */
	if (pmod->list[j] == 0) {
	    list[j] = 0;
	    continue;
	}
	list[j] = k;
	start = 0;
	s = 0;
	for (i=0; i<nunits; i++) { 
	    /* the observations */
	    if (unit_obs[i] == 0) {
		if (pdinfo->time_series == STACKED_TIME_SERIES) {
		    start += T;
		} else {
		    start++;
		}
		continue;
	    }
	    xx = 0.0;
	    if (pdinfo->time_series == STACKED_TIME_SERIES) {
		for (t=start; t<start+T; t++) {
		    if (!na(pmod->uhat[t])) {
			xx += Z[pmod->list[j]][t];
		    }
		}
		start += T;
	    } else {
		for (t=start; t<pdinfo->n; t += nunits) {
		    if (!na(pmod->uhat[t])) {
			xx += Z[pmod->list[j]][t];
		    }
		}
		start++;
	    }
	    (*groupZ)[k][s] = xx / (double) effT;
#ifdef PDEBUG
	    fprintf(stderr, "Set groupZ[%d][%d] = %g\n", k, i, (*groupZ)[k][i]);
#endif
	    s++;
	}
	k++;
    }

#ifdef PDEBUG
    fprintf(stderr, "group_means_variance: about to run OLS\n");
    fprintf(stderr, "*groupZ=%p, ginfo=%p\n", (void *)*groupZ, (void *)*ginfo);
#endif

    meanmod = lsq(list, groupZ, *ginfo, OLS, OPT_A, 0.0);

#ifdef PDEBUG
    fprintf(stderr, "gmv: lsq errcode was %d\n", meanmod.errcode);
    fprintf(stderr, "meanmod.sigma = %g\n", meanmod.sigma);
#endif

    if (meanmod.errcode) {
	xx = NADBL;
    } else {
	xx = meanmod.sigma * meanmod.sigma;
    }

    clear_model(&meanmod);
    free(list);

#ifdef PDEBUG
    fprintf(stderr, "gmv: done freeing stuff\n");
#endif

    return xx;
}

/* .................................................................. */

static void 
vcv_slopes (hausman_t *haus, MODEL *pmod, int nunits, int subt)
{
    int i, j, k, idx;

    k = 0;
    for (i=0; i<haus->ns; i++) {
	for (j=i; j<haus->ns; j++) {
	    idx = ijton(i+1, j+1, pmod->ncoeff);
#ifdef PDEBUG
	    fprintf(stderr, "setting sigma[%d] using vcv[%d] (%d, %d) = %g\n",
		    k, idx, i+2, j+2, pmod->vcv[idx]);
#endif
	    if (subt) {
		haus->sigma[k++] -= pmod->vcv[idx];
	    } else {
		haus->sigma[k++] = pmod->vcv[idx];
	    }
	}
    }
}

/* .................................................................. */

static int bXb (hausman_t *haus)
{
    char uplo = 'L'; 
    integer nrhs = 1;
    integer n = haus->ns;
    integer ldb = n;
    integer info = 0;
    integer *ipiv;
    double *x;
    int i;

    /* make a copy of haus->bdiff first */
    x = malloc(n * sizeof *x);
    if (x == NULL) {
	return E_ALLOC;
    }

    ipiv = malloc(n * sizeof *ipiv);
    if (ipiv == NULL) {
	free(x);
	return E_ALLOC;
    }

    for (i=0; i<haus->ns; i++) {
	x[i] = haus->bdiff[i];
    }

    /* solve for X-inverse * b */
    dspsv_(&uplo, &n, &nrhs, haus->sigma, ipiv, x, &ldb, &info);

    if (info > 0) {
	fprintf(stderr, "Hausman sigma matrix is singular\n");
    } else if (info < 0) {
	fprintf(stderr, "Illegal entry in Hausman sigma matrix\n");
    } else {
	haus->H = 0.0;
	for (i=0; i<haus->ns; i++) {
	    haus->H += x[i] * haus->bdiff[i];
	}
    }

    free(x);
    free(ipiv);

    return info;
}

/* .................................................................. */

static double LSDV (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		    int nunits, int *unit_obs, int T,
		    hausman_t *haus, PRN *prn) 
{
    int i, j, t, start, oldv = pdinfo->v;
    int dvlen, ndum = 0;
    int *dvlist;
    double var, F;
    MODEL lsdv;

    for (i=0; i<nunits; i++) {
	if (unit_obs[i] > 1) {
	    ndum++;
	}
    }

    dvlen = pmod->list[0] + ndum;
    ndum--; 

#ifdef PDEBUG
    fprintf(stderr, "ndum = %d, dvlen = %d\n", ndum, dvlen);
#endif

    /* We can be assured there's an intercept in the original
       regression */

    dvlist = malloc(dvlen * sizeof *dvlist);
    if (dvlist == NULL) {
	return NADBL;
    }

    if (dataset_add_vars(ndum, pZ, pdinfo)) {
	free(dvlist);
	return NADBL;
    }

    start = 0;
    j = 0;
    for (i=0; i<nunits; i++) {
	int dv = oldv + j;

	if (unit_obs[i] < 2) {
	    if (pdinfo->time_series == STACKED_TIME_SERIES) {
		start += T;
	    } else {
		start++;
	    }
	    continue;
	}
	sprintf(pdinfo->varname[dv], "unit_%d", i + 1);
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[dv][t] = 0.0;
	}
	if (pdinfo->time_series == STACKED_TIME_SERIES) {
	    for (t=start; t<start+T; t++) {
		(*pZ)[dv][t] = 1.0;
	    }
	    start += T;
	} else {
	    for (t=start; t<pdinfo->n; t += nunits) {
		(*pZ)[dv][t] = 1.0;
	    }
	    start++;
	}
	j++;
	if (j == ndum) {
	    break;
	}
    }

    dvlist[0] = dvlen - 1;

    for (i=1; i<=pmod->list[0]; i++) {
	dvlist[i] = pmod->list[i];
    }
    for (i=0; i<ndum; i++) {
	dvlist[pmod->list[0] + i + 1] = oldv + i;
    }

#ifdef PDEBUG
    printlist(dvlist, "dvlist");
#endif

    lsdv = lsq(dvlist, pZ, pdinfo, OLS, OPT_A, 0.0);

#ifdef PDEBUG
    printmodel(&lsdv, pdinfo, OPT_NONE, prn);
#endif

    if (lsdv.errcode) {
	var = NADBL;
	pputs(prn, _("Error estimating fixed effects model\n"));
	errmsg(lsdv.errcode, prn);
    } else {
	var = lsdv.sigma * lsdv.sigma;
	pputs(prn, 
	      _("                          Fixed effects estimator\n"
		"          allows for differing intercepts by cross-sectional "
		"unit\n"
		"         (slope standard errors in parentheses, a_i = "
		"intercepts)\n\n"));

	/* skip the constant in this printing */
	for (i=1; i<pmod->list[0] - 1; i++) {
	    print_panel_coeff(&lsdv, &lsdv, pdinfo, i, prn);
	    if (haus != NULL) {
		haus->bdiff[i-1] = lsdv.coeff[i];
	    }
	}

	j = 0;
	for (i=0; i<nunits; i++) {
	    /* print per-unit intercept estimates */
	    char dumstr[VNAMELEN];
	    double x;

	    if (unit_obs[i] < 2) {
		continue;
	    }

	    if (j == ndum) {
		x = lsdv.coeff[0];
	    } else {
		x = lsdv.coeff[j + pmod->list[0] - 1] + lsdv.coeff[0];
	    }
	    sprintf(dumstr, "a_%d", i + 1);
	    pprintf(prn, "%*s: %14.4g\n", VNAMELEN, dumstr, x);
	    j++;
	}

	pprintf(prn, _("\nResidual variance: %g/(%d - %d) = %g\n"), 
		lsdv.ess, lsdv.nobs, lsdv.ncoeff, var);
	F = (pmod->ess - lsdv.ess) * lsdv.dfd / (lsdv.ess * ndum);
	pprintf(prn, _("Joint significance of unit dummy variables:\n"
		       " F(%d, %d) = %g with p-value %g\n"), ndum,
		lsdv.dfd, F, fdist(F, ndum, lsdv.dfd));
	pputs(prn, _("(A low p-value counts against the null hypothesis that "
		     "the pooled OLS model\nis adequate, in favor of the fixed "
		     "effects alternative.)\n\n"));

	if (haus != NULL) {
	    makevcv(&lsdv);
	    haus->sigma_e = lsdv.sigma;
	    vcv_slopes(haus, &lsdv, nunits, 0);
	}
    }

    clear_model(&lsdv);
    dataset_drop_vars(pdinfo->v - oldv, pZ, pdinfo);
    free(dvlist);

    return var;
}

/* .................................................................. */

static int random_effects (MODEL *pmod, double **Z, DATAINFO *pdinfo, 
			   double **groupZ, double theta, 
			   int nunits, int effn, int *unit_obs, int T, 
			   hausman_t *haus, PRN *prn)
{
    double **reZ;
    DATAINFO *reinfo;
    MODEL remod;
    int *relist;
    int re_n = T * effn;
    int i, j, k, t, bigt;
    int err = 0;

    /* regression list */
    relist = malloc((pmod->list[0] + 1) * sizeof *relist);
    if (relist == NULL) {
	return E_ALLOC;
    }

    reinfo = create_new_dataset(&reZ, pmod->list[0], re_n, 0);
    if (reinfo == NULL) {
	free(relist);
	return E_ALLOC;
    }

#ifdef PDEBUG
    fprintf(stderr, "reZ: series length = T * effn = %d * %d = %d\n",
	    T, effn, T * effn);
#endif

    relist[0] = pmod->list[0];

    /* create transformed variables: original data minus theta
       times the appropriate group or unit mean 
    */

    k = 1;
    for (i=1; i<=relist[0]; i++) {
	const double *xi = Z[pmod->list[i]];
	double *gm = groupZ[k];
	int u = 0;

	if (pmod->list[i] == 0) {
	    relist[i] = 0;
	    continue;
	}

	relist[i] = k;

	if (pdinfo->time_series == STACKED_TIME_SERIES) {
	    for (j=0; j<nunits; j++) {
		if (unit_obs[j] == 0) {
		    continue;
		}
		for (t=0; t<T; t++) {
		    int rt = u * T + t;

		    bigt = j * T + t;
		    if (na(pmod->uhat[bigt])) {
			reZ[k][rt] = NADBL;
		    } else {
			reZ[k][rt] = xi[bigt] - theta * gm[u];
		    }
#ifdef PDEBUG
		    fprintf(stderr, "set reZ[%d][%d]=Z[%d][%d]-theta*"
			    "groupZ[%d][%d]=%g\n", k, rt, 
			    pmod->list[i], bigt,
			    k, u, reZ[k][rt]);
#endif
		}
		u++;
	    }
	} else {
	    int rt = 0;

	    for (t=0; t<T; t++) {
		u = 0;
		for (j=0; j<nunits; j++) {
		    bigt = t * nunits + j;
		    if (unit_obs[j] == 0) {
			continue;
		    }
		    if (na(pmod->uhat[bigt])) {
			reZ[k][rt] = NADBL;
		    } else {
			reZ[k][rt] = xi[bigt] - theta * gm[u];
		    }
#ifdef PDEBUG
		    fprintf(stderr, "set reZ[%d][%d]=Z[%d][%d]-theta*"
			    "groupZ[%d][%d]=%g\n", k, rt, 
			    pmod->list[i], bigt,
			    k, u, reZ[k][rt]);
#endif
		    u++;
		    rt++;
		}
	    }
	}
	k++;
    }

    for (t=0; t<re_n; t++) {
	reZ[0][t] -= theta;
    }

#ifdef PDEBUG
    fprintf(stderr, "random_effects: about to run OLS\n");
#endif

    remod = lsq(relist, &reZ, reinfo, OLS, OPT_A, 0.0);
    if ((err = remod.errcode)) {
	pputs(prn, _("Error estimating random effects model\n"));
	errmsg(err, prn);
    } else {
	pputs(prn,
	      _("                         Random effects estimator\n"
		"           allows for a unit-specific component to the "
		"error term\n"
		"                     (standard errors in parentheses)\n\n"));
	for (i=0; i<relist[0] - 1; i++) {
	    print_panel_coeff(pmod, &remod, pdinfo, i, prn);
	    if (i > 0) {
		haus->bdiff[i-1] -= remod.coeff[i];
	    }
	}
	makevcv(&remod);
	vcv_slopes(haus, &remod, nunits, 1);
    }

    clear_model(&remod);

    free_Z(reZ, reinfo);
    clear_datainfo(reinfo, CLEAR_FULL);
    free(reinfo);

    free(relist);    

    return err;
}

/* .................................................................. */

static void 
unit_error_variances (double *uvar, const MODEL *pmod, const DATAINFO *pdinfo,
		      int nunits, int T, const int *unit_obs)
{
    int i, t, start = 0;

    for (i=0; i<nunits; i++) {
	uvar[i] = 0.0;
	if (pdinfo->time_series == STACKED_TIME_SERIES) {
	    for (t=start; t<start+T; t++) {
		if (!na(pmod->uhat[t])) {
		    uvar[i] += pmod->uhat[t] * pmod->uhat[t];
		}
	    }
	    start += T;
	} else {
	    for (t=start; t<pdinfo->n; t += nunits) {
		if (!na(pmod->uhat[t])) {
		    uvar[i] += pmod->uhat[t] * pmod->uhat[t];
		}
	    }
	    start++;
	}
	if (unit_obs[i] > 1) {
	    uvar[i] /= (double) unit_obs[i]; 
	}
    }
}

/* .................................................................. */

static int breusch_pagan_LM (const MODEL *pmod, const DATAINFO *pdinfo, 
			     int nunits, const int *unit_obs,
			     int T, int effT, PRN *prn)
{
    double *ubar, LM, eprime = 0.0;
    int i, t, start = 0;

    ubar = malloc(nunits * sizeof *ubar);
    if (ubar == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<nunits; i++) {
	ubar[i] = 0.0;
	if (pdinfo->time_series == STACKED_TIME_SERIES) {
	    for (t=start; t<start+T; t++) {
		if (!na(pmod->uhat[t])) {
		    ubar[i] += pmod->uhat[t];
		}
	    }
	    start += T;
	} else {
	    for (t=start; t<pdinfo->n; t += nunits) {
		if (!na(pmod->uhat[t])) {
		    ubar[i] += pmod->uhat[t];
		}
	    }
	    start++;
	}
	ubar[i] /= (double) unit_obs[i]; 
	eprime += ubar[i] * ubar[i];
    }

#ifdef PDEBUG
    fprintf(stderr,  "breusch_pagan: found ubars\n");
#endif

    pputs(prn, _("\nMeans of pooled OLS residuals for cross-sectional "
		 "units:\n\n"));

    for (i=0; i<nunits; i++) {
	if (unit_obs[i] == 0) {
	    continue;
	}
	pprintf(prn, _(" unit %2d: %13.5g\n"), i + 1, ubar[i]);
    }

    free(ubar);

    LM = (double) pmod->nobs / (2.0 * (effT - 1.0)) * 
	pow((effT * effT * eprime / pmod->ess) - 1.0, 2);

    pprintf(prn, _("\nBreusch-Pagan test statistic:\n"
		   " LM = %g with p-value = prob(chi-square(1) > %g) = %g\n"), 
	    LM, LM, chisq(LM, 1));
    pputs(prn, _("(A low p-value counts against the null hypothesis that "
		 "the pooled OLS model\nis adequate, in favor of the random "
		 "effects alternative.)\n\n"));

    return 0;
}

/* .................................................................. */

static int do_hausman_test (hausman_t *haus, PRN *prn)
{
#ifdef PDEBUG
    int i, ns = haus->ns;
    int nterms = (ns * ns + ns) / 2;

    for (i=0; i<ns; i++) {
 	fprintf(stderr, "b%d_FE - beta%d_RE = %g\n", i, i, haus->bdiff[i]);
    }
    fputc('\n', stderr);

    for (i=0; i<nterms; i++) {
 	fprintf(stderr, "vcv_diff[%d] = %g\n", i, haus->sigma[i]);
    }
#endif

    if (bXb(haus)) { 
	pputs(prn, _("Error attempting to invert vcv difference matrix\n"));
	return 1;
    }

    if (haus->H < 0) {
	pputs(prn, _("\nHausman test matrix is not positive definite (this "
		     "result may be treated as\n\"fail to reject\" the random effects "
		     "specification).\n"));
    } else {
	pprintf(prn, _("\nHausman test statistic:\n"
		       " H = %g with p-value = prob(chi-square(%d) > %g) = %g\n"),
		haus->H, haus->ns, haus->H, chisq(haus->H, haus->ns));
	pputs(prn, _("(A low p-value counts against the null hypothesis that "
		     "the random effects\nmodel is consistent, in favor of the fixed "
		     "effects model.)\n"));
    }

    return 0;
}

int get_maj_min (const DATAINFO *pdinfo, int *maj, int *min)
{
    int startmaj, startmin;
    int endmaj, endmin;

    if (sscanf(pdinfo->stobs, "%d:%d", &startmaj, &startmin) != 2) {
	return 1;
    }

    if (sscanf(pdinfo->endobs, "%d:%d", &endmaj, &endmin) != 2) {
	return 1;
    } 

    *maj = endmaj - startmaj + 1;
    *min = endmin - startmin + 1;

    return 0;
}

int n_included_units (const MODEL *pmod, const DATAINFO *pdinfo,
		      int *unit_obs)
{
    int nmaj, nmin;
    int k, ninc = 0;

    if (get_maj_min(pdinfo, &nmaj, &nmin)) {
	return -1;
    }

    if (pdinfo->time_series == STACKED_TIME_SERIES) {
	int nunits = nmaj;
	int nperiods = nmin;
	int i, j;

	for (i=0; i<nunits; i++) {
	    unit_obs[i] = 0;
	    for (j=0; j<nperiods; j++) {
		k = i * nperiods + j;
		if (!na(pmod->uhat[k])) {
		    unit_obs[i] += 1;
		}
	    }
	    if (unit_obs[i] > 0) {
		ninc++;
	    }
	}
    } else {
	/* stacked cross sections */
	int nunits = nmin;
	int nperiods = nmaj;
	int i, j;

	for (i=0; i<nunits; i++) {
	    unit_obs[i] = 0;
	    for (j=0; j<nperiods; j++) {
		k = j * nunits + i;
		if (!na(pmod->uhat[k])) {
		    unit_obs[i] += 1;
		}
	    }
	    if (unit_obs[i] > 0) {
		ninc++;
	    }	    
	}
    }

    return ninc;
}

static int effective_T (int *unit_obs, int nunits)
{
    int i, effT = 0;

    for (i=0; i<nunits; i++) {
	if (unit_obs[i] > effT) {
	    effT = unit_obs[i];
	}
    }

    return effT;
}

/* .................................................................. */

int panel_diagnostics (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		       PRN *prn)
{
    int unbal = gretl_model_get_int(pmod, "unbalanced");
    int nunits, ns, T;
    int effn, effT;
    int *unit_obs = NULL;
    double var1, var2;
    hausman_t haus;
    int err = 0;

#ifdef PDEBUG
    fputs("\n*** Starting panel_diagnostics ***\n", stderr);
#endif

    if (pmod->ifc == 0) {
	return 1;
    }

    if (get_panel_structure(pdinfo, &nunits, &T)) {
	return E_DATA;
    }

    haus_init(&haus);

    unit_obs = malloc(nunits * sizeof *unit_obs);
    if (unit_obs == NULL) {
	return E_ALLOC;
    }

    effn = n_included_units(pmod, pdinfo, unit_obs);
    if (effn < nunits) {
	fprintf(stderr, "number of units included = %d\n", effn);
    }
    effT = effective_T(unit_obs, nunits);

#ifdef PDEBUG
    fprintf(stderr, "nunits=%d, T=%d, effn=%d, effT=%d\n",
	    nunits, T, effn, effT);
#endif

    if (!unbal && effn > pmod->ncoeff) {
	ns = pmod->ncoeff - 1;
	err = haus_alloc(&haus, ns);
	if (err) {
	    goto bailout;
	}
    }   

    if (!unbal) {
	pprintf(prn, _("      Diagnostics: assuming a balanced panel with %d "
		       "cross-sectional units\n "
		       "                        observed over %d periods\n\n"), 
		effn, effT);
    }

    var2 = LSDV(pmod, pZ, pdinfo, nunits, unit_obs, T, 
		(unbal)? NULL : &haus, prn);

#ifdef PDEBUG
    fprintf(stderr, "panel_diagnostics: LSDV gave variance = %g\n", var2);
#endif

    if (unbal) {
	pprintf(prn, "Omitting random effects model since "
		"panel is unbalanced\n");
	goto bailout;
    }

    breusch_pagan_LM(pmod, pdinfo, nunits, unit_obs, T, effT, prn);

#ifdef PDEBUG
    fprintf(stderr, "panel_diagnostics: done breusch_pagan_LM()\n");
#endif
    
    if (effn > pmod->ncoeff && !na(var2)) {
	double **groupZ = NULL;
	DATAINFO *ginfo = NULL;
	double theta;
	
	var1 = group_means_variance(pmod, *pZ, pdinfo, &groupZ, &ginfo, 
				    nunits, effn, unit_obs, T, effT);

	if (na(var1)) { 
	    pputs(prn, _("Couldn't estimate group means regression\n"));
	} else {
	    pprintf(prn, _("Residual variance for group means "
			   "regression: %g\n\n"), var1);    
	    theta = 1.0 - sqrt(var2 / (effT * var1));
	    random_effects(pmod, *pZ, pdinfo, groupZ, theta, nunits, 
			   effn, unit_obs, T, &haus, prn);
	    do_hausman_test(&haus, prn);
	}

	free_Z(groupZ, ginfo);
	clear_datainfo(ginfo, CLEAR_FULL);
	free(ginfo);
    }

 bailout:

    free(unit_obs);
    haus_free(&haus);

    return err;
}

/* .................................................................. */

static int 
write_uvar_to_dataset (double *uvar, int nunits, int T,
		       double **Z, DATAINFO *pdinfo)
{
    int i, t, bigt;
    int uv = pdinfo->v - 1;

    if (pdinfo->time_series == STACKED_TIME_SERIES) {
	for (i=0; i<nunits; i++) {
	    for (t=0; t<T; t++) {
		bigt = i * T + t;
		if (uvar[i] <= 0.0) {
		    Z[uv][bigt] = 0.0;
		} else {
		    Z[uv][bigt] = 1.0 / sqrt(uvar[i]);
		}
	    }
	}
    } else {
	for (t=0; t<T; t++) {
	    for (i=0; i<nunits; i++) {
		bigt = t * nunits + i;
		if (uvar[i] <= 0.0) {
		    Z[uv][bigt] = 0.0;
		} else {
		    Z[uv][bigt] = 1.0 / sqrt(uvar[i]);
		}
	    }
	}
    }

    return 0;
}

static int
allocate_weight_var (double ***pZ, DATAINFO *pdinfo)
{
    if (dataset_add_vars(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    strcpy(pdinfo->varname[pdinfo->v - 1], "unit_wt");

    return 0;
}

static double max_coeff_diff (const MODEL *pmod, const double *bvec)
{
    double diff, maxdiff = 0.0;
    int i;

    for (i=0; i<pmod->ncoeff; i++) {
	diff = fabs(pmod->coeff[i] - bvec[i]);
	if (diff > maxdiff) {
	    maxdiff = diff;
	}
    }

    return maxdiff;
}

static int
groupwise_hetero_test (MODEL *pmod, double s2, const double *uvar, 
		       int nunits, const int *unit_obs)
{
    GRETLTEST test;
    double x2, s2h = 0.0;
    int i, df = 0;

    for (i=0; i<nunits; i++) {
	int uT = unit_obs[i];

	if (uT > 0) {
	    s2h += uT * log(uvar[i]);
	    df++;
	}
    }

    x2 = pmod->nobs * log(s2) - s2h;
    df--;

    gretl_test_init(&test);
    strcpy(test.type, 
	   N_("Likelihood ratio test for groupwise heteroskedasticity"));
    strcpy(test.h_0, N_("the units have a common error variance"));
    test.teststat = GRETL_TEST_LR;
    test.dfn = df;
    test.value = x2;
    test.pvalue = chisq(x2, df);

    return add_test_to_model(pmod, &test);
}

#define LN_2_PI 1.837877066409345

static double pooled_ll (const MODEL *pmod)
{
    double n = pmod->nobs;

    return -(n / 2.0) * (1.0 + LN_2_PI - log(n) + log(pmod->ess));
}

static double real_ll (const MODEL *pmod, const double *uvar, 
		       int nunits, const int *unit_obs)
{
    double ll = -(pmod->nobs / 2.0) * LN_2_PI;
    int i;

    for (i=0; i<nunits; i++) {
	if (unit_obs[i] > 0) {
	    ll -= (unit_obs[i] / 2.0) * (1.0 + log(uvar[i]));
	}
    }

    return ll;
}

/* we can't estimate a group-specific variance based on just one
   observation */

static int singleton_check (const int *unit_obs, int nunits)
{
    int i;

    for (i=0; i<nunits; i++) {
	if (unit_obs[i] == 1) {
	    return 1;
	}
    }

    return 0;
}

#define SMALLDIFF 0.0001
#define WLS_MAX   20

MODEL panel_wls_by_unit (int *list, double ***pZ, DATAINFO *pdinfo,
			 gretlopt opt, PRN *prn)
{
    MODEL mdl;
    gretlopt wlsopt = OPT_A;
    double *uvar = NULL;
    double *bvec = NULL;
    double s2, diff = 1.0;
    int *unit_obs = NULL;
    int *wlist = NULL;
    int nunits, effn, T, effT;
    int orig_v = pdinfo->v;
    int i, iter = 0;

    gretl_errmsg_clear();

    if (opt & OPT_T) {
	/* iterating: no degrees-of-freedom correction */
	wlsopt |= OPT_N; 
    } 

    gretl_model_init(&mdl);

    if (get_panel_structure(pdinfo, &nunits, &T)) {
	mdl.errcode = E_DATA;
	return mdl;
    }

    unit_obs = malloc(nunits * sizeof *unit_obs);
    if (unit_obs == NULL) {
	mdl.errcode = E_ALLOC;
	return mdl;
    }

    uvar = malloc(nunits * sizeof *uvar);
    if (unit_obs == NULL) {
	free(unit_obs);
	mdl.errcode = E_ALLOC;
	return mdl;
    }    
    
    mdl = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
    if (mdl.errcode) {
	goto bailout;
    }

    effn = n_included_units(&mdl, pdinfo, unit_obs);
    effT = effective_T(unit_obs, nunits);  

    if (opt & OPT_T) {
	if (singleton_check(unit_obs, nunits)) {
	    gretl_errmsg_set(_("Can't produce ML estimates: "
			       "some units have only one observation"));
	    mdl.errcode = E_DF;
	    goto bailout;
	}
    }

    s2 = mdl.ess / mdl.nobs;

    if ((opt & OPT_V) && (opt & OPT_T)) {
	pprintf(prn, "\nOLS error variance = %g\n", s2);
	pprintf(prn, "log-likelihood = %g\n", pooled_ll(&mdl));
    }

    if (allocate_weight_var(pZ, pdinfo)) {
	mdl.errcode = E_ALLOC;
	goto bailout;
    }

    if (opt & OPT_T) {
	bvec = malloc(mdl.ncoeff * sizeof *bvec);
	if (bvec == NULL) {
	    mdl.errcode = E_ALLOC;
	    goto bailout;
	}
    }

    /* allocate and construct WLS regression list */
    wlist = malloc((mdl.list[0] + 2) * sizeof *wlist);
    if (wlist == NULL) {
	mdl.errcode = E_ALLOC;
	goto bailout;
    }
    wlist[0] = mdl.list[0] + 1;
    wlist[1] = pdinfo->v - 1; /* weight variable: the last var added */
    for (i=2; i<=wlist[0]; i++) {
	wlist[i] = mdl.list[i-1];
    }

    /* if wanted, iterate to ML solution; otherwise just do
       one-step FGLS estimation 
    */

    while (diff > SMALLDIFF) {

	iter++;

	unit_error_variances(uvar, &mdl, pdinfo, nunits, T, unit_obs);

	if (opt & OPT_V) {
	    if (opt & OPT_T) {
		pprintf(prn, "\n*** %s %d ***\n", _("iteration"), 
			iter);
	    } else {
		pputc(prn, '\n');
	    }
	    pputs(prn, " unit    variance\n");
	    for (i=0; i<nunits; i++) {
		if (unit_obs[i] > 0) {
		    pprintf(prn, "%5d%12g\n", i + 1, uvar[i]);
		}
	    }
	}

	write_uvar_to_dataset(uvar, nunits, T, *pZ, pdinfo);

	if (opt & OPT_T) {
	    /* save coefficients for comparison */
	    for (i=0; i<mdl.ncoeff; i++) {
		bvec[i] = mdl.coeff[i];
	    }
	}

	clear_model(&mdl);

	mdl = lsq(wlist, pZ, pdinfo, WLS, wlsopt, 0.0);

	if (mdl.errcode || iter > WLS_MAX) {
	    mdl.errcode = E_NOCONV;
	    break;
	}

	if (opt & OPT_T) {
	    diff = max_coeff_diff(&mdl, bvec);
	    if (opt & OPT_V && iter == 1) {
		pprintf(prn, "\nFGLS pooled error variance = %g\n",
			mdl.ess / mdl.nobs);
	    }
	} else {
	    /* one-step FGLS */
	    break;
	} 
    }

    if (!mdl.errcode) {
	set_model_id(&mdl);
	gretl_model_set_int(&mdl, "n_included_units", effn);
	gretl_model_set_int(&mdl, "unit_weights", 1);
	mdl.nwt = 0;
	if (opt & OPT_T) {
	    gretl_model_set_int(&mdl, "iters", iter);
	    groupwise_hetero_test(&mdl, s2, uvar, nunits, unit_obs);
	    unit_error_variances(uvar, &mdl, pdinfo, nunits, T, unit_obs);
	    mdl.lnL = real_ll(&mdl, uvar, nunits, unit_obs);
	    if (opt & OPT_V) {
		pputc(prn, '\n');
	    }
	}
    }    

 bailout:

    free(unit_obs);
    free(uvar);
    free(wlist);
    free(bvec);

    dataset_drop_vars(pdinfo->v - orig_v, pZ, pdinfo);
    
    return mdl;
}

/* .................................................................. */

static void panel_copy_var (double **targZ, DATAINFO *targinfo, int targv,
			    double *src, DATAINFO *srcinfo, int srcv,
			    int order)
{
    int t, j = 0;

    for (t=srcinfo->t1; t<=srcinfo->t2; t++) {
	if (t % srcinfo->pd >= order) {
	    targZ[targv][j++] = src[t];
	} 
    }

    if (srcv == -1) {
	strcpy(targinfo->varname[targv], "uhat");
	strcpy(VARLABEL(targinfo, targv), _("residual"));
    } else {
	strcpy(targinfo->varname[targv], srcinfo->varname[srcv]);
	strcpy(VARLABEL(targinfo, targv), VARLABEL(srcinfo, srcv));
    }
}

static void make_reduced_data_info (DATAINFO *targ, DATAINFO *src, int order)
{
    targ->pd = src->pd - order;
    ntodate(targ->stobs, src->t1 + order, src);
    targ->sd0 = obs_str_to_double(targ->stobs); 
    targ->time_series = src->time_series;
}

static void panel_lag (double **tmpZ, DATAINFO *tmpinfo, 
		       double *src, DATAINFO *srcinfo, 
		       int v, int order, int lag)
{
    int t, j = 0;

    for (t=srcinfo->t1; t<=srcinfo->t2; t++) {
	if (t % srcinfo->pd >= order) {
	    tmpZ[v][j++] = src[t - lag];
	}
    }

    sprintf(tmpinfo->varname[v], "uhat_%d", lag);
    *VARLABEL(tmpinfo, v) = 0;
}

/* - do some sanity checks
   - create a local copy of the required portion of the data set,
     skipping the obs that will be missing
   - copy in the lags of uhat
   - estimate the aux model
   - destroy the temporary data set
*/

int panel_autocorr_test (MODEL *pmod, int order, 
			 double **Z, DATAINFO *pdinfo, 
			 PRN *prn, GRETLTEST *test)
{
    int *aclist;
    double **tmpZ;
    DATAINFO *tmpinfo;
    MODEL aux;
    double trsq, LMF;
    int i, nv, nunits, nobs, err = 0;
    int sn = pdinfo->t2 - pdinfo->t1 + 1;

    /* basic checks */
    if (order <= 0) order = 1;
    if (order > pdinfo->pd - 1) return E_DF;
    if (pmod->ncoeff + order >= sn) return E_DF;

    if (pdinfo->time_series != STACKED_TIME_SERIES ||
	!balanced_panel(pdinfo)) { 
        return E_DATA;
    }

    if (pmod->missmask != NULL) {
	return E_MISSDATA;
    }

    /* get number of cross-sectional units */
    nunits = sn / pdinfo->pd;

    /* we lose "order" observations for each unit */
    nobs = sn - nunits * order;

    /* the required number of variables */
    nv = pmod->list[0] + order;

    /* create temporary reduced dataset */
    tmpinfo = create_new_dataset(&tmpZ, nv, nobs, 0);
    if (tmpinfo == NULL) return E_ALLOC;

    make_reduced_data_info(tmpinfo, pdinfo, order);

#ifdef PDEBUG
    fprintf(stderr, "Created data set, n=%d, pd=%d, vars=%d, stobs='%s'\n", 
	    tmpinfo->n, tmpinfo->pd, tmpinfo->v, tmpinfo->stobs);
#endif

    /* allocate the auxiliary regression list */
    aclist = malloc((nv + 1) * sizeof *aclist);
    if (aclist == NULL) {
	err = E_ALLOC;
    } 

    if (!err) {
	int k;

	aclist[0] = pmod->list[0] + order;
	/* copy model uhat to position 1 in temp data set */
	aclist[1] = 1;
	panel_copy_var(tmpZ, tmpinfo, 1,
		       &pmod->uhat[0], pdinfo, -1,
		       order);
	/* copy across the original indep vars, making
	   the new regression list while we're at it */
	k = 2;
	for (i=2; i<=pmod->list[0]; i++) {
	    if (pmod->list[i] == 0) { /* the constant */
		aclist[i] = 0;
	    } else {
		aclist[i] = k;
		panel_copy_var(tmpZ, tmpinfo, k, 
			       &Z[pmod->list[i]][0], pdinfo, pmod->list[i], 
			       order);
		k++;
	    }
	}
    }

    if (!err) {
	int v = pmod->list[0] - 1;

	/* add lags of uhat to temp data set */
	for (i=1; i<=order; i++) {
	    panel_lag(tmpZ, tmpinfo, &pmod->uhat[0], 
		      pdinfo, v + i, order, i);
	    aclist[v + i + 1] = v + i;
	}
    }

    if (!err) {
	aux = lsq(aclist, &tmpZ, tmpinfo, OLS, OPT_A, 0.0);
	err = aux.errcode;
	if (err) {
	    errmsg(aux.errcode, prn);
	}
    } 

    if (!err) {
	aux.aux = AUX_AR;
	aux.order = order;
	printmodel(&aux, tmpinfo, OPT_NONE, prn);
	trsq = aux.rsq * aux.nobs;
	LMF = (aux.rsq/(1.0 - aux.rsq)) * 
	    (aux.nobs - pmod->ncoeff - order)/order; 

	pprintf(prn, "\n%s: LMF = %f,\n", _("Test statistic"), LMF);
	pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"), 
		order, aux.nobs - pmod->ncoeff - order, LMF,
		fdist(LMF, order, aux.nobs - pmod->ncoeff - order));

	pprintf(prn, "\n%s: TR^2 = %f,\n", 
		_("Alternative statistic"), trsq);
	pprintf(prn, "%s = P(%s(%d) > %g) = %.3g\n\n", 	_("with p-value"), 
		_("Chi-square"), order, trsq, chisq(trsq, order));

	if (test != NULL) {
	    strcpy(test->type, N_("LM test for autocorrelation up to order %s"));
	    strcpy(test->h_0, N_("no autocorrelation"));
	    sprintf(test->param, "%d", order);
	    test->teststat = GRETL_TEST_LMF;
	    test->dfn = order;
	    test->dfd = aux.nobs - pmod->ncoeff - order;
	    test->value = LMF;
	    test->pvalue = fdist(LMF, test->dfn, test->dfd);
	}
    }

    free(aclist);
    clear_model(&aux); 

    free_Z(tmpZ, tmpinfo);
    clear_datainfo(tmpinfo, CLEAR_FULL);
    free(tmpinfo);

    return err;
}

int switch_panel_orientation (double **Z, DATAINFO *pdinfo)
{
    double **tmpZ;
    int i, j, k, t, nvec;
    int nunits = pdinfo->pd;
    int nperiods = pdinfo->n / nunits;
    double pdx;
    char **markers = NULL;

    tmpZ = malloc((pdinfo->v - 1) * sizeof *tmpZ);
    if (tmpZ == NULL) {
	return E_ALLOC;
    }

    /* allocate temporary data matrix */
    j = 0;
    for (i=1; i<pdinfo->v; i++) {
	if (pdinfo->vector[i]) {
	    tmpZ[j] = malloc(pdinfo->n * sizeof **tmpZ);
	    if (tmpZ[j] == NULL) {
		for (i=0; i<j; i++) {
		    free(tmpZ[i]);
		}
		free(tmpZ);
		return E_ALLOC;
	    }
	    j++;
	} 
    }
    nvec = j;

    /* allocate marker space if relevant */
    if (pdinfo->S != NULL) {
	markers = malloc(pdinfo->n * sizeof *markers);
	if (markers != NULL) {
	    for (t=0; t<pdinfo->n; t++) {
		markers[t] = malloc(OBSLEN);
		if (markers[t] == NULL) {
		    free(markers);
		    markers = NULL;
		    break;
		} else {
		    strcpy(markers[t], pdinfo->S[t]);
		}
	    }
	}
    }

    /* copy the data (vectors) across */
    j = 0;
    for (i=1; i<pdinfo->v; i++) {
	if (pdinfo->vector[i]) {
	    for (t=0; t<pdinfo->n; t++) {
		tmpZ[j][t] = Z[i][t];
	    }
	    j++;
	} 
    }

    /* copy the data back in transformed order: construct a set of
       time series for each unit in turn -- and do markers if present */
    for (k=0; k<nunits; k++) {
	j = 0;
	for (i=1; i<pdinfo->v; i++) {
	    if (!pdinfo->vector[i]) continue;
	    for (t=0; t<nperiods; t++) {
		Z[i][k * nperiods + t] = tmpZ[j][k + nunits * t];
	    }
	    j++;
	}
	if (markers != NULL) {
	    for (t=0; t<nperiods; t++) {
		strcpy(pdinfo->S[k * nperiods + t], markers[k + nunits * t]);
	    }
	}
    }

    /* change the datainfo setup */
    pdinfo->time_series = STACKED_TIME_SERIES;
    pdinfo->pd = nperiods;

    pdinfo->sd0 = 1.0;
    pdx = 0.1;
    while (nperiods /= 10) {
	pdx *= 0.1;
    }
    pdinfo->sd0 += pdx;

    ntodate(pdinfo->stobs, 0, pdinfo);
    ntodate(pdinfo->endobs, pdinfo->n - 1, pdinfo);

    /* clean up */

    for (i=0; i<nvec; i++) {
	free(tmpZ[i]);
    }
    free(tmpZ);

    if (markers != NULL) {
	for (t=0; t<pdinfo->n; t++) {
	    free(markers[t]);
	}
	free(markers);
    }

    return 0;
}



