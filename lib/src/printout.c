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
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

/*  printout.c - simple text print routines for some gretl structs */ 

#include "libgretl.h"
#include "internal.h"
#include "version.h"
#include <time.h>

#ifdef OS_WIN32
#define isnan(x) ((x) != (x))
#endif

static void print_float_10 (const double x, print_t *prn);
static void print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			 const int c, print_t *prn);
static void _depvarstats (const MODEL *pmod, print_t *prn);
static int _essline (const MODEL *pmod, print_t *prn, int wt);
static void _rsqline (const MODEL *pmod, print_t *prn);
static void _Fline (const MODEL *pmod, print_t *prn);
static void _dwline (const MODEL *pmod, print_t *prn);
static void print_discrete_stats (const MODEL *pmod, 
				  const DATAINFO *pdinfo, 
				  print_t *prn);
static void print_coeff_interval (const DATAINFO *pdinfo, const MODEL *pmod, 
				  const int c, const double t, print_t *prn);
void _putxx (const double xx);
void _mxout (const double *rr, const int *list, const int ci,
	     const DATAINFO *pdinfo, const int batch, print_t *prn);


/* ......................................................... */ 

static void noconst (print_t *prn)
{
    pprintf(prn, "The model has no constant term.\n"  
	    "F is calculated as in Sect. 4.4 of Ramanathan's Introductory "
	    "Econometrics.\n"
	    "R-squared is the square of the correlation between the "
	    "observed and fitted\n values of the dependent variable.\n\n");
}

/* ......................................................... */ 

static void _depvarstats (const MODEL *pmod, print_t *prn)
{
    pprintf(prn, "Mean of dep. var. %17.3f  S.D. of dep. variable %17.3f\n", 
	    pmod->ybar, pmod->sdy);
}

/* ........................................................... */

static void printxx (const double xx, char *str, const int ci)
{
#define BIG 1.0e+9
#define SMALL 1.0e-6
    static char word[32];
    double xf, xxabs;
    long xi;
    int d1, d2;

    switch (ci) {
    case PRINT:
	d1 = 8;  
	d2 = 6;
	break;
    case STORE:
	d1 = 10;
	d2 = 12;
	break;
    case SUMMARY:
	d1 = d2 = 6;
	break;
    default:
	d1 = d2 = 8;
	break;
    }

    if (xx < BIG) {
	xi = (long) xx;
	xf = xx - xi;
    }
    else xf = 0.5;

    xxabs = fabs(xx);
    if (floateq(xf, 0.0)) sprintf(word, "%.0f", xx); 
    else if (xxabs < SMALL) sprintf(word, "%.4g", xx);
    else if (xxabs < 1.0) sprintf(word, "%.*g", d1, xx);
    else sprintf(word, "%.*g", d2, xx);
/*      if (haschar('.', word) >= 0 && strlen(word) < (d2 + 1)) */
/*  	strcat(word, "0"); */
    strcpy(str, word);	
}	

/* ......................................................... */ 

static int _essline (const MODEL *pmod, print_t *prn, int wt)
{
    if ((wt && pmod->ess_wt < 0) || (!wt && pmod->ess < 0)) {
	pprintf(prn, "Error sum of squares (%g) is not > 0\n\n", 
		(wt)? pmod->ess_wt : pmod->ess);
	return 1;
    }

    pprintf(prn, "Error Sum of Sq (ESS) ");
    space(3, prn);
    print_float_10(wt? pmod->ess_wt : pmod->ess, prn);
    pprintf(prn, "  Std Err of Resid. (sgmahat) ");
    space(1, prn);
    print_float_10(wt? pmod->sigma_wt : pmod->sigma, prn);
    pprintf(prn, "\n");
    return 0;
}

/* ......................................................... */ 

static void _rsqline (const MODEL *pmod, print_t *prn)
{
    double xx = pmod->rsq;

    if (pmod->rsq > .999 && pmod->rsq < .999999) xx = .999;

    pprintf(prn, "Unadjusted R-squared %14.3f  Adjusted R-squared ", xx);
    if (na(pmod->adjrsq))
	pprintf(prn, "%21s", "undefined\n");
    else {
	xx = pmod->adjrsq;
	if (pmod->adjrsq > .999 && pmod->rsq < .999999) xx = .999;
	pprintf(prn, "%20.3f\n", xx);
    }
}

/* ......................................................... */ 

static void _Fline (const MODEL *pmod, print_t *prn)
{
    char tmp[32];

    sprintf(tmp, "F-statistic (%d, %d)", pmod->dfn, pmod->dfd);
    pprintf(prn, "%s", tmp);
    space(24 - strlen(tmp), prn);
    if (na(pmod->fstt))
	pprintf(prn, "%11s  p-value for F() %23s\n", "undefined", "undefined");
    else pprintf(prn, "%11g  p-value for F() %23f\n", pmod->fstt,
	   fdist(pmod->fstt, pmod->dfn, pmod->dfd));
}

/* ......................................................... */ 

static void _dwline (const MODEL *pmod, print_t *prn)
{
    if (na(pmod->dw))
    pprintf(prn, "Durbin-Watson stat. %15s  First-order autocorr. "
	    "coeff %11s\n", "undefined", "undefined");
    else 
	pprintf(prn, "Durbin-Watson stat. %15.3f  First-order autocorr. "
		"coeff %11.3f\n", 
		pmod->dw, pmod->rho);
}

/* ......................................................... */ 

static void dhline (const MODEL *pmod, print_t *prn)
{
    double sderr, h = 0.0;
    int i = pmod->ldepvar, T = pmod->nobs - 1;
    char hstring[20];

    sderr = pmod->sderr[i-1];
    if ((T * sderr * sderr) > 1.0) 
	strcpy(hstring, "         undefined");
    else {
	h = pmod->rho * sqrt(T/(1 - T * sderr * sderr));
	sprintf(hstring, "    %14.3f", h);
    }
    pprintf(prn, "Durbin's h stat. %s  First-order autocorr. coeff %11.3f\n", 
	   hstring, pmod->rho);
    if (floatneq(h, 0.0)) 
	pprintf(prn, "(Using variable %d for h stat, with T' = %d)\n", 
	       pmod->list[i], T);
}

/* ......................................................... */ 

static int _pmax (const MODEL *pmod)
{
    int i, k = 0;
    double tstat, tmin = 4.0;
    
    for (i=1; i <= pmod->ncoeff - pmod->ifc; i++) {
	tstat = fabs(pmod->coeff[i] / pmod->sderr[i]);
	if (tstat < tmin) {
	    tmin = tstat;
	    k = i;
	}
    }
    if (tprob(tmin, pmod->dfd) > .10) return pmod->list[k+1];
    return 0;
}

/* ......................................................... */ 

static void _pmax_line (const MODEL *pmod, const DATAINFO *pdinfo, 
			print_t *prn)
{
    int k = pmod->ncoeff - pmod->ifc;

    if (k < 2) return;
    if ((k = _pmax(pmod)))
	pprintf(prn, "Excluding the constant, p-value was highest "
		"for variable %d (%s)\n\n", k, pdinfo->varname[k]);
}

/* ......................................................... */ 

static void covhdr (print_t *prn)
{
    pprintf(prn, "\nCOVARIANCE MATRIX OF REGRESSION COEFFICIENTS\n\n");
}

/* ......................................................... */ 

void _putxx (const double xx)
{
    if (xx < 0.0001) puts("< 0.0001");
    else printf("%g\n", xx);
}

/* ......................................................... */

void session_time (void)
{
    time_t runtime = time(NULL);

    printf("Current session: %s", ctime(&runtime));
}

/* ......................................................... */

void logo (void)
{
    printf("gretl client, for library version %s,\n", version_string);
    puts("copyright Ramu Ramanathan and Allin Cottrell.");
    puts("This is free software with ABSOLUTELY NO WARRANTY.");
}

/* ......................................................... */

void gui_logo (void)
{
    printf("gretl: gui client for gretl version %s,\n", version_string);
    puts("copyright Allin Cottrell.");
    puts("This is free software with ABSOLUTELY NO WARRANTY.");
}

/* ......................................................... */

void print_model_confints (const MODEL *pmod, const DATAINFO *pdinfo, 
			   print_t *prn)
{
    int i, ncoeff = pmod->list[0];
    double t = tcrit95(pmod->dfd);

    pprintf(prn, "t(%d, .025) = %.3f\n\n", pmod->dfd, t);
    pprintf(prn, "      VARIABLE      COEFFICIENT      95%% CONFIDENCE "
	    "INTERVAL\n\n");      

    if (pmod->ifc) {
	print_coeff_interval(pdinfo, pmod, ncoeff, t, prn);
	ncoeff--;
    }
    for (i=2; i<=ncoeff; i++) 
	print_coeff_interval(pdinfo, pmod, i, t, prn);
    pprintf(prn, "\n");
}

/* ......................................................... */

static void print_model_tests (const MODEL *pmod, print_t *prn)
{
    int i;

    for (i=0; i<pmod->ntests; i++) {
	pprintf(prn, "%s -\n"
		"  Null hypothesis: %s\n"
		"  Test statistic: %s\n"
		"  with p-value = %s\n\n",
		(pmod->tests[i]).type, (pmod->tests[i]).h_0, 
		(pmod->tests[i]).teststat, (pmod->tests[i]).pvalue);
    }
}

/* ......................................................... */ 

void printmodel (const MODEL *pmod, const DATAINFO *pdinfo, print_t *prn)
{
    int i, ncoeff;
    char startdate[8];
    char enddate[8];
    int t1 = pmod->t1, t2 = pmod->t2;

    if (pmod->ci == CORC || pmod->ci == HILU) t1 += 1;

    ncoeff = pmod->list[0];
    ntodate(startdate, t1, pdinfo);
    ntodate(enddate, t2, pdinfo);

    switch(pmod->aux) {
    case AUX_AR:
	pprintf(prn, "\nTest for autocorrelation up to the periodicity\n");
	break;	
    case AUX_ARCH:
	pprintf(prn, "\nTest for ARCH of order %d\n", 
		pmod->list[0] - 2);
	break;	
    case AUX_SQ:
	pprintf(prn, "\nAuxiliary regression for non-linearity test "
		"(squared terms)\n");
	break;
    case AUX_LOG:
	pprintf(prn, "\nAuxiliary regression for non-linearity test "
		"(log terms)\n");
	break;	
    case AUX_WHITE:
	pprintf(prn, "\nWhite's test for heteroskedasticity\n");
	break;	
    case AUX_CHOW:
	pprintf(prn, "\nAugmented regression for Chow test\n");
	break;
    case AUX_COINT:
	pprintf(prn, "\nCointegrating regression - \n");
	break;
    case AUX_ADF:
	pprintf(prn, "\nAugmented Dickey-Fuller regression\n");
	break;
    case VAR:
	break;
    case AUX_ADD:
    default:
	if (pmod->ID < 0) pprintf(prn, "\n");
	if (pmod->name) pprintf(prn, "\n%s:\n", pmod->name);
	else pprintf(prn, "\nMODEL %d: ", pmod->ID);
	break;
    }

    if (pmod->ci == OLS || pmod->ci == VAR) pprintf(prn, "OLS ");
    else if (pmod->ci == WLS) pprintf(prn, "WLS "); 
    else if (pmod->ci == ARCH) pprintf(prn, "WLS (ARCH) ");
    else if (pmod->ci == CORC) pprintf(prn, "Cochrane-Orcutt ");
    else if (pmod->ci == HILU) pprintf(prn, "Hildreth-Lu ");
    else if (pmod->ci == TSLS) pprintf(prn, "TSLS ");
    else if (pmod->ci == HSK) pprintf(prn, "Heteroskedasticity ");
    else if (pmod->ci == AR) pprintf(prn, "AR ");
    else if (pmod->ci == HCCM) pprintf(prn, "HCCM ");
    else if (pmod->ci == PROBIT) pprintf(prn, "Probit ");
    else if (pmod->ci == LOGIT) pprintf(prn, "Logit ");
    else if (pmod->ci == POOLED) pprintf(prn, "Pooled OLS ");
    pprintf(prn, "estimates using the %d observations %s-%s\n",
	   t2-t1+1, startdate, enddate);
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG)
	pprintf(prn, "Dependent variable: uhat");
    else pprintf(prn, "Dependent variable: %s\n", 
		 pdinfo->varname[pmod->list[1]]);
    if (pmod->ci == WLS || pmod->ci == ARCH) 
	pprintf(prn, "Variable used as weight: %s\n", 
		pdinfo->varname[pmod->nwt]);
    if (pmod->infomsg[0] != '\0') pprintf(prn, "%s\n", pmod->infomsg);
    if (pmod->wt_dummy) 
	pprintf(prn, "Weight var is a dummy variable, effective obs = %d\n\n",
		pmod->nobs);
    else pprintf(prn, "\n");

    if (pmod->ci == PROBIT || pmod->ci == LOGIT) {
	print_discrete_stats(pmod, pdinfo, prn);
	return;
    }

    pprintf(prn, "      VARIABLE      COEFFICIENT      STDERROR       "
	    "T STAT    2Prob(t > |T|)\n\n");

    if (pmod->ifc) {
	print_coeff(pdinfo, pmod, ncoeff, prn);
	ncoeff--;
    }
    for (i=2; i<=ncoeff; i++) print_coeff(pdinfo, pmod, i, prn);
    pprintf(prn, "\n");

    if (pmod->aux == AUX_ARCH || pmod->aux == AUX_ADF)
	return;
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG) {
	_rsqline(pmod, prn);
	return;
    }

    if (!pmod->ifc) noconst(prn);
    
    if (pmod->aux == AUX_WHITE) {
	_rsqline(pmod, prn);
	pprintf(prn, "\nTest statistic: TR^2 = %f,\n", 
		pmod->rsq * pmod->nobs);
	pprintf(prn, "with p-value = prob(Chi-square(%d) > %f) = %f\n\n", 
		pmod->ncoeff - 1, pmod->rsq * pmod->nobs,
		chisq(pmod->rsq * pmod->nobs, pmod->ncoeff - 1)); 
	return;
    }

    if (pmod->ci == OLS || pmod->ci == VAR || pmod->ci == TSLS
	|| pmod->ci == HCCM || pmod->ci == POOLED ||
	(pmod->ci == WLS && pmod->wt_dummy)) {
	_depvarstats(pmod, prn);
	if (_essline(pmod, prn, 0)) return;
	_rsqline(pmod, prn);
	_Fline(pmod, prn);
	if (pmod->ci == OLS || (pmod->ci == WLS && pmod->wt_dummy)) {
	    if (pmod->ldepvar) dhline(pmod, prn);
	    else _dwline(pmod, prn);
	}
	/* FIXME -- check output below */
	if (pmod->ci == HCCM || pmod->ci == TSLS) _dwline(pmod, prn);
	if (pmod->ci == TSLS) pprintf(prn, "\n"
	       "R-squared is computed as the square of the correlation "
	       "between observed and\nfitted values of the dependent "
	       "variable.\n");
	print_aicetc(pmod, prn);
	_pmax_line(pmod, pdinfo, prn);
    }
    else if ((pmod->ci == WLS && !(pmod->wt_dummy)) || 
	     pmod->ci == HSK || pmod->ci == ARCH) {
	pprintf(prn, "Statistics based on the weighted data:\n\n"
	       "R-squared is suppressed as it is not meaningful.  The "
	       "F-statistic tests\nthe hypothesis that all parameters "
	       "including the constant term are zero.\n\n");
	if (_essline(pmod, prn, 1)) return;
	_Fline(pmod, prn);
	_dwline(pmod, prn);
	pprintf(prn, "\nStatistics based on the original data:\n\n"
	       "R-squared is computed as the square of the correlation "
	       "between observed and\nfitted values of the dependent "
	       "variable.\n\n");
	_depvarstats(pmod, prn);
	if (_essline(pmod, prn, 0)) return;
	_rsqline(pmod, prn); 
	print_aicetc(pmod, prn);
	_pmax_line(pmod, pdinfo, prn);
    }
    else if (pmod->ci == CORC || pmod->ci == HILU) {
	pprintf(prn, "Statistics based on the rho-differenced data:\n\n"
	       "R-squared is computed as the square of the correlation "
	       "between observed and\nfitted values of the dependent "
	       "variable.\n\n");	
	if (_essline(pmod, prn, 0)) return;
	_rsqline(pmod, prn);
	_Fline(pmod, prn);
	_dwline(pmod, prn);
	print_aicetc(pmod, prn);
	_pmax_line(pmod, pdinfo, prn);
    }
    print_model_tests(pmod, prn);
}

/* ........................................................... */

void print_add (const COMPARE *add, const int *addvars, 
		const DATAINFO *pdinfo, const int aux_code, print_t *prn)
{
    int i;
    char spc[3];

    if (aux_code != AUX_SQ && aux_code != AUX_LOG) {
	strcpy(spc, "  ");
	pprintf(prn, "Comparison of Model %d and Model %d:\n", 
		add->m1, add->m2);
    } else spc[0] = '\0';

    if (aux_code == AUX_ADD && addvars[0] > 1 && add->ols) {
	pprintf(prn, "\n%sNull hypothesis: the regression parameters are "
		"zero for the added variables\n\n", spc);
	for (i = 1; i<=addvars[0]; i++) 
	    pprintf(prn, "%s  %s\n", spc, pdinfo->varname[addvars[i]]);	
	pprintf(prn, "\n  Test statistic: F(%d, %d) = %f, ", 
		add->dfn, add->dfd, add->F);
	pprintf(prn, "with p-value = %f\n", 
		fdist(add->F, add->dfn, add->dfd));
    }
    else if (aux_code == AUX_ADD && addvars[0] > 1 && add->discrete) {
	pprintf(prn, "\n%sNull hypothesis: the regression parameters are "
		"zero for the added variables\n\n", spc);
	for (i = 1; i<=addvars[0]; i++) 
	    pprintf(prn, "%s  %s\n", spc, pdinfo->varname[addvars[i]]);	
	pprintf(prn, "\n  Test statistic: Chi-square(%d) = %f, ", 
		add->dfn, add->chisq);
	pprintf(prn, "with p-value = %f\n\n", 
		chisq(add->chisq, add->dfn));
	return;
    }
    else if (aux_code == AUX_SQ || aux_code == AUX_LOG) {
	pprintf(prn, "\nTest statistic: "
		"TR^2 = %f,\n", add->trsq);
	pprintf(prn, "with p-value = prob(Chi-square(%d) > %f) = %f\n\n", 
		add->dfn, add->trsq, chisq(add->trsq, add->dfn));
	return;
    }
    pprintf(prn, "%sOf the 8 model selection statistics, %d ", 
	    spc, add->score);
    if (add->score == 1) pprintf(prn, "has improved.\n");
    else pprintf(prn, "have improved.\n\n");
}

/* ........................................................... */

void print_omit (const COMPARE *omit, const int *omitvars, 
		 const DATAINFO *pdinfo, print_t *prn)
{
    int i;

    pprintf(prn, "Comparison of Model %d and Model %d:\n\n", 
	    omit->m1, omit->m2);
    if (omit->ols && omit->dfn > 0 && omitvars[0] > 1) {
	pprintf(prn, "  Null hypothesis: the regression parameters "
		"are zero for the variables\n\n");
	for (i = 1; i<=omitvars[0]; i++) {
	    pprintf(prn, "    %s\n", pdinfo->varname[omitvars[i]]);	
	} 
	pprintf(prn, "\n  Test statistic: F(%d, %d) = %f, ", 
		omit->dfn, omit->dfd, omit->F);
	pprintf(prn, "with p-value = %f\n", 
		fdist(omit->F, omit->dfn, omit->dfd));
    }
    else if (omit->discrete && omit->dfn > 0 && omitvars[0] > 1) {
	pprintf(prn, "  Null hypothesis: the regression parameters "
		"are zero for the variables\n\n");
	for (i = 1; i<=omitvars[0]; i++) {
	    pprintf(prn, "    %s\n", pdinfo->varname[omitvars[i]]);	
	} 
	pprintf(prn, "\n  Test statistic: Chi-square(%d) = %f, ", 
		omit->dfn, omit->chisq);
	pprintf(prn, "with p-value = %f\n\n", 
		chisq(omit->chisq, omit->dfn));
	return;
    } 
    pprintf(prn, "  Of the 8 model selection statistics, %d ", omit->score);
    if (omit->score == 1) pprintf(prn, "has improved.\n");
    else pprintf(prn, "have improved.\n\n");    
}

/* ....................................................... */

void print_aicetc (const MODEL *pmod, print_t *prn)
{
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG ||
	pmod->aux == AUX_COINT || pmod->aux == AUX_WHITE ||
	pmod->aux == AUX_AR) return;

    if (pmod->dfd == 0) {
	pprintf(prn, "\n");
	return;
    }

    pprintf(prn, "\nMODEL SELECTION STATISTICS\n\n");	
    pprintf(prn, "SGMASQ    %13g     AIC       %13g     FPE       %12g\n"
	    "HQ        %13g     SCHWARZ   %13g     SHIBATA   %12g\n"
	    "GCV       %13g",
	    pmod->criterion[0], pmod->criterion[1], 
	    pmod->criterion[2], pmod->criterion[3], 
	    pmod->criterion[4], pmod->criterion[5], pmod->criterion[6]);
    if (pmod->criterion[7] > 0.0) pprintf(prn, "     RICE      %13g\n", 
					  pmod->criterion[7]);
    else pprintf(prn, "     RICE          undefined\n");
    pprintf(prn, "\n");
}

/* ........................................................ */

static int make_list (int **plist, const DATAINFO *pdinfo)
{
    int i, n = 1;
    int *trylist;

    trylist = malloc(pdinfo->v * sizeof *trylist);
    if (trylist == NULL) return 1;
    for (i=1; i<pdinfo->v; i++) {
	if (hidden_var(i, pdinfo)) continue;
	trylist[n++] = i;
    }
    trylist[0] = n - 1;
    *plist = trylist;
    return 0;
}

/* ......................................................... */ 

void printcorr (int *list, const CORRMAT corrmat, 
		const DATAINFO *pdinfo, print_t *prn)
{
    int i = 1, j, k = 0, m, ncoeffs, freelist = 0;
    char corrstring[25];
    int *tmplist;

    if (list == NULL) {
	make_list(&tmplist, pdinfo);
	list = tmplist;
	freelist = 1;
    }

    m = corrmat.nvar;
    ncoeffs = (m * m - m)/2;

    pprintf(prn, "\nPairwise correlation coefficients:\n\n");
    while (k < ncoeffs) {
        for (i=1; i<=m; i++) {
	    for (j=i+1; j<=m; j++) {
		sprintf(corrstring, "corr(%s, %s)", 
			pdinfo->varname[list[i]], 
			pdinfo->varname[list[j]]);
		if (na(corrmat.r[k]))
		    pprintf(prn, "  %-24s is undefined\n", corrstring);
		else if (corrmat.r[k] < 0.) 
		    pprintf(prn, "  %-24s = %.3f\n", corrstring, 
			    corrmat.r[k]);
		else 
		    pprintf(prn, "  %-24s =  %.3f\n", corrstring, 
			    corrmat.r[k]);
		k++;
	    }
        }
    }
    pprintf(prn, "\n");
    if (freelist) free(list);
}

/* ......................................................... */ 

void printfreq (FREQDIST *freq, print_t *prn)
{
    int i, k, nlw, K = freq->numbins - 1;
    char word[32];

    pprintf(prn, "\nFrequency distribution for %s, obs %d-%d\n", 
	   freq->varname, freq->t1 + 1, freq->t2 + 1);
    pprintf(prn, "number of bins = %d, mean = %.3f, sd = %.3f\n", 
	   freq->numbins, freq->xbar, freq->sdx);
    pprintf(prn, "\n       interval          midpt      frequency\n\n");

    for (k=0; k<=K; k++) {
	*word = '\0';
	if (k == 0) pprintf(prn, "          <  ");
	else if (k == K) pprintf(prn, "          >= ");
	else pprintf(prn, "%10.3g - ", freq->endpt[k]);
	if (k == K) sprintf(word, "%.3g", freq->endpt[k]);
	else sprintf(word, "%.3g", freq->endpt[k+1]);
	pprintf(prn, "%s", word);
	nlw = 10 - strlen(word);
	space(nlw, prn);
	sprintf(word, " %.3g", freq->midpt[k]);
	pprintf(prn, "%s", word);
	nlw = 10 - strlen(word);
	space(nlw, prn);
	pprintf(prn, "%6d  ", freq->f[k]);
	i = 36.0 * freq->f[k]/freq->n;
	while (i--) pprintf(prn, "*");
	pprintf(prn, "\n");
    }

    pprintf(prn, "\nTest for null hypothesis of normal distribution:\n");
    pprintf(prn, "Chi-squared(2) = %.3f with pvalue %.5f\n", 
	   freq->chisqu, chisq(freq->chisqu, 2));
    if (strlen(freq->errmsg) > 2) 
	pprintf(prn, "%s", freq->errmsg);
}

/* ......................................................... */ 

void print_smpl (const DATAINFO *pdinfo, int fulln, print_t *prn)
{
    char date1[8], date2[8];

    if (fulln) {
	pprintf(prn, "Full data set: %d observations\n"
		"Current sample: %d observations\n", 
		fulln, pdinfo->n);
	return;
    }

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);
    pprintf(prn, "Full data range: %s - %s (n = %d)\n",
	    pdinfo->stobs, pdinfo->endobs, pdinfo->n);
    pprintf(prn, "Current sample:  %s - %s", date1, date2);
    if (pdinfo->t1 == 0 && pdinfo->t2 == pdinfo->n - 1) 
	pprintf(prn, "\n");
    else pprintf(prn, " (n = %d)\n", pdinfo->t2 - pdinfo->t1 + 1);  
}

/* ......................................................... */ 

static void print_float_10 (const double x, print_t *prn)
{
    int pad;
    char numstr[10];
    double xx = x;

    if (fabs(x) < 1.0e-14) xx = 0;  /* is this wise? */

    if (xx == 0.) {
	pprintf(prn, "%10.3f", xx);
	return;
    }
    if (fabs(xx) >= 1000000) {
	if (xx < 0.0) sprintf(numstr, "%.4g", xx);
	else sprintf(numstr, "%.5g", xx);
	pad = (10 - strlen(numstr));
	if (pad > 0) space(pad, prn);
	pprintf(prn, "%s", numstr);
	return;
    }
    if (fabs(xx) >= 100000) {
	pprintf(prn, "%10.0f", xx);
	return;
    }
    if (fabs(xx) < .001 && fabs(xx) >= .00001) {
	pprintf(prn, "%10.7f", xx);
	return;
    }
    if (fabs(xx) < .00001) {
	if (xx < 0.0) sprintf(numstr, "%.4g", xx);
	else sprintf(numstr, "%.5g", xx);
	pad = (10 - strlen(numstr));
	if (pad > 0) space(pad, prn);
	pprintf(prn, "%s", numstr);
	return;
    } 
    if (fabs(xx) >= 10000 && xx < 0.) {
	pprintf(prn, "%10.3f", xx);
	return;
    }
    pprintf(prn, "%10.4f", xx);
}

/* ......................................................... */ 

static void print_coeff_interval (const DATAINFO *pdinfo, const MODEL *pmod, 
				  const int c, const double t, print_t *prn)
{
    double maxerr;

    pprintf(prn, " %3d) %8s ", pmod->list[c], 
	   pdinfo->varname[pmod->list[c]]);
    space(6, prn);
    if (isnan(pmod->coeff[c-1]))
	pprintf(prn, "%10s", "undefined");
    else print_float_10 (pmod->coeff[c-1], prn);
    space(4, prn);
    if (isnan(pmod->sderr[c-1])) {
	pprintf(prn, "%10s", "undefined");
    } else {
	if (pmod->sderr[c-1] > 0)
	    maxerr = t * pmod->sderr[c-1];
	else maxerr = 0;
        printxs(pmod->coeff[c-1] - maxerr, 15, PRINT, prn);
        pprintf(prn, " to");
        printxs(pmod->coeff[c-1] + maxerr, 10, PRINT, prn);
    }
    pprintf(prn, "\n");
}

/* ......................................................... */ 

static void print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			 const int c, print_t *prn)
{
    double t, pvalue;

    pprintf(prn, " %3d) %8s ", pmod->list[c], 
	   pdinfo->varname[pmod->list[c]]);
    space(6, prn);
    if (isnan(pmod->coeff[c-1]))
	pprintf(prn, "%10s", "undefined");
    else print_float_10 (pmod->coeff[c-1], prn);
    space(4, prn);
    if (isnan(pmod->sderr[c-1])) {
	pprintf(prn, "%10s", "undefined");
	pprintf(prn, "%12s", "undefined");
	pprintf(prn, "%14s", "undefined");
	pvalue = 999.0;
    } else {
	print_float_10 (pmod->sderr[c-1], prn); 
	if (pmod->sderr[c-1] > 0.) {
	    t = pmod->coeff[c-1]/pmod->sderr[c-1];
	    if (pmod->aux == AUX_ADF) {
		pvalue = 1.;
		pprintf(prn, " %12.3f %13s", t, "unknown");
	    } else {
		pvalue = tprob(t, pmod->dfd);
		pprintf(prn, " %12.3f %14f", t, pvalue);
	    }
	} 
	else {
	    pvalue = 1.;
	    pprintf(prn, "     %12s", "undefined");
	}
    }
    if (pvalue < 0.01) pprintf(prn, " ***");
    else if (pvalue < 0.05) pprintf(prn, " **");
    else if (pvalue < 0.10) pprintf(prn, " *");
    pprintf(prn, "\n");
}

/* ......................................................... */ 

void print_rho (int *arlist, const MODEL *pmod, 
		const int c, print_t *prn)
{
    char ustr[5];
    
    sprintf(ustr, "u_%d", arlist[c]);
    pprintf(prn, "%14s ", ustr); 
    space(6, prn);
    print_float_10 (pmod->coeff[c], prn);
    space(4, prn);
    print_float_10 (pmod->sderr[c], prn); 
    pprintf(prn, " %12.3f %14f\n",
	   pmod->coeff[c]/pmod->sderr[c],
	   tprob(pmod->coeff[c]/pmod->sderr[c], 
		 pmod->dfd));	
}

/* ......................................................... */ 

int outcovmx (MODEL *pmod, const DATAINFO *pdinfo, const int batch, 
	      print_t *prn)
{
    int k, nbetas;
    int *tmplist = NULL;

    nbetas = pmod->list[0] - 1;
    if (copylist(&tmplist, pmod->list)) return E_ALLOC;
    for (k=1; k<=nbetas; k++) tmplist[k] = pmod->list[k+1];
    tmplist[0] = nbetas;

    if (pmod->vcv == NULL && makevcv(pmod)) return E_ALLOC;
    _mxout(pmod->vcv, tmplist, pmod->ci, pdinfo, batch, prn);  

    free(tmplist);
    return 0;
}

/* ......................................................... */ 

void print_white_vcv (const MODEL *pmod, print_t *prn)
{
    int i, j, index, ncoeff;

    ncoeff = pmod->list[0] - 1;
    covhdr(prn);
    index = 1;
    for (i=1; i<=ncoeff; i++) { 
	for (j=i; j<=ncoeff; j++) {
	    pprintf(prn, "\tCov(%3d, %3d) = %15g\n",
                   pmod->list[i+1], pmod->list[j+1], pmod->vcv[index]);
	    index++;
	}
    }
    pprintf(prn, "\n");
}

/* ......................................................... */ 

static void outxx (const double xx, const int ci, print_t *prn)
{
    if (ci == CORR) {
	if (na(xx)) pprintf(prn, " %13s", "undefined");
	else pprintf(prn, " %13.3f", xx);
    } else {
	if (xx > -0.001 && xx < 0.001)
	    pprintf(prn, " %13e", xx);
	else pprintf(prn, " %13g", xx);
    }
}

/* ........................................................... */
  
void takenotes (const int batch, const int runit)
{
    char s[4];

    if (batch || runit == 1) return;
    puts("\nTake notes and then press return key to continue (q to quit)");
    fflush(stdout);
    fgets(s, 3, stdin);
}

/* ........................................................... */
  
int takenotes_quit (const int batch, const int runit)
{
    char s[4];

    if (batch || runit == 1) return 0;
    puts("\nTake notes and then press return key to continue (q to quit)");
    fflush(stdout);
    s[0] = '\0';
    fgets(s, 3, stdin);
    if (s[0] == 'q') return 1;
    return 0;
}

/* ........................................................... */
                
int _pgbreak_quit (const int n, int *lineno, const int batch)
{
    int runit = 0; /* FIXME */

    if (batch || (*lineno + n <= 20)) return 0;
    if (takenotes_quit(batch, runit)) return 1;
    *lineno = 1;
    return 0;
}
/* ........................................................... */
                
void _pgbreak (const int n, int *lineno, const int batch)
{
    int runit = 0; /* FIXME */

    if (batch || (*lineno + n <= 20)) return;
    takenotes(batch, runit);
    *lineno = 1;
}

/* ......................................................... */ 

void _mxout (const double *rr, const int *list, const int ci,
	     const DATAINFO *pdinfo, const int batch, print_t *prn)
/*  Given a single dimensional array, which represents a
    symmetric matrix, prints out an upper triangular matrix
    of any size. However due to screen and printer column
    limitations the program breaks up a large upper triangular
    matrix into 5 variables at a time. For example, if there
    were 10 variables the program would first print an upper
    triangular matrix of the first 5 rows and columns, then
    it would print a rectangular matrix of the first 5 rows
    but now columns 6 - 10, and finally an upper triangular
    matrix of rows 6 -10 and columns 6 - 10
*/
{
    register int i, j;
    int lo, ljnf, nf, li2, p, k, index, ij2, lineno = 0;
    char s[16];
    enum { FIELDS = 5 };

    if (ci != CORR) covhdr(prn);

    lo = list[0];
    for (i=0; i<=lo/FIELDS; i++) {
	nf = i * FIELDS;
	li2 = lo - nf;
	p = (li2 > FIELDS) ? FIELDS : li2;
	if (p == 0) break;
	_pgbreak(3, &lineno, batch);

	/* print the varname headings */
	for(j=1; j<=p; ++j)  {
	    ljnf = list[j + nf];
	    strcpy(s, pdinfo->varname[ljnf]);
	    space(9 - strlen(s), prn);
	    pprintf(prn, "%3d) %s", ljnf, s);
	}
	pprintf(prn, "\n");
	lineno += 2;;

	/* print rectangular part, if any, of matrix */
	for (j=1; j<=nf; j++) {
	    _pgbreak(1, &lineno, batch);
	    lineno++;
	    for(k=1; k<=p; k++) {
		index = ijton(j, nf+k, lo);
		outxx(rr[index], ci, prn);
	    }
	    pprintf(prn, "   (%d\n", list[j]);
	}

	/* print upper triangular part of matrix */
	for (j=1; j<=p; ++j) {
	    _pgbreak(1, &lineno, batch);
	    lineno++;
	    ij2 = nf + j;
	    space(14 * (j - 1), prn);
	    for (k=j; k<=p; k++) {
		index = ijton(ij2, nf+k, lo);
		outxx(rr[index], ci, prn);
	    }
	    pprintf(prn, "   (%d\n", list[ij2]);
	}
	pprintf(prn, "\n");
    }
}

/* ........................................................ */
  
void space (int n, print_t *prn)
{
    if (n > 0) while (n--) pprintf(prn, " ");
}

/* ........................................................ */

static void printgx (const double xx, print_t *prn)
{
    static char word[32];
    int lw;

    sprintf(word, "%11g", xx);
    lw = strlen(word);
    pprintf(prn, "%s", word);
    space(13 - lw, prn);
} 

/* ........................................................ */

void graphyzx (const int *list, const double *zy1, const double *zy2, 
	       const double *zx, int n, const char *yname, 
	       const char *xname, const DATAINFO *pdinfo, 
	       int oflag, print_t *prn)
/*
  if n > 0 graphs zy1 against zx, otherwise
  graphs zy1[i] and zy2[i] against zx[i] for i = 1, 2, .... n
  no of rows = 40 if oflag = 1, else it is = 18 or 16
*/
{
    register int i, j;
    int ix, iy1, iy2, lx, ly, xzero, yzero, nrows, nr2, ncols, nc2,
	ls, lw, t1, t2, option = 0;
    double xmin, xmax, xrange, ymin, ymax, yrange, y1min, y1max; 
    double xx, y2min, y2max;
    char p[41][132];
    static char word[32];

    if (pdinfo != NULL) {
	t1 = pdinfo->t1;
	t2 = pdinfo->t2;
    } else {
	t1 = 0;
	t2 = (n < 0)? -n - 1 : n - 1;
    }

    if (n < 0) {
	n = -n;
	option = 1;
	minmax(t1, t2, zy1, &y1min, &y1max);
	minmax(t1, t2, zy2, &y2min, &y2max);
	ymin = (y1min < y2min)? y1min : y2min;
	ymax = (y1max > y2max)? y1max : y2max;
    }
    else minmax(t1, t2, zy1, &ymin, &ymax);
    yrange = ymax - ymin;
    xzero = yzero = 0;
    /* setting the number of columns and rows to be used */
    ncols = 60;
    if (oflag == OPT_O) nrows = 40;
    else nrows = option ? 16 : 18 ;
    nr2 = nrows/2;
    nc2 = ncols/2;
    minmax(t1, t2, zx, &xmin, &xmax);
    xrange = xmax - xmin;

    /* Initialize picture matrix */
    for (i=0; i<=nrows; ++i) {
	p[i][0] = (i%5 == 0)? '+' : '|'; 
	for (j=1; j<=ncols+1; j++) p[i][j] = ' ';
    }
    /*
      if min is < 0 and max > 0, draw line at zero value
    */
    if (xmin <0 && xmax >0) {
	xzero = 0.5 -1.0*xmin*ncols/xrange;
	for (i=0; i<=nrows; i++) p[i][xzero+1] = '|';
    }
    if (ymin <0 && ymax >0) {
	yzero = 0.5 -1.0*ymin*nrows/yrange;
	for (j=0; j<=ncols; j++) p[yzero][j+1] = '-';
    }
    /*  loop replaces blanks in PICTURE with o's that correspond to the
	scaled values of the specified variables */
    if (option) for (i=0; i<n; ++i) {
	ix = (floatneq(xrange, 0.0))? 
	    ((zx[i] - xmin)/xrange)*ncols : nc2;
	iy1 = (floatneq(yrange, 0.0))? 
	    ((zy1[i] - ymin)/yrange)*nrows : nr2;
	iy2 = (floatneq(yrange, 0.0))? 
	    ((zy2[i] - ymin)/yrange)*nrows : nr2;
	if (iy1 != iy2) {
	    p[iy1][ix+1] = 'o';
	    p[iy2][ix+1] = 'x';
	}
	else p[iy1][ix+1] = '+';
    }
    else for (i=0; i<n; ++i) {
	ix = (floatneq(xrange, 0.0))? 
	    ((zx[i] - xmin)/xrange)*ncols : nc2;
	iy1 = (floatneq(yrange, 0.0))? 
	    ((zy1[i] - ymin)/yrange)*nrows : nr2;
	p[iy1][ix+1] = 'o';
    }

    /* loop prints out the matrix PICTURE that is stored in the
       2-dimensional p matrix. */
    if (!option) pprintf(prn, "%14s\n", yname);
    else if (list) 
	pprintf(prn, "%7co stands for %s and x stands for %s (+ means they "
		"are equal)\n\n%9s, %s\n", ' ', 
		yname, pdinfo->varname[list[2]], yname, 
		pdinfo->varname[list[2]]);
    for (i=nrows; i>=0; --i) {
	if (i && i == yzero) pprintf(prn, "        0.0  ");
	else if (i == nrows || i%5 == 0) {
	    xx = ymin + ((ymax-ymin) * i/nrows);
	    printgx(xx, prn);
	}
	else space(13, prn);
	for (j=0; j<=ncols+1; ++j) pprintf(prn, "%c", p[i][j]);
	pprintf(prn, "\n");
    }
    space(13, prn);
    pprintf(prn, "|");
    for (j=0; j<=ncols; j++) if (j%10 == 0) pprintf(prn, "+");
    else pprintf(prn, "-");
    pprintf(prn, "\n");
    space(14, prn);
    sprintf(word, "%g", xmin);
    lx = strlen(word);
    lw = 13 + lx;
    pprintf(prn, "%s", word);
    sprintf(word, "%s", xname);
    ly = strlen(word);
    ls = 30 - lx - ly/2;
    space(ls, prn);
    pprintf(prn, "%s", word);
    lw = lw + ls + ly; 
    sprintf(word, "%g", xmax);

    ls = strlen(word);
    if (ls < 7) space(73 - lw, prn);
    else { 
	lw = lw + ls;
	space(79 - lw, prn);
    }
    pprintf(prn, "%s\n\n", word);
}

/* ........................................................... */

static void fit_resid_head (const MODEL *pmod, const DATAINFO *pdinfo, 
			    print_t *prn)
{
    int i;
    char label[9], date1[8], date2[8]; 

    ntodate(date1, pmod->t1, pdinfo);
    ntodate(date2, pmod->t2, pdinfo);
    pprintf(prn, "\nFull data range: %s - %s (n = %d)\n",
	    pdinfo->stobs, pdinfo->endobs, pdinfo->n);
    pprintf(prn, "Model estimation range:  %s - %s", date1, date2);
    if (pmod->t1 == 0 && pmod->t2 == pdinfo->n - 1) 
	pprintf(prn, "\n");
    else pprintf(prn, " (n = %d)\n", pmod->t2 - pmod->t1 + 1); 

    pprintf(prn, "Standard error of residuals = %f\n", pmod->sigma);
    
    if (pdinfo->pd == 1) pprintf(prn, "\n Obs ");
    else pprintf(prn, "\n\n     Obs ");
    for (i=1; i<4; i++) {
	if (i == 1) strcpy(label, pdinfo->varname[pmod->list[1]]);
	if (i == 2) strcpy(label, "fitted");
	if (i == 3) strcpy(label, "residual");
	pprintf(prn, "%13s", label);
    }
    pprintf(prn, "\n");
}

/* ........................................................... */

static void varheading (const int v1, const int v2, 
			const DATAINFO *pdinfo, const int *list,
			print_t *prn)
/*  skips to new page and prints names of variables
    from v1 to v2 */
{
    int mv;
        
    if (pdinfo->pd == 1) pprintf(prn, "\n Obs ");
    else pprintf(prn, "\n\n     Obs ");
    for (mv=v1; mv<=v2; ++mv) 
	pprintf(prn, "%13s", pdinfo->varname[list[mv]]);
    pprintf(prn, "\n");
    pprintf(prn, "\n");
}

/* ........................................................... */

void printxs (double xx, int n, int ci, print_t *prn)
{
    int ls;
    char s[32];

    printxx(xx, s, ci);
    ls = strlen(s);
    pprintf(prn, " ");
    space(n-3-ls, prn);
    pprintf(prn, "%s", s);
}

/* ........................................................ */

void _printstr (print_t *prn, const double xx, int *ls)
{
    int lwrd;
    char str[32];

    printxx(xx, str, 0);
    strcat(str,"  ");
    lwrd = strlen(str);
    if (*ls+lwrd > 78) {
	*ls = 0;
	pprintf(prn, "\n");
    }
    pprintf(prn, "%s", str);
    *ls += lwrd;
}

/* ........................................................... */

static void printz (const double *z, const DATAINFO *pdinfo, 
		    print_t *prn)
/* prints series z from current sample t1 to t2 */
{
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2, ls = 0;
    double xx;

    if (isconst(t1, t2, z)) _printstr(prn, z[t1], &ls);
    else for (t=t1; t<=t2; t++) {
	xx = z[t];
	_printstr(prn, xx, &ls);
    }
    pprintf(prn, "\n");
}

/* ........................................................... */

static void bufspace (char *buf, int n)
{
    if (n > 0) while (n--) strcat(buf, " ");
}

#ifdef OLD_PRINT 
/* ........................................................... */

static void bufprintxx (const double xx, char *str)
{
#define BIG 1.0e+9
#define SMALL 1.0e-6
    char word[32];
    double xf, xxabs;
    long xi;

    if (xx < BIG) {
	xi = xx;
	xf = xx - xi;
    }
    else xf = 0.5;

    xxabs = fabs(xx);
    if (xf == 0.0) sprintf(word, "%.0f", xx); 
    else if (xxabs < SMALL) sprintf(word, "%.4g", xx);
    else if (xxabs < 1.0) sprintf(word, "%.*g", 8, xx);
    else sprintf(word, "%.*g", 6, xx);
    strcpy(str, word);	
}	 

/* ........................................................... */

static void bufprintxs (char *buf, double xx)
{
    int ls;
    char s[32];

    bufprintxx(xx, s);
    ls = strlen(s);
    strcat(buf, " ");
    bufspace(buf, 12-ls);
    strcat(buf, s);
}
#endif /* OLD_PRINT */

/* ........................................................... */

int get_signif (double *x, int n)
     /* return either (a) the number of significant digits in
	a data series (+), or (b) the number of decimal places to
	use when printing the series (-) */
{
    static char numstr[48];
    int i, j, s, smax = 0;
    int globalsmax = 6; /* stipulated max. signif. digits */
    int lead, leadmax = 0, leadmin = 99;
    double xx;

    for (i=0; i<n; i++) {
	xx = fabs(x[i]);
	sprintf(numstr, "%.12f", xx);
	s = strlen(numstr) - 1;
	for (j=s; j>0; j--) {
	    if (numstr[j] == '0') s--;
	    else if (numstr[j] == '.') continue;
	    else break;
	}
	if (s > smax) smax = s;
	lead = 0;
	for (j=0; j<=s; j++) {
	    if (xx >= 1.0 && numstr[j] != '.') lead++;
	    else break;
	}
	if (lead > leadmax) leadmax = lead;
	if (lead < leadmin) leadmin = lead;
    }
    if (smax > globalsmax) smax = globalsmax;
    if ((leadmin < leadmax) && (leadmax < smax)) 
	smax = -1 * (smax - leadmax); /* # of decimal places */
    else if (leadmax == smax)
	smax = 0;
    else if (leadmax == 0)
	smax = -1 * (smax - 1);
    return smax;
}

/* ........................................................... */

#ifdef notdef
static int get_leading_zeros (double x)
{
    x = fabs(x);

    if (x < .000001) return 6;
    if (x < .00001) return 5;
    if (x < .0001) return 4;
    if (x < .001) return 3;
    if (x < .01) return 2;
    if (x < .1) return 1;
    return 0;
}
#endif

/* ........................................................... */

/* #define PRN_DEBUG 1 */

int bufprintnum (char *buf, double x, int signif, int width)
{
    static char numstr[24];
    int i, l;

    if (signif < 0) {
#ifdef PRN_DEBUG
	    fprintf(stderr, "got %d for signif: "
		    "printing with %%.%df\n", signif, signif);
#endif
	sprintf(numstr, "%.*f", -1 * signif, x);
    } else if (signif == 0) {
#ifdef PRN_DEBUG
	    fprintf(stderr, "got 0 for signif: "
		    "printing with %%.0f\n");
#endif
	sprintf(numstr, "%.0f", x);
    } else {
	double z = fabs(x);

	if (z < 1) l = 0;
	else if (z < 10) l = 1;
	else if (z < 100) l = 2;
	else if (z < 1000) l = 3;
	else if (z < 10000) l = 4;
	else if (z < 100000) l = 5;
	else if (z < 1000000) l = 6;
	else l = 7;
	if (l >= signif) { 
#ifdef PRN_DEBUG
	    fprintf(stderr, "got %d for leftvals, %d for signif: "
		    "printing with %%.%dG\n", l, signif, signif);
#endif
	    sprintf(numstr, "%.*G", signif, x);
	} else if (z >= .10) {
	    /* was if (!get_leading_zeros(x)) */
#ifdef PRN_DEBUG
	    fprintf(stderr, "got %d for leftvals, %d for signif: "
		    "printing with %%.%df\n", l, signif, signif-l);
#endif
	    sprintf(numstr, "%.*f", signif - l, x);
	} else {
#ifdef PRN_DEBUG
	    fprintf(stderr, "got %d for leftvals, %d for signif: "
		    "printing with %%.%dG\n", l, signif, signif);
#endif
	    sprintf(numstr, "%.*G", signif, x);
	}
    }

    l = width - strlen(numstr);
    for (i=0; i<l; i++)
	strcat(buf, " ");
    strcat(buf, numstr);

    return 0;
}

/* ........................................................... */

int printdata (int *list, double **pZ, const DATAINFO *pdinfo, 
	       int batch, int byobs, print_t *prn)
/*  prints the data for the variables in list, from obs t1 to
    obs t2.
*/
{
    int idate, l0, j, v, v1, v2, j5, nvj5, lineno, ncol;
    register int t;
    int gui, isconst; 
    int *pmax = NULL; 
    int pd = pdinfo->pd, t1 = pdinfo->t1, t2 = pdinfo->t2,
	n = pdinfo->n;
    double xx, xdate, sd0 = pdinfo->sd0;
    int *tmplist = NULL, freelist = 0;
    char line[96];

    if (prn->buf != NULL) gui = 1;
    else gui = 0;

    lineno = 1;
    if (list == NULL) {
	if (make_list(&tmplist, pdinfo)) return 1;
	list = tmplist;
	freelist = 1;
    }
    l0 = list[0];

    /* special case: all vars have constant value over sample */
    isconst = 1;
    for (j=1; j<=list[0]; j++) {
	for (t=t1+1; t<=t2; t++) {
	    if (floatneq((*pZ)[n*list[j] + t], (*pZ)[n*list[j] + t1])) {
		isconst = 0;
		break;
	    }
	}
	if (!isconst) break;
    }
    if (isconst) {
	for (j=1; j<=list[0]; j++) 
	    pprintf(prn, "%8s = %10g\n", pdinfo->varname[list[j]], 
		    (*pZ)[n*list[j] + t1]);
	if (freelist) free(list);
	return 0;
    }

    if (!byobs) {
	/* print data by variables */
	for (j=1; j<=list[0]; j++) {
	    pprintf(prn, "Varname: %s\n", pdinfo->varname[list[j]]);
	    pprintf(prn, "Data frequency: %d\n", pdinfo->pd);
	    print_smpl (pdinfo, 0, prn);
	    pprintf(prn, "\n");
	    printz(&((*pZ)[n*list[j]]), pdinfo, prn);
	    pprintf(prn, "\n");
	}
	return 0;
    }

    /* experimental */
    pmax = malloc(l0 * sizeof *pmax);
    if (pmax == NULL) return 1;
    for (j=1; j<=l0; j++) {
	/* this runs fairly quickly, even for large dataset */
	pmax[j-1] = get_signif(&(*pZ)[n*list[j] + t1], t2-t1+1);
    }

    /* print data by observations */
    ncol = 5;
    for (j=0; j<=l0/ncol; j++) {
	j5 = j * ncol;
	nvj5 = l0 - j5;
	v1 = j5 +1;
	if (nvj5) {
	    v2 = (ncol > nvj5)? nvj5 : ncol;
	    v2 += j5;
	    varheading(v1, v2, pdinfo, list, prn);
	    if (!gui && _pgbreak_quit(1, &lineno, batch)) return 0;
	    lineno++;
	    for (t=t1; t<=t2; t++)   {
		if (pdinfo->markers) { /* data marker strings present */
		    sprintf(line, "%8s ", pdinfo->S[t]);
		} else {
		    xdate = date(t, pd, sd0);
		    idate = (int) xdate;
		    if (pd == 1) 
			sprintf(line, "%4d ", idate);
		    else if (pd < 10) 
			sprintf(line, "%8.1f ", xdate);
		    else 
			sprintf(line, "%8.2f ", xdate);
		} /* end print obs marker */
		for (v=v1; v<=v2; v++) {
		    xx = (*pZ)[n*list[v] + t];
		    if (na(xx)) {
			bufspace(line, 13);
		    } else { 
			/*  bufprintxs(line, xx); */ 
			bufprintnum(line, xx, pmax[v-1], 13);
		    }
		}
		if (pprintf(prn, "%s\n", line))
		    return 1;
		if (!gui && _pgbreak_quit(1, &lineno, batch)) return 0;
		lineno++;
		if (!batch) {
		    if ((t-t1+1) % 21 == 0) {
			varheading(v1, v2, pdinfo, list, prn);
			if (!gui && _pgbreak_quit(1, &lineno, batch)) return 0;
			lineno++;
		    }
		}
	    } /* end of t loop */
	} /* end if nvj5 */
    } /* end for j loop */
    pprintf(prn, "\n");
    lineno++;
    if (freelist) free(list);
    free(pmax);
    return 0;
}

/* ........................................................... */

int print_fit_resid (const MODEL *pmod, double **pZ, 
		     DATAINFO *pdinfo, print_t *prn)
{
    int pmax, idate, depvar, t, nfit, ast = 0;
    int pd = pdinfo->pd, t1 = pmod->t1, t2 = pmod->t2, n = pdinfo->n;
    double xx, xdate, sd0 = pdinfo->sd0;
    char fcastline[32];

    depvar = pmod->list[1];

    sprintf(fcastline, "fcast %s %s fitted", pdinfo->stobs, 
	    pdinfo->endobs);
    nfit = fcast(fcastline, pmod, pdinfo, pZ, NULL);
    if (nfit < 0) return 1;

    if (isdummy(depvar, t1, t2, *pZ, n) > 0)
	pmax = get_precision(&(*pZ)[n*nfit], n);
    else
	pmax = get_precision(&(*pZ)[n*depvar], n);

    fit_resid_head(pmod, pdinfo, prn);
    for (t=0; t<n; t++) {
	if (t == t1 && t) pprintf(prn, "\n");
	if (t == t2 + 1) pprintf(prn, "\n");
	if (pdinfo->markers) { /* data marker strings present */
	    pprintf(prn, "%8s ", pdinfo->S[t]); 
	} else {
	    xdate = date(t, pd, sd0);
	    idate = (int) xdate;
	    if (pd == 1) pprintf(prn, "%4d ", idate);
	    else if (pd < 10) pprintf(prn, "%8.1f ", xdate);
	    else pprintf(prn, "%8.2f ", xdate);
	}
/*  	for (i=1; i<4; i++) { */
/*  	    if (i == 1) xx = (*pZ)[n*depvar + t]; */
/*  	    if (i == 2) xx = (*pZ)[n*nfit + t]; */
/*  	    if (i == 3) xx = (*pZ)[n*depvar + t] - (*pZ)[n*nfit + t]; */
/*  	    printxs(xx, 15, PRINT, prn); */
/*  	} */
	xx = (*pZ)[n*depvar + t] - (*pZ)[n*nfit + t];
	if (fabs(xx) > 2.5 * pmod->sigma) ast = 1;
	pprintf(prn, "%12.*f%12.*f%12.*f", 
		pmax, (*pZ)[n*depvar + t],
		pmax, (*pZ)[n*nfit + t], pmax, xx);
	pprintf(prn, "%s\n", (fabs(xx) > 2.5 * pmod->sigma)? " *" : "");
    }
    pprintf(prn, "\n");
    if (ast) pprintf(prn, "Note: * denotes a residual in excess of "
		     "2.5 standard errors\n");
    return 0;
}

/* ........................................................... */

void print_ar (MODEL *pmod, print_t *prn)
{
    pprintf(prn, "Statistics based on the rho-differenced data\n"
           "(R-squared is computed as the square of the correlation "
           "between observed and\nfitted values of the dependent "
           "variable):\n\n");
    if (_essline(pmod, prn, 0)) return;
    _rsqline(pmod, prn);
    _Fline(pmod, prn);
    _dwline(pmod, prn);
    print_aicetc(pmod, prn); 
}

/* ........................................................... */

static void print_discrete_coeff (const DATAINFO *pdinfo, 
				  const MODEL *pmod, 
				  const int c, print_t *prn)
{
    double tstat;

    pprintf(prn, " %3d) %8s ", pmod->list[c], 
	   pdinfo->varname[pmod->list[c]]);
    space(6, prn);
    if (isnan(pmod->coeff[c-1]))
	pprintf(prn, "%10s", "undefined");
    else print_float_10 (pmod->coeff[c-1], prn);
    space(4, prn);
    print_float_10 (pmod->sderr[c-1], prn);
    tstat = pmod->coeff[c-1]/pmod->sderr[c-1];
    pprintf(prn, " %12.3f  ", tstat);
    if (pmod->list[c] != 0)
	print_float_10 (pmod->slope[c-1], prn); 
    pprintf(prn, "\n");
}

/* ........................................................... */

static void print_discrete_stats (const MODEL *pmod, 
				  const DATAINFO *pdinfo, 
				  print_t *prn)
{
    int i, ncoeff = pmod->list[0];

    pprintf(prn, "      VARIABLE      COEFFICIENT      STDERROR       "
	    "T STAT       SLOPE\n");
    pprintf(prn, "                                                    "
	    "           (at mean)\n");

    if (pmod->ifc) {
	print_discrete_coeff(pdinfo, pmod, ncoeff, prn);
	ncoeff--;
    }
    for (i=2; i<=ncoeff; i++) 
	print_discrete_coeff(pdinfo, pmod, i, prn);
    pprintf(prn, "\n");
    pprintf(prn, "Mean of %s = %.3f\n", 
	    pdinfo->varname[pmod->list[1]], pmod->ybar);
    pprintf(prn, "Number of cases 'correctly predicted' = %d (%.1f%%)\n", 
	    pmod->correct, 100 * (double) pmod->correct / (double) pmod->nobs);
    pprintf(prn, "f(beta'x) at mean of independent vars = %.3f\n", pmod->sdy);
    pprintf(prn, "Log-likelihood = %.3f\n", pmod->lnL);
    if (pmod->aux != AUX_OMIT && pmod->aux != AUX_ADD) {
	i = pmod->ncoeff - 1;
	pprintf(prn, "Likelihood ratio test: "
		"Chi-square(%d) = %.3f (p-value %f)\n\n",
		i, pmod->chisq, chisq(pmod->chisq, i));
    } else pprintf(prn, "\n");
}
