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

/*  modelprint.c */ 

#include "libgretl.h"
#include "internal.h"

#ifdef OS_WIN32
#define isnan(x) ((x) != (x))
#endif

#define PLAIN_FORMAT(c) (c == GRETL_PRINT_FORMAT_PLAIN)
#define RTF_FORMAT(c) (c == GRETL_PRINT_FORMAT_RTF)
#define TEX_FORMAT(c) (c == GRETL_PRINT_FORMAT_TEX || c == GRETL_PRINT_FORMAT_TEX_DOC)
#define STANDALONE(c) (c == GRETL_PRINT_FORMAT_TEX_DOC)


static int print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			int c, PRN *prn);
static int rtf_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			    int c, PRN *prn);
static void depvarstats (const MODEL *pmod, PRN *prn);
static int essline (const MODEL *pmod, PRN *prn, int wt);
static void rsqline (const MODEL *pmod, PRN *prn);
static void Fline (const MODEL *pmod, PRN *prn);
static void dwline (const MODEL *pmod, PRN *prn);
static void print_rho_terms (const MODEL *pmod, PRN *prn);
static void print_discrete_statistics (const MODEL *pmod, 
				       const DATAINFO *pdinfo,
				       PRN *prn);
static void print_aicetc (const MODEL *pmod, PRN *prn);
static void tex_print_aicetc (const MODEL *pmod, PRN *prn);
static void rtf_print_aicetc (const MODEL *pmod, PRN *prn);

/* ......................................................... */ 

static void noconst (PRN *prn)
{
    pprintf(prn, _("The model has no constant term.\n"  
	    "F is calculated as in Sect. 4.4 of Ramanathan's Introductory "
	    "Econometrics.\n"
	    "R-squared is the square of the correlation between the "
	    "observed and fitted\n values of the dependent variable.\n\n"));
}

#define RTFTAB "\\par \\ql \\tab "

/* ......................................................... */ 

static void depvarstats (const MODEL *pmod, PRN *prn)
{
    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "  %s = %g\n", _("Mean of dependent variable"), 
		pmod->ybar);
	pprintf(prn, "  %s = %g\n", _("Standard deviation of dep. var."), 
		pmod->sdy);
    }
    
    else if (TEX_FORMAT(prn->format)) {
	char x1str[32], x2str[32];

	tex_dcolumn_double(pmod->ybar, x1str);
	tex_dcolumn_double(pmod->sdy, x2str);
	pprintf(prn, "%s & %s \\\\\n %s & %s \\\\\n",
		I_("Mean of dependent variable"), x1str,
		I_("S.D. of dependent variable"), x2str);
    }

    else if (RTF_FORMAT(prn->format)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_("Mean of dependent variable"), 
		pmod->ybar);
	pprintf(prn, RTFTAB "%s = %g\n", I_("Standard deviation of dep. var."), 
		pmod->sdy);
    }
}

/* ......................................................... */ 

static int essline (const MODEL *pmod, PRN *prn, int wt)
{
    if (PLAIN_FORMAT(prn->format)) {    
	if ((wt && pmod->ess_wt < 0) || (!wt && pmod->ess < 0)) {
	    char tmp[128];

	    sprintf(tmp, _("Error sum of squares (%g) is not > 0"), 
		    (wt)? pmod->ess_wt : pmod->ess), 
		pprintf(prn, "%s\n\n", tmp);
	    return 1;
	}

	pprintf(prn, "  %s = %g\n", _("Sum of squared residuals"), 
		wt? pmod->ess_wt : pmod->ess);
	pprintf(prn, "  %s = %g\n", _("Standard error of residuals"), 
		wt? pmod->sigma_wt : pmod->sigma);
	return 0;
    }

    if (RTF_FORMAT(prn->format)) {
	if ((wt && pmod->ess_wt < 0) || (!wt && pmod->ess < 0)) {
	    char tmp[128];

	    sprintf(tmp, I_("Error sum of squares (%g) is not > 0"), 
		    (wt)? pmod->ess_wt : pmod->ess), 
		pprintf(prn, "\\par \\ql %s\\par\n\n", tmp);
	    return 1;
	}

	pprintf(prn, RTFTAB "%s = %g\n", I_("Sum of squared residuals"), 
		wt? pmod->ess_wt : pmod->ess);
	pprintf(prn, RTFTAB "%s = %g\n", I_("Standard error of residuals"), 
		wt? pmod->sigma_wt : pmod->sigma);
	return 0;
    }

    if (TEX_FORMAT(prn->format)) {
	char x1str[32], x2str[32];

	tex_dcolumn_double(pmod->ess, x1str);
	tex_dcolumn_double(pmod->sigma, x2str);
	pprintf(prn, "%s & %s \\\\\n %s ($\\hat{\\sigma}$) & %s \\\\\n",
		I_("Sum of squared residuals"), x1str,
		I_("Standard error of residuals"), x2str);
	return 0;
    }

    return 0;
}

/* ......................................................... */ 

static void rsqline (const MODEL *pmod, PRN *prn)
{
    if (PLAIN_FORMAT(prn->format)) {  
	pprintf(prn, "  %s = %g\n", _("Unadjusted R-squared"), pmod->rsq);
	if (!na(pmod->adjrsq)) {
	    pprintf(prn, "  %s = %g\n", _("Adjusted R-squared"),  
		    pmod->adjrsq);
	}
    }

    else if (RTF_FORMAT(prn->format)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_("Unadjusted R{\\super 2}"), pmod->rsq);
	if (!na(pmod->adjrsq)) {
	    pprintf(prn, RTFTAB "%s = %g\n", I_("Adjusted R{\\super 2}"),  
		    pmod->adjrsq);
	}	
    }

    else if (TEX_FORMAT(prn->format)) {  
	char x1str[32], x2str[32];

	tex_dcolumn_double(pmod->rsq, x1str);
	tex_dcolumn_double(pmod->adjrsq, x2str);
	pprintf(prn, "%s & %s \\\\\n %s & %s \\\\\n",
		I_("Unadjusted $R^2$"), x1str, 
		I_("Adjusted $\\bar{R}^2$"), x2str);
    }
}

/* ......................................................... */ 

static void ladstats (const MODEL *pmod, PRN *prn)
{
    int utf = PLAIN_FORMAT(prn->format);

    if (TEX_FORMAT(prn->format)) {  
	char x1str[32];

	tex_dcolumn_double(pmod->ess, x1str);
	pprintf(prn, "%s & %s \\\\\n",
		I_("Sum of absolute residuals"), x1str); 
    }
    
    else {
	pprintf(prn, "  %s = %g\n", 
		(utf)? _("Sum of absolute residuals") :
		I_("Sum of absolute residuals"),
		pmod->ess);
    }
}

/* ......................................................... */ 

static void print_f_pval_str (double pval, PRN *prn)
{
    int utf = PLAIN_FORMAT(prn->format);

    if (utf || RTF_FORMAT(prn->format)) {
	if (pval < .00001) {
	    pprintf(prn, " (%s < %.5f)\n", 
		    (utf)? _("p-value") : I_("p-value"), 0.00001);
	} else {
	    pprintf(prn, " (%s = %.3g)\n", 
		    (utf)? _("p-value") : I_("p-value"), pval);
	}
    }

    else if (TEX_FORMAT(prn->format)) {
	if (pval < .00001) {
	    return;
	} else {
	    char pstr[32];

	    tex_dcolumn_double(pval, pstr);
	    pprintf(prn, "%s $F()$ & %s \\\\\n", I_("p-value for"), pstr);
	}
    }
}

/* ......................................................... */ 

static void Fline (const MODEL *pmod, PRN *prn)
{
    if (PLAIN_FORMAT(prn->format)) {
	char tmp[32];

	sprintf(tmp, "%s (%d, %d)", _("F-statistic"), pmod->dfn, pmod->dfd);
	if (na(pmod->fstt)) {
	    pprintf(prn, "  %s = %s", tmp, _("undefined"));
	} else {
	    pprintf(prn, "  %s = %g", tmp, pmod->fstt);
	    print_f_pval_str(fdist(pmod->fstt, pmod->dfn, pmod->dfd), prn);
	}
    }

    else if (TEX_FORMAT(prn->format)) {
	char x1str[32];

	tex_dcolumn_double(pmod->fstt, x1str);
	pprintf(prn, "$F(%d, %d)$ & %s \\\\\n", pmod->dfn, pmod->dfd, x1str);
	print_f_pval_str(fdist(pmod->fstt, pmod->dfn, pmod->dfd), prn);
    }

    else if (RTF_FORMAT(prn->format)) {
	char tmp[32];

	sprintf(tmp, "%s (%d, %d)", I_("F-statistic"), pmod->dfn, pmod->dfd);
	if (na(pmod->fstt)) {
	    pprintf(prn, RTFTAB "%s = %s", tmp, I_("undefined"));
	} else {
	    pprintf(prn, RTFTAB "%s = %g", tmp, pmod->fstt);
	    print_f_pval_str(fdist(pmod->fstt, pmod->dfn, pmod->dfd), prn);
	}
    }
}

/* ......................................................... */ 

static void dwline (const MODEL *pmod, PRN *prn)
{
    if (PLAIN_FORMAT(prn->format)) {
	if (!na(pmod->dw)) {
	    pprintf(prn, "  %s = %g\n", _("Durbin-Watson statistic"), 
		    pmod->dw);
	    pprintf(prn, "  %s = %g\n", _("First-order autocorrelation coeff."), 
		    pmod->rho);
	} 
    }

    else if (TEX_FORMAT(prn->format)) {
	char x1str[32], x2str[32];
	tex_dcolumn_double(pmod->dw, x1str);
	tex_dcolumn_double(pmod->rho, x2str);
	pprintf(prn, "%s & %s \\\\\n %s ($\\hat{\\rho}$) & %s \n",
		I_("Durbin--Watson statistic"), x1str, 
		I_("First-order autocorrelation coeff."), x2str);
    }

    else if (RTF_FORMAT(prn->format)) {
	if (!na(pmod->dw)) {
	    pprintf(prn, RTFTAB "%s = %g\n", I_("Durbin-Watson statistic"), 
		    pmod->dw);
	    pprintf(prn, RTFTAB "%s = %g\n", I_("First-order autocorrelation coeff."), 
		    pmod->rho);
	} 
    }
}

/* ......................................................... */ 

static void dhline (const MODEL *pmod, PRN *prn)
{
    double sderr, h = 0.0;
    int i = pmod->ldepvar, T = pmod->nobs - 1;

    sderr = pmod->sderr[i-1];

    if ((T * sderr * sderr) >= 1.0) return;

    h = pmod->rho * sqrt(T/(1 - T * sderr * sderr));

    if (PLAIN_FORMAT(prn->format)) {
	char tmp[128];

	sprintf(tmp, _("Durbin's h stat. %g  First-order autocorr. coeff %g"), 
		h, pmod->rho);
	pprintf(prn, "  %s\n", tmp);

	sprintf(tmp, _("(Using variable %d for h stat, with T' = %d)"), 
		pmod->list[i], T);
	pprintf(prn, "  %s\n", tmp);
    }

    else if (RTF_FORMAT(prn->format)) {
	char tmp[128];

	h = pmod->rho * sqrt(T/(1 - T * sderr * sderr));
	sprintf(tmp, I_("Durbin's h stat. %g  First-order autocorr. coeff %g"), 
		h, pmod->rho);
	pprintf(prn, RTFTAB "%s\n", tmp);

	sprintf(tmp, I_("(Using variable %d for h stat, with T' = %d)"), 
		pmod->list[i], T);
	pprintf(prn, RTFTAB "%s\n", tmp);
    }

    else if (TEX_FORMAT(prn->format)) {
	char x1str[32];

	tex_dcolumn_double(h, x1str);
	pprintf(prn, "%s & %s \\\\\n",
		I_("Durbin's $h$ statistic"), x1str);
    }
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
    if (tprob(tmin, pmod->dfd) > .10) {
	return pmod->list[k+1];
    }
    return 0;
}

/* ......................................................... */ 

static void pmax_line (const MODEL *pmod, const DATAINFO *pdinfo, 
		       PRN *prn)
{
    if (!TEX_FORMAT(prn->format)) {
	int k = pmod->ncoeff - pmod->ifc;

	if (k < 2) return;

	if ((k = _pmax(pmod))) {
	    char tmp[128];

	    if (PLAIN_FORMAT(prn->format)) {
		sprintf(tmp, _("Excluding the constant, p-value was highest "
			"for variable %d (%s)"), k, pdinfo->varname[k]);
		pprintf(prn, "%s\n\n", tmp);
	    }
	    else if (RTF_FORMAT(prn->format)) {
		sprintf(tmp, I_("Excluding the constant, p-value was highest "
			"for variable %d (%s)"), k, pdinfo->varname[k]);
		pprintf(prn, "\\par %s\\par\n", tmp);
	    }
	}
    }
}

/* ......................................................... */

static const char *aux_string (int aux, int format)
{
    if (aux == AUX_SQ)
	return N_("Auxiliary regression for non-linearity test "
		 "(squared terms)");
    else if (aux == AUX_LOG)
	return N_("Auxiliary regression for non-linearity test "
		 "(log terms)");
    else if (aux == AUX_WHITE)
	return N_("White's test for heteroskedasticity");
    else if (aux == AUX_CHOW)
	return N_("Augmented regression for Chow test");
    else if (aux == AUX_COINT) {
	if (TEX_FORMAT(format)) return N_("Cointegrating regression -- ");
	else return N_("Cointegrating regression - ");
    }
    else if (aux == AUX_ADF) {
	if (TEX_FORMAT(format)) return N_("Augmented Dickey--Fuller regression");
	return N_("Augmented Dickey-Fuller regression");
    }
    else return "";
}

/* ......................................................... */

static const char *estimator_string (int ci, int format)
{
    if (ci == OLS || ci == VAR) return N_("OLS");
    else if (ci == WLS) return N_("WLS"); 
    else if (ci == ARCH) return N_("WLS (ARCH)");
    else if (ci == TSLS) return N_("TSLS");
    else if (ci == HSK) return N_("Heteroskedasticity");
    else if (ci == AR) return N_("AR");
    else if (ci == LAD) return N_("LAD");
    else if (ci == HCCM) return N_("HCCM");
    else if (ci == PROBIT) return N_("Probit");
    else if (ci == LOGIT) return N_("Logit");
    else if (ci == POOLED) return N_("Pooled OLS");
    else if (ci == CORC) {
	if (TEX_FORMAT(format)) return N_("Cochrane--Orcutt");
	else return N_("Cochrane-Orcutt");
    }
    else if (ci == HILU) {
	if (TEX_FORMAT(format)) return N_("Hildreth--Lu");
	else return N_("Hildreth-Lu");
    }

    else return "";
}

/* ......................................................... */

void get_test_stat_string (GRETLTEST *test, char *str, int format)
{
    int tex = TEX_FORMAT(format);

    switch (test->teststat) {
    case GRETL_TEST_TR2:
	if (tex) {
	    sprintf(str, "$TR^2$ = %g", test->value);
	} else if (RTF_FORMAT(format)) {
	    sprintf(str, "TR{\\super 2} = %g", test->value);
	} else {
	    sprintf(str, "TR^2 = %g", test->value);
	}
	break;
    case GRETL_TEST_F:
	if (tex) 
	    sprintf(str, "$F(%d, %d)$ = %g", test->dfn, test->dfd, test->value);
	else 
	    sprintf(str, "F(%d, %d) = %g", test->dfn, test->dfd, test->value);
	break;
    case GRETL_TEST_LMF:
	sprintf(str, "LMF = %g", test->value);
	break;
    case GRETL_TEST_HARVEY_COLLIER:
	if (tex)
	    sprintf(str, "Harvey--Collier $t(%d)$ = %g", test->dfn, test->value);
	else
	    sprintf(str, "Harvey-Collier t(%d) = %g", test->dfn, test->value);
	break;
    default:
	*str = 0;
    }
}

void get_test_pval_string (GRETLTEST *test, char *str, int format)
{
    int tex = TEX_FORMAT(format);

    switch (test->teststat) {
    case GRETL_TEST_TR2:
	if (tex) sprintf(str, "$P$($\\chi^2_{%d} >$ %g) = %g", 
			 test->dfn, test->value, test->pvalue);
	else sprintf(str, "P(Chi-Square(%d) > %g) = %g", 
		     test->dfn, test->value, test->pvalue);
	break;
    case GRETL_TEST_F:
	if (tex) 
	    sprintf(str, "$P$($F(%d, %d) >$ %g)$ = %g", 
		    test->dfn, test->dfd, test->value, test->pvalue);
	else 
	    sprintf(str, "P(F(%d, %d) > %g) = %g", 
		    test->dfn, test->dfd, test->value, test->pvalue);
	break;
    case GRETL_TEST_LMF:
	if (tex) 
	    sprintf(str, "$P$($\\chi^2_{%d} >$ > %g) = %g", 
		    test->dfn, test->value, test->pvalue);
	else
	    sprintf(str, "P(Chi-Square(%d) > %g) = %g", 
		    test->dfn, test->value, test->pvalue);
	break;
    case GRETL_TEST_HARVEY_COLLIER:
	if (tex)
	    sprintf(str, "$P$($t_{%d} >$ %g)  = %g", 
		    test->dfn, test->value, test->pvalue);
	else
	    sprintf(str, "P(t(%d) > %g) = %g", 
		    test->dfn, test->value, test->pvalue);
	break;
    default:
	*str = 0;
    }
}

/* ......................................................... */

static void print_model_tests (const MODEL *pmod, PRN *prn)
{
    int i;
    char test_str[64], pval_str[64];

    if (PLAIN_FORMAT(prn->format)) {
	for (i=0; i<pmod->ntests; i++) {
	    get_test_stat_string(&pmod->tests[i], test_str, prn->format);
	    get_test_pval_string(&pmod->tests[i], pval_str, prn->format);
	    pprintf(prn, "%s -\n"
		    "  %s: %s\n"
		    "  %s: %s\n"
		    "  %s = %s\n\n",
		    (pmod->tests[i]).type, 
		    _("Null hypothesis"), (pmod->tests[i]).h_0, 
		    _("Test statistic"), test_str, 
		    _("with p-value"), pval_str);
	}
    }

    else if (TEX_FORMAT(prn->format)) {
	if (pmod->ntests > 0) {
	    pprintf(prn, "\\vspace{1em}\n\\begin{raggedright}\n");
	    for (i=0; i<pmod->ntests; i++) {
		if (i > 0) {
		    pprintf(prn, "\\vspace{1ex}\n");;
		}
		get_test_stat_string(&pmod->tests[i], test_str, prn->format);
		get_test_pval_string(&pmod->tests[i], pval_str, prn->format);
		pprintf(prn, "%s --\\\\\n"
			"\\quad %s: %s\\\\\n"
			"\\quad %s: %s\\\\\n"
			"\\quad %s = %s\\\\\n",
			(pmod->tests[i]).type, 
			I_("Null hypothesis"), (pmod->tests[i]).h_0, 
			I_("Test statistic"), test_str, 
			I_("with p-value"), pval_str);
	    }
	    pprintf(prn, "\\end{raggedright}\n");
	}
    }

    else if (RTF_FORMAT(prn->format)) {
	for (i=0; i<pmod->ntests; i++) {
	    get_test_stat_string(&pmod->tests[i], test_str, prn->format);
	    get_test_pval_string(&pmod->tests[i], pval_str, prn->format);
	    pprintf(prn, "\\par \\ql %s -\\par\n"
		    " %s: %s\\par\n"
		    " %s: %s\\par\n"
		    " %s = %s\\par\n\n",
		    (pmod->tests[i]).type, 
		    I_("Null hypothesis"), (pmod->tests[i]).h_0,
		    I_("Test statistic"), test_str, 
		    I_("with p-value"), pval_str);
	}
    }
}

/* ......................................................... */

static void modelprint_setup_obs (const MODEL *pmod, int *t1, int *t2)
{
    if (pmod->ci == CORC || pmod->ci == HILU) {
	*t1 += 1;
    }

    if (pmod->data != NULL) {
	*t2 += get_misscount(pmod);
    }
}

/* ......................................................... */

static void print_model_heading (const MODEL *pmod, 
				 const DATAINFO *pdinfo, 
				 PRN *prn)
{
    char startdate[9], enddate[9], vname[16];
    int t1 = pmod->t1, t2 = pmod->t2;
    int tex = TEX_FORMAT(prn->format);
    int utf = PLAIN_FORMAT(prn->format);

    modelprint_setup_obs(pmod, &t1, &t2);

    ntodate(startdate, t1, pdinfo);
    ntodate(enddate, t2, pdinfo);

    switch (pmod->aux) {
    case AUX_SQ:
    case AUX_LOG:
    case AUX_WHITE:
    case AUX_CHOW:
    case AUX_COINT:
    case AUX_ADF:
	if (utf) {
	    pprintf(prn, "\n%s\n", _(aux_string(pmod->aux, prn->format)));
	} else if (tex) {
	    pprintf(prn, "\n%s\n", I_(aux_string(pmod->aux, prn->format)));
	} else { /* RTF */
	    pprintf(prn, "%s\\par\n", I_(aux_string(pmod->aux, prn->format)));
	}
	break;
    case AUX_AR:
	if (utf) { 	
	    pprintf(prn, "\n%s ", _("Breusch-Pagan test for"));
	} else {
	    pprintf(prn, "\n%s ", I_("Breusch-Pagan test for"));
	} 
	if (pmod->order > 1) {
	    pprintf(prn, "%s %d\n", (utf)? _("autocorrelation up to order") :
		    I_("autocorrelation up to order"), 
		    pmod->order);
	} else {
	    pprintf(prn, "%s\n", (utf)? _("first-order autocorrelation") :
		    I_("first-order autocorrelation"));
	}
	break;	
    case AUX_ARCH:
	pprintf(prn, "\n%s %d\n", 
		(utf)? _("Test for ARCH of order") : 
		I_("Test for ARCH of order"), 
		pmod->order);
	break;	
    case VAR:
	break;
    case AUX_ADD:
    default:
	if (pmod->ID < 0) pprintf(prn, "\n");
	else if (pmod->name) {
	    pprintf(prn, "\n%s:\n", pmod->name);
	} else {
	    pprintf(prn, "\n%s %d: ", (utf)? _("Model") : I_("Model"), pmod->ID);
	}
	break;
    }

    pprintf(prn, (utf)?
	    _("%s estimates using the %d observations %s%s%s") :
	    I_("%s estimates using the %d observations %s%s%s"),
	    _(estimator_string(pmod->ci, prn->format)), 
	    pmod->nobs, startdate, (tex)? "--" : "-", enddate);

    if (tex) pprintf(prn, "\\\\\n");
    else if (RTF_FORMAT(prn->format)) pprintf(prn, "\\par\n");
    else pprintf(prn, "\n");

    /* special names for dependent variable in cases of certain sorts
       of auxiliary regressions */
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG) {
	pprintf(prn, "%s: %s\n", _("Dependent variable"),
		(tex)? "$\\hat{u}$" : "uhat");
    }
    else if (pmod->aux == AUX_WHITE) {
	pprintf(prn, "%s: uhat^2\n", 
		(utf)? _("Dependent variable") : I_("Dependent variable"),
		(tex)? "$\\hat{u}^2$" : "uhat^2");
    }
    else if (pmod->aux == AUX_ARCH) {
	pprintf(prn, "%s: ut^2\n", 
		(utf)? _("Dependent variable") : I_("Dependent variable"),
		(tex)? "$u_t^2$" : "ut^2");
    }
    else { /* ordinary dependent variable */
	if (tex) tex_escape(vname, pdinfo->varname[pmod->list[1]]);
	pprintf(prn, "%s: %s\n", 
		(utf)? _("Dependent variable") : I_("Dependent variable"),
		(tex)? vname : pdinfo->varname[pmod->list[1]]);
    }

    if (pmod->ci == WLS || pmod->ci == ARCH) {
	if (tex) {
	    tex_escape(vname, pdinfo->varname[pmod->nwt]);
	    pprintf(prn, "\\\\\n");
	}
	pprintf(prn, "%s: %s\n", 
		(utf)? _("Variable used as weight") : I_("Variable used as weight"), 
		(tex)? vname : pdinfo->varname[pmod->nwt]);
    }

    if (PLAIN_FORMAT(prn->format) && gretl_msg[0] != '\0') {
	pprintf(prn, "%s\n", gretl_msg);
    }

    if (pmod->wt_dummy) { /* FIXME alt formats */
	pprintf(prn, "%s %d\n", 
		(utf)? _("Weight var is a dummy variable, effective obs =") :
		I_("Weight var is a dummy variable, effective obs ="),
		pmod->nobs);
    }

    if (RTF_FORMAT(prn->format)) pprintf(prn, "\\par\n");
    else pprintf(prn, "\n");
}

static void model_format_start (PRN *prn)
{
    if (TEX_FORMAT(prn->format)) {
	if (STANDALONE(prn->format)) {
	    pprintf(prn, "\\documentclass{article}\n"
		    "\\usepackage{dcolumn}\n");
#ifdef ENABLE_NLS
	    pprintf(prn, "\\usepackage[latin1]{inputenc}\n");
#endif
	    pprintf(prn, "\\begin{document}\n\n"
		    "\\thispagestyle{empty}\n");
	}
	pprintf(prn, "\\begin{center}\n");
    }

    else if (RTF_FORMAT(prn->format)) {
	pprintf(prn, "{\\rtf1\\par\n\\qc ");
	return;
    }
}

#define RTF_COEFF_ROW  "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                       "\\cellx500\\cellx1900\\cellx3300\\cellx4700\\cellx6100" \
                       "\\cellx7500\\cellx8000\n\\intbl"

#define RTF_DISCRETE_ROW  "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                       "\\cellx500\\cellx1900\\cellx3300\\cellx4700\\cellx6100" \
                       "cellx8000\n\\intbl"

#define RTF_SELST_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                       "\\cellx1333\\cellx2666\\cellx4000\\cellx5333" \
                       "\\cellx6666\\cellx8000\n\\intbl"


static void print_coeff_table_start (PRN *prn, int discrete)
{
    if (PLAIN_FORMAT(prn->format)) {
	if (discrete) {
	    pprintf(prn, _("      VARIABLE      COEFFICIENT        STDERROR"
			   "       T STAT       SLOPE\n"));
	    pprintf(prn, "                                                 "
		    "                 %s\n", _("(at mean)"));
	} else {
	    pprintf(prn, _("      VARIABLE      COEFFICIENT        STDERROR"
			   "       T STAT   2Prob(t > |T|)\n\n"));
	}
    }

    else if (TEX_FORMAT(prn->format)) {
	char pt = get_local_decpoint();

	pprintf(prn, "\\vspace{1em}\n\n"
		"\\begin{tabular*}{\\textwidth}"
		"{@{\\extracolsep{\\fill}}\n"
		"l%% col 1: varname\n"
		"  D{%c}{%c}{-1}%% col 2: coeff\n"
		"    D{%c}{%c}{-1}%% col 3: sderr\n"
		"      D{%c}{%c}{-1}%% col 4: t-stat\n"
		"        D{%c}{%c}{4}}%% col 5: p-value (or slope)\n"
		"%s &\n"
		"  \\multicolumn{1}{c}{%s} &\n"
		"    \\multicolumn{1}{c}{%s} &\n"
		"      \\multicolumn{1}{c}{%s} &\n"
		"        \\multicolumn{1}{c}{%s} \\\\[1ex]\n",
		pt, pt, pt, pt, pt, pt, pt, pt, I_("Variable"),
		I_("Coefficient"), I_("Std.\\ Error"), 
		I_("$t$-statistic"), 
		(discrete)? I_("Slope at mean"): I_("p-value"));
    }   

    else if (RTF_FORMAT(prn->format)) {
	if (discrete) {
	    pprintf(prn, "{" RTF_DISCRETE_ROW
		    " \\qr \\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\intbl \\row\n",
		    I_("Variable"), I_("Coefficient"), I_("Std. Error"), 
		    I_("t-statistic"), I_("Slope at mean"));
	} else {
	    pprintf(prn, "{" RTF_COEFF_ROW
		    " \\qr \\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\ql \\cell"
		    " \\intbl \\row\n",
		    I_("Variable"), I_("Coefficient"), I_("Std. Error"), 
		    I_("t-statistic"), I_("p-value"));
	}
    } 
}

static void print_coeff_table_end (PRN *prn)
{
    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "\n");
    }
    else if (TEX_FORMAT(prn->format)) {
	pprintf(prn, "\\end{tabular*}\n\n");
    }
    else if (RTF_FORMAT(prn->format)) {
	pprintf(prn, "}\n\n");
    }
}

static void model_format_end (PRN *prn)
{
    if (TEX_FORMAT(prn->format)) {
	pprintf(prn, "\n\\end{center}\n");
	if (STANDALONE(prn->format)) {
	    pprintf(prn, "\n\\end{document}\n"); 
	}
    }

    else if (RTF_FORMAT(prn->format)) {
	pprintf(prn, "\n}\n");
    }
}  

static int print_coefficients (const MODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    int i, err = 0, gotnan = 0;
    int ncoeff = pmod->list[0];

    if (pmod->ifc) {
	if (PLAIN_FORMAT(prn->format)) {
	    err = print_coeff(pdinfo, pmod, ncoeff, prn);
	}
	else if (TEX_FORMAT(prn->format)) {
	    err = tex_print_coeff(pdinfo, pmod, ncoeff, prn);
	}
	else if (RTF_FORMAT(prn->format)) {
	    err = rtf_print_coeff(pdinfo, pmod, ncoeff, prn);
	}
	if (err) gotnan = 1;
	ncoeff--;
    }

    for (i=2; i<=ncoeff; i++) {
	if (PLAIN_FORMAT(prn->format)) {
	    err = print_coeff(pdinfo, pmod, i, prn);
	}
	else if (TEX_FORMAT(prn->format)) {
	    err = tex_print_coeff(pdinfo, pmod, i, prn);
	}
	else if (RTF_FORMAT(prn->format)) {
	    err = rtf_print_coeff(pdinfo, pmod, i, prn);
	}
	if (err) gotnan = 1;
    }

    return gotnan;
} 

static void print_middle_table_start (PRN *prn)
{
    if (TEX_FORMAT(prn->format)) {
	char pt = get_local_decpoint();

	pprintf(prn, 
		"\\vspace{1em}\n\n"
		"\\begin{tabular}{lD{%c}{%c}{-1}}\n",
		pt, pt, pt, pt);
    }
}

static void print_middle_table_end (PRN *prn)
{
    if (TEX_FORMAT(prn->format)) {
	pprintf(prn, "\\end{tabular}\n\n");
    }
    else if (RTF_FORMAT(prn->format)) {
	pprintf(prn, "\\par\n");
    }
    else {
	pprintf(prn, "\n");
    }
}

static void r_squared_message (PRN *prn)
{
    pprintf(prn, "\n%s.\n",    
	    _("R-squared is computed as the square of the correlation "
	      "between observed and\nfitted values of the dependent variable"));
}

static void weighted_stats_message (PRN *prn)
{
    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "%s:\n\n", _("Statistics based on the weighted data"));
    }
    else if (TEX_FORMAT(prn->format)) {
	pprintf(prn, "\\vspace{1em}%s:\n\n", 
		I_("Statistics based on the weighted data"));
    } 
    else { /* RTF */
	pprintf(prn, "\\par \\qc\n%s:\n\n", 
		I_("Statistics based on the weighted data"));	
    }
}

static void original_stats_message (PRN *prn)
{
    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "%s:\n\n", _("Statistics based on the original data"));
    } 
    else if (TEX_FORMAT(prn->format)) {
	pprintf(prn, "\\vspace{1em}\n%s:\n\n", 
		I_("Statistics based on the original data"));
    } 
    else { /* RTF */
	pprintf(prn, "\\par \\qc\n%s:\n\n", 
		I_("Statistics based on the original data"));
    }
}

static void rho_differenced_stats_message (PRN *prn)
{
    if (PLAIN_FORMAT(prn->format)) {    
	pprintf(prn, "%s:\n\n", _("Statistics based on the rho-differenced data"));
    } 
    else if (TEX_FORMAT(prn->format)) {
	pprintf(prn, "\\vspace{1em}\n%s:\n\n", 
		I_("Statistics based on the rho-differenced data"));
    } 
    else { /* RTF */
	pprintf(prn, "\\par \\qc\n%s:\n\n", 
		I_("Statistics based on the rho-differenced data"));
    }	
}

static void print_whites_results (const MODEL *pmod, PRN *prn)
{
    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "\n%s: TR^2 = %f,\n", _("Test statistic"), 
		pmod->rsq * pmod->nobs);
	pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n", 
		_("with p-value"), _("Chi-square"), 
		pmod->ncoeff - 1, pmod->rsq * pmod->nobs,
		chisq(pmod->rsq * pmod->nobs, pmod->ncoeff - 1)); 
    }

    if (RTF_FORMAT(prn->format)) { /* FIXME */
	pprintf(prn, "\\par \\ql\n%s: TR{\\super 2} = %f,\n", I_("Test statistic"), 
		pmod->rsq * pmod->nobs);
	pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n", 
		I_("with p-value"), I_("Chi-square"), 
		pmod->ncoeff - 1, pmod->rsq * pmod->nobs,
		chisq(pmod->rsq * pmod->nobs, pmod->ncoeff - 1)); 
    }

    else if (TEX_FORMAT(prn->format)) {
	pprintf(prn, "\n%s: $TR^2$ = %f,\n", I_("Test statistic"), 
		pmod->rsq * pmod->nobs);
	pprintf(prn, "%s = $P$($\\chi^2(%d)$ > %f) = %f\n\n",
		I_("with p-value"), 
		pmod->ncoeff - 1, pmod->rsq * pmod->nobs,
		chisq(pmod->rsq * pmod->nobs, pmod->ncoeff - 1)); 
    }
}

/**
 * printmodel:
 * @pmod: pointer to gretl model.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Print to @prn the estimates in @pmod plus associated statistics.
 * 
 * Returns: 0 on success, 1 if some of the values to print were NAN.
 */

int printmodel (const MODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    int gotnan = 0;
    int is_discrete = (pmod->ci == PROBIT || pmod->ci == LOGIT);

    if (prn->format != GRETL_PRINT_FORMAT_PLAIN) {
	model_format_start(prn);
    }

    print_model_heading (pmod, pdinfo, prn);

    print_coeff_table_start (prn, is_discrete);

    gotnan = print_coefficients (pmod, pdinfo, prn);

    if (pmod->ci == AR) {
	print_rho_terms (pmod, prn); 
    }

    print_coeff_table_end (prn);

    if (pmod->aux == AUX_ARCH || pmod->aux == AUX_ADF) {
	goto close_format;
    }

    if (is_discrete) {
	print_discrete_statistics(pmod, pdinfo, prn);
	goto close_format;
    }

    if (pmod->ci == LAD) {
	print_middle_table_start(prn);
	depvarstats(pmod, prn);
	ladstats(pmod, prn);
	print_middle_table_end(prn);
	goto close_format;
    }

    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG || pmod->aux == AUX_AR) {
	print_middle_table_start(prn);
	rsqline(pmod, prn);
	print_middle_table_end(prn);
	goto close_format;
    }

    if (!pmod->ifc && PLAIN_FORMAT(prn->format)) noconst(prn);
    
    if (pmod->aux == AUX_WHITE) { 
	rsqline(pmod, prn);
	print_whites_results(pmod, prn);
	goto close_format;
    }

    if (pmod->ci == OLS || pmod->ci == VAR || pmod->ci == TSLS
	|| pmod->ci == HCCM || pmod->ci == POOLED  
	|| (pmod->ci == WLS && pmod->wt_dummy)) {

	print_middle_table_start(prn);
	depvarstats(pmod, prn);
	if (essline(pmod, prn, 0)) {
	    print_middle_table_end(prn);
	    goto close_format;
	}
	rsqline(pmod, prn);
	Fline(pmod, prn);
	if (pmod->ci == OLS || (pmod->ci == WLS && pmod->wt_dummy)) {
	    if (pmod->ldepvar) dhline(pmod, prn);
	    else dwline(pmod, prn);
	}
	/* FIXME -- check output below */
	if (pmod->ci == HCCM || pmod->ci == TSLS) dwline(pmod, prn);
	print_middle_table_end(prn);

	if (pmod->ci == TSLS && PLAIN_FORMAT(prn->format)) {
	    r_squared_message(prn);
	}
    }

    else if (pmod->ci == HSK || pmod->ci == ARCH ||
	     (pmod->ci == WLS && pmod->wt_dummy == 0)) {

	weighted_stats_message(prn);
	print_middle_table_start(prn);
	if (essline(pmod, prn, 1)) {
	    print_middle_table_end(prn);
	    goto close_format;
	}
	rsqline(pmod, prn);
	Fline(pmod, prn);
	dwline(pmod, prn);
	print_middle_table_end(prn);

	original_stats_message(prn);
	print_middle_table_start(prn);
	depvarstats(pmod, prn);
	if (essline(pmod, prn, 0)) {
	    print_middle_table_end(prn);
	    goto close_format;
	}
	print_middle_table_end(prn);
    }

    else if (pmod->ci == AR || pmod->ci == CORC || pmod->ci == HILU) {

	rho_differenced_stats_message(prn);

	print_middle_table_start(prn);
	if (essline(pmod, prn, 0)) {
	    print_middle_table_end(prn);
	    goto close_format;
	}
	rsqline(pmod, prn);
	Fline(pmod, prn);
	dwline(pmod, prn);
	print_middle_table_end(prn);
    }

    if (PLAIN_FORMAT(prn->format)) print_aicetc(pmod, prn);
    else if (TEX_FORMAT(prn->format)) tex_print_aicetc(pmod, prn);
    else if (RTF_FORMAT(prn->format)) rtf_print_aicetc(pmod, prn);

    pmax_line(pmod, pdinfo, prn);
    
    print_model_tests(pmod, prn);

 close_format:
    if (!PLAIN_FORMAT(prn->format)) {
	model_format_end(prn);
    }

    return gotnan;
}

/* ....................................................... */

static void print_aicetc (const MODEL *pmod, PRN *prn)
{
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG ||
	pmod->aux == AUX_COINT || pmod->aux == AUX_WHITE ||
	pmod->aux == AUX_AR) return;

    if (pmod->dfd == 0) {
	pprintf(prn, "\n");
	return;
    }

    pprintf(prn, "  %s\n\n", _("MODEL SELECTION STATISTICS"));	
    pprintf(prn, "  SGMASQ    %#11g     AIC       %#11g     FPE       %#11g\n"
	    "  HQ        %#11g     SCHWARZ   %#11g     SHIBATA   %#11g\n"
	    "  GCV       %#11g",
	    pmod->criterion[0], pmod->criterion[1], 
	    pmod->criterion[2], pmod->criterion[3], 
	    pmod->criterion[4], pmod->criterion[5], pmod->criterion[6]);
    if (pmod->criterion[7] > 0.0) pprintf(prn, "     RICE      %#11g\n", 
					  pmod->criterion[7]);
    else pprintf(prn, "     RICE        %s\n", _("undefined"));
    pprintf(prn, "\n");
}

/* ......................................................... */ 

static void make_cname (const char *orig, char *cname)
{
    char *p;
    unsigned char c;

    if (orig == NULL || strlen(orig) == 0) return;

    p = strrchr(orig, '_');
    if (p == NULL) {
	strcpy(cname, orig);
	return;
    }

    c = (unsigned char) *(p + 1);

    if (isdigit(c)) {
	int lag = atoi(++p);

	sprintf(cname, "ut^2(-%d)", lag);
    }
}

/* ......................................................... */ 

static void print_pval_str (double pval, char *str)
{
    if (pval < .00001) {
	sprintf(str, "< %.5f", 0.00001);
    } else {
	sprintf(str, "%f", pval);
    }
}

/* ......................................................... */ 

static int print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			int c, PRN *prn)
{
    double t, pvalue = 999;
    int gotnan = 0;
    int do_pval = (pmod->ci != LOGIT && pmod->ci != PROBIT);
    char varname[12];

    /* special treatment for ARCH model coefficients */
    if (pmod->aux == AUX_ARCH) {
	make_cname(pdinfo->varname[pmod->list[c]], varname);
    } else {
	strcpy(varname, pdinfo->varname[pmod->list[c]]);
    }

    pprintf(prn, " %3d) %8s ", pmod->list[c], varname);

    _bufspace(3, prn);

    /* print coeff value if well-defined */
    if (isnan(pmod->coeff[c-1])) {
	pprintf(prn, "%16s", _("undefined"));
	gotnan = 1;
    } else {
	gretl_print_value (pmod->coeff[c-1], prn);
    }

    _bufspace(2, prn);

    /* get out if std error is undefined */
    if (isnan(pmod->sderr[c-1])) {
	pprintf(prn, "%16s\n", _("undefined"));
	return 1;
    }

    gretl_print_value (pmod->sderr[c-1], prn); 

    /* std error is well-defined, but is it positive? */
    if (pmod->sderr[c-1] > 0.) {
	t = pmod->coeff[c-1] / pmod->sderr[c-1];
	pprintf(prn, " %7.3f", t);
	if (pmod->aux == AUX_ADF) {
	    do_pval = 0;
	    pprintf(prn, "%12s", _("unknown"));
	}
	if (do_pval) {
	    char pvalstr[16];

	    pvalue = tprob(t, pmod->dfd);
	    print_pval_str(pvalue, pvalstr);
	    pprintf(prn, "%12s", pvalstr);
	}
    } else if (do_pval) { /* zero standard error */
	do_pval = 0;
	pprintf(prn, "     %12s", _("undefined"));
    }

    if (do_pval) {
	if (pvalue < 0.01) pprintf(prn, " ***");
	else if (pvalue < 0.05) pprintf(prn, " **");
	else if (pvalue < 0.10) pprintf(prn, " *");
    } 
    else if (pmod->list[c] != 0 && pmod->slope != NULL) { /* LOGIT, PROBIT */
	gretl_print_value (pmod->slope[c-1], prn);
    }

    pprintf(prn, "\n");

    return gotnan;
}

/* ......................................................... */ 

static void rtf_print_double (double xx, PRN *prn)
{
    pprintf(prn, " \\qc %#.*g\\cell", GRETL_DIGITS, xx);
}

/* ......................................................... */ 

static int rtf_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			    int c, PRN *prn)
{
    double t, pvalue = 999;
    int gotnan = 0;
    int do_pval = (pmod->ci != LOGIT && pmod->ci != PROBIT);
    char varname[12];

    /* special treatment for ARCH model coefficients */
    if (pmod->aux == AUX_ARCH) {
	make_cname(pdinfo->varname[pmod->list[c]], varname);
    } else {
	strcpy(varname, pdinfo->varname[pmod->list[c]]);
    }    

    pprintf(prn, RTF_COEFF_ROW);

    pprintf(prn, " \\qr %d\\cell \\ql %s\\cell", pmod->list[c], 
	    varname);

    if (isnan(pmod->coeff[c-1])) {
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	gotnan = 1;
    } else {
	rtf_print_double(pmod->coeff[c-1], prn);
    }

    if (isnan(pmod->sderr[c-1])) {
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	goto rtf_coeff_getout;
    } 

    rtf_print_double(pmod->sderr[c-1], prn); 

    if (pmod->sderr[c-1] > 0.) {
	t = pmod->coeff[c-1] / pmod->sderr[c-1];
	pprintf(prn, " \\qc %.4f\\cell", t);
	if (pmod->aux == AUX_ADF) {
	    do_pval = 0;
	    pprintf(prn, " \\qc %s\\cell", I_("unknown"));
	}
	if (do_pval) {
	    char pvalstr[16];

	    pvalue = tprob(t, pmod->dfd);
	    print_pval_str(pvalue, pvalstr);
	    pprintf(prn, " \\qc %s\\cell", pvalstr);
	}	
    } else if (do_pval) { /* zero standard error */
	do_pval = 0;
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
    }

    if (do_pval) {
	if (pvalue < 0.01) 
	    pprintf(prn, " \\ql ***\\cell");
	else if (pvalue < 0.05) 
	    pprintf(prn, " \\ql **\\cell");
	else if (pvalue < 0.10) 
	    pprintf(prn, " \\ql *\\cell");
	else 
	    pprintf(prn, " \\ql \\cell");
    } else if (pmod->list[c] != 0 && pmod->slope != NULL) { /* LOGIT, PROBIT */
	rtf_print_double(pmod->slope[c-1], prn);
    }

 rtf_coeff_getout:
    pprintf(prn, " \\intbl \\row\n");

    return gotnan;
}

/* ......................................................... */ 

static void print_rho (const ARINFO *arinfo, int c, int dfd, PRN *prn)
{
    char ustr[32];
    double xx = arinfo->rho[c] / arinfo->sderr[c];

    if (PLAIN_FORMAT(prn->format)) {
	sprintf(ustr, "u_%d", arinfo->arlist[c]);
	pprintf(prn, "%14s ", ustr); 
	_bufspace(3, prn);
	gretl_print_value (arinfo->rho[c], prn);
	_bufspace(2, prn);
	gretl_print_value (arinfo->sderr[c], prn); 
	pprintf(prn, " %7.3f %11f\n", xx, tprob(xx, dfd));	
    }

    else if (TEX_FORMAT(prn->format)) {
	char coeff[32], sderr[32];

	tex_dcolumn_double(arinfo->rho[c], coeff);
	tex_dcolumn_double(arinfo->sderr[c], sderr);

	sprintf(ustr, "$\\hat{u}_{t-%d}$", arinfo->arlist[c]);

	pprintf(prn, "%s &\n"
		"  %s &\n"
		"    %s &\n"
		"      %.4f &\n"
		"        %.4f \\\\\n",  
		ustr,
		coeff,
		sderr,
		arinfo->rho[c] / arinfo->sderr[c],
		tprob(xx, dfd));
    }
}

/* ......................................................... */ 

static void print_rho_terms (const MODEL *pmod, PRN *prn)
{
    int i, dfd;
    double xx = 0.0;

    if (pmod->arinfo == NULL || 
	pmod->arinfo->arlist == NULL ||
	pmod->arinfo->rho == NULL ||
	pmod->arinfo->sderr == NULL) {
	return;
    }

    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "\n%s:\n\n", _("Estimates of the AR coefficients"));
    }

    dfd = pmod->dfd + (pmod->ncoeff - pmod->arinfo->arlist[0]);

    for (i=1; i<=pmod->arinfo->arlist[0]; i++) {
	print_rho(pmod->arinfo, i, dfd, prn);
	xx += pmod->arinfo->rho[i]; 
    }

    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "\n%s = %#g\n\n", _("Sum of AR coefficients"), xx);
    }
}

static void tex_float_str (double x, char *str)
{
    if (x < 0.) {
	strcpy(str, "$-$");
	sprintf(str + 3, "%#.*g", GRETL_DIGITS, -x);
    } else {
	sprintf(str, "%#.*g", GRETL_DIGITS, x);
    }
}

static void print_discrete_statistics (const MODEL *pmod, 
				       const DATAINFO *pdinfo,
				       PRN *prn)
{
    int i;
    double pc_correct = 100 * (double) pmod->correct / pmod->nobs;

    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "\n");
	pprintf(prn, "%s %s = %.3f\n", _("Mean of"), 
		pdinfo->varname[pmod->list[1]], pmod->ybar);
	pprintf(prn, "%s = %d (%.1f%%)\n", _("Number of cases 'correctly predicted'"), 
		pmod->correct, pc_correct);
	pprintf(prn, "f(beta'x) %s = %.3f\n", _("at mean of independent vars"), 
		pmod->sdy);
	pprintf(prn, "%s = %.3f\n", _("Log-likelihood"), pmod->lnL);
	if (pmod->aux != AUX_OMIT && pmod->aux != AUX_ADD) {
	    i = pmod->ncoeff - 1;
	    pprintf(prn, "%s: %s(%d) = %.3f (%s %f)\n",
		    _("Likelihood ratio test"), _("Chi-square"), 
		    i, pmod->chisq, _("p-value"), chisq(pmod->chisq, i));
	}
	pprintf(prn, "\n");
    }

    else if (RTF_FORMAT(prn->format)) {
	pprintf(prn, "\n");
	pprintf(prn, "\\par %s %s = %.3f\n", I_("Mean of"), 
		pdinfo->varname[pmod->list[1]], pmod->ybar);
	pprintf(prn, "\\par %s = %d (%.1f%%)\n", I_("Number of cases 'correctly predicted'"), 
		pmod->correct, pc_correct);
	pprintf(prn, "\\par f(beta'x) %s = %.3f\n", I_("at mean of independent vars"), 
		pmod->sdy);
	pprintf(prn, "\\par %s = %.3f\n", I_("Log-likelihood"), pmod->lnL);
	if (pmod->aux != AUX_OMIT && pmod->aux != AUX_ADD) {
	    i = pmod->ncoeff - 1;
	    pprintf(prn, "\\par %s: %s(%d) = %.3f (%s %f)\\par\n",
		    I_("Likelihood ratio test"), I_("Chi-square"), 
		    i, pmod->chisq, I_("p-value"), chisq(pmod->chisq, i));
	}
	pprintf(prn, "\n");
    }

    else if (TEX_FORMAT(prn->format)) {
	char lnlstr[16];

	pprintf(prn, "\\vspace{1em}\n\\begin{raggedright}\n");
	pprintf(prn, "%s %s = %.3f\\\\\n", I_("Mean of"), 
		pdinfo->varname[pmod->list[1]], pmod->ybar);
	pprintf(prn, "%s = %d (%.1f percent)\\\\\n", 
		I_("Number of cases `correctly predicted'"), 
		pmod->correct, pc_correct);
	pprintf(prn, "$f(\\beta'x)$ %s = %.3f\\\\\n", I_("at mean of independent vars"), 
		pmod->sdy);
	tex_float_str(pmod->lnL, lnlstr);
	pprintf(prn, "%s = %s\\\\\n", I_("Log-likelihood"), lnlstr);
	if (pmod->aux != AUX_OMIT && pmod->aux != AUX_ADD) {
	    i = pmod->ncoeff - 1;
	    pprintf(prn, "%s: $\\chi^2_{%d}$ = %.3f (%s %f)\\\\\n",
		    I_("Likelihood ratio test"), 
		    i, pmod->chisq, I_("p-value"), chisq(pmod->chisq, i));
	}
	pprintf(prn, "\\end{raggedright}\n");
    }
}

static void tex_print_aicetc (const MODEL *pmod, PRN *prn)
{
    pprintf(prn, 
	    "\\vspace{1em}\n\n"
	    "%s\n\n"
	    "\\vspace{1em}\n\n"
	    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrlrlr}\n",
	    I_("Model selection statistics"));
    pprintf(prn, 
	    "\\textsc{sgmasq}  & %#g & "  
	    "\\textsc{aic}     & %#g & "  
	    "\\textsc{fpe}     & %#g \\\\\n"
	    "\\textsc{hq}      & %#g & "
	    "\\textsc{schwarz} & %#g & "  
	    "\\textsc{shibata} & %#g \\\\\n"
	    "\\textsc{gcv}     & %#g & "  
	    "\\textsc{rice}    & %#g\n",
	    pmod->criterion[0], pmod->criterion[1], pmod->criterion[2],
	    pmod->criterion[3], pmod->criterion[4], pmod->criterion[5],
	    pmod->criterion[6], pmod->criterion[7]);
    pprintf(prn, "\\end{tabular*}\n\n");
}

/* ....................................................... */

static void rtf_print_aicetc (const MODEL *pmod, PRN *prn)
{
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG ||
	pmod->aux == AUX_COINT || pmod->aux == AUX_WHITE ||
	pmod->aux == AUX_AR) return;

    if (pmod->dfd == 0) return;

    pprintf(prn, "\\par \\qc %s\\par\n\n", I_("Model Selection Statistics"));
    pprintf(prn, "{" RTF_SELST_ROW
	    " \\ql SGMASQ\\cell"
	    " \\qc %#g\\cell"
	    " \\ql AIC\\cell"
	    " \\qc %#g\\cell"
	    " \\ql FPE\\cell"
	    " \\qc %#g\\cell"
	    " \\intbl \\row\n"
	    RTF_SELST_ROW
	    " \\ql HQ\\cell"
	    " \\qc %#g\\cell"
	    " \\ql SCHWARZ\\cell"
	    " \\qc %#g\\cell"
	    " \\ql SHIBATA\\cell"
	    " \\qc %#g\\cell"
	    " \\intbl \\row\n"
	    RTF_SELST_ROW
	    " \\ql GCV\\cell"
	    " \\qc %#g\\cell"
	    " \\ql RICE\\cell",
	    pmod->criterion[0], pmod->criterion[1], 
	    pmod->criterion[2], pmod->criterion[3], 
	    pmod->criterion[4], pmod->criterion[5], pmod->criterion[6]);
    if (pmod->criterion[7] > 0.0) 
	pprintf(prn, " \\qc %#g\\cell", 
		pmod->criterion[7]);
    else
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
    pprintf(prn, " \\qr \\cell \\qr \\cell");

    pprintf(prn, " \\intbl \\row}\n\n");
}

/* ....................................................... */

static void mp_other_stats (const mp_results *mpvals, PRN *prn)
{
    char fstr[16];
    int len = 24;

    if (doing_nls()) len = 36;
    
    pprintf(prn, "%-*s", len, _("Standard error"));
    gretl_print_fullwidth_double(mpvals->sigma, GRETL_MP_DIGITS, prn);
    pprintf(prn, "\n");

    pprintf(prn, "%-*s", len, _("Error Sum of Squares"));
    gretl_print_fullwidth_double(mpvals->ess, GRETL_MP_DIGITS, prn);
    pprintf(prn, "\n");

    pprintf(prn, "%-*s", len, _("Unadjusted R-squared"));
    gretl_print_fullwidth_double(mpvals->rsq, GRETL_MP_DIGITS, prn);
    pprintf(prn, "\n");

    pprintf(prn, "%-*s", len, _("Adjusted R-squared"));
    gretl_print_fullwidth_double(mpvals->adjrsq, GRETL_MP_DIGITS, prn);
    pprintf(prn, "\n");

    sprintf(fstr, "F(%d, %d)", mpvals->dfn, mpvals->dfd);
    pprintf(prn, "%-*s", len, fstr);
    gretl_print_fullwidth_double(mpvals->fstt, GRETL_MP_DIGITS, prn);
    pprintf(prn, "\n");
}

/* ....................................................... */

static void print_mpvals_coeff (const mp_results *mpvals, 
				int c, PRN *prn)
{
    int cnum = c - 1;

    pprintf(prn, " %3d) %8s ", mpvals->varlist[c], mpvals->varnames[c]);

    if (mpvals->ifc && c == mpvals->ncoeff) cnum = 0;

    gretl_print_fullwidth_double(mpvals->coeff[cnum], 
				 GRETL_MP_DIGITS, prn);
    gretl_print_fullwidth_double(mpvals->sderr[cnum], 
				 GRETL_MP_DIGITS, prn); 

    pprintf(prn, "\n");
}

/* ....................................................... */

void print_mpols_results (const mp_results *mpvals, DATAINFO *pdinfo,
			  PRN *prn)
{
    int i, ncoeff;
    char startdate[9], enddate[9];
    int nobs = mpvals->t2 - mpvals->t1 + 1;

    ncoeff = mpvals->ncoeff;
    ntodate(startdate, mpvals->t1, pdinfo);
    ntodate(enddate, mpvals->t2, pdinfo);

    if (!PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "FIXME: this is still to be implemented!\n\n");
    }

    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, _("Multiple-precision OLS estimates using "
		       "the %d observations %s-%s\n"),
		nobs, startdate, enddate);
	pprintf(prn, "%s: %s\n\n", _("Dependent variable"),
		mpvals->varnames[0]);

	pprintf(prn, _("      VARIABLE         COEFFICIENT          "
		       "        STD. ERROR\n"));
    }

    if (mpvals->ifc) {
	print_mpvals_coeff(mpvals, ncoeff, prn);
	ncoeff--;
    }
    for (i=2; i<=ncoeff; i++) {
	print_mpvals_coeff(mpvals, i, prn);
    }
    pprintf(prn, "\n");

    mp_other_stats (mpvals, prn);
}
