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
#include "gretl_private.h"
#include "libset.h"

#ifndef cmplx
typedef struct _cmplx cmplx;
struct _cmplx {
    double r;
    double i;
};
#endif

#define PLAIN_FORMAT(c) (c == GRETL_PRINT_FORMAT_PLAIN)
#define RTF_FORMAT(c) (c == GRETL_PRINT_FORMAT_RTF)
#define TEX_FORMAT(c) (c == GRETL_PRINT_FORMAT_TEX || c == GRETL_PRINT_FORMAT_TEX_DOC)
#define STANDALONE(c) (c == GRETL_PRINT_FORMAT_TEX_DOC)

#define NO_RBAR_SQ(a) (a == AUX_SQ || a == AUX_LOG || a == AUX_WHITE || a == AUX_AR)

static int print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			int c, int longnames, PRN *prn);
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
static void print_arma_stats (const MODEL *pmod, PRN *prn);
static void print_arma_roots (const MODEL *pmod, PRN *prn);
static void print_tobit_stats (const MODEL *pmod, PRN *prn);

/* ......................................................... */ 

static void noconst (const MODEL *pmod, PRN *prn)
{
    if (na(pmod->rsq)) return;

    pputs(prn, _("The model has no constant term.\n"  
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
	pprintf(prn, "  %s = %.*g\n", _("Mean of dependent variable"), 
		GRETL_DIGITS, pmod->ybar);
	pprintf(prn, "  %s = %.*g\n", _("Standard deviation of dep. var."), 
		GRETL_DIGITS, pmod->sdy);
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

	if (pmod->ci == ARMA) {
	    pprintf(prn, "  %s = %.*g\n", _("Mean of residuals"), 
		    GRETL_DIGITS, gretl_model_get_double(pmod, "mean_error"));
	}

	pprintf(prn, "  %s = %.*g\n", _("Sum of squared residuals"), 
		GRETL_DIGITS, wt? pmod->ess_wt : pmod->ess);
	pprintf(prn, "  %s = %.*g\n", _("Standard error of residuals"), 
		GRETL_DIGITS, wt? pmod->sigma_wt : pmod->sigma);
    }

    else if (RTF_FORMAT(prn->format)) {
	if ((wt && pmod->ess_wt < 0) || (!wt && pmod->ess < 0)) {
	    char tmp[128];

	    sprintf(tmp, I_("Error sum of squares (%g) is not > 0"), 
		    (wt)? pmod->ess_wt : pmod->ess); 
	    pprintf(prn, "\\par \\ql %s\\par\n\n", tmp);
	    return 1;
	}

	if (pmod->ci == ARMA) {
	    pprintf(prn, RTFTAB "%s = %g\n", I_("Mean of residuals"), 
		    gretl_model_get_double(pmod, "mean_error"));
	}

	pprintf(prn, RTFTAB "%s = %g\n", I_("Sum of squared residuals"), 
		wt? pmod->ess_wt : pmod->ess);
	pprintf(prn, RTFTAB "%s = %g\n", I_("Standard error of residuals"), 
		wt? pmod->sigma_wt : pmod->sigma);
    }

    else if (TEX_FORMAT(prn->format)) {
	char x1str[32], x2str[32];

	if (pmod->ci == ARMA) {
	    tex_dcolumn_double(gretl_model_get_double(pmod, "mean_error"), x1str);
	    pprintf(prn, "%s & %s \\\\\n", I_("Mean of residuals"), 
		    x1str);
	}

	tex_dcolumn_double(pmod->ess, x1str);
	tex_dcolumn_double(pmod->sigma, x2str);
	pprintf(prn, "%s & %s \\\\\n%s ($\\hat{\\sigma}$) & %s \\\\\n",
		I_("Sum of squared residuals"), x1str,
		I_("Standard error of residuals"), x2str);
    }

    return 0;
}

/* ......................................................... */ 

static void rsqline (const MODEL *pmod, PRN *prn)
{
    if (PLAIN_FORMAT(prn->format)) { 
	if (!na(pmod->rsq)) {
	    pprintf(prn, "  %s = %.*g\n", _("Unadjusted R-squared"), 
		    GRETL_DIGITS, pmod->rsq);
	}
	if (!NO_RBAR_SQ(pmod->aux) && !na(pmod->adjrsq)) {
	    pprintf(prn, "  %s = %.*g\n", _("Adjusted R-squared"),  
		    GRETL_DIGITS, pmod->adjrsq);
	}
    }

    else if (RTF_FORMAT(prn->format)) {
	if (!na(pmod->rsq)) {
	    pprintf(prn, RTFTAB "%s = %g\n", I_("Unadjusted R{\\super 2}"), pmod->rsq);
	}
	if (!NO_RBAR_SQ(pmod->aux) && !na(pmod->adjrsq)) {
	    pprintf(prn, RTFTAB "%s = %g\n", I_("Adjusted R{\\super 2}"),  
		    pmod->adjrsq);
	}	
    }

    else if (TEX_FORMAT(prn->format)) {  
	char r2[32];

	if (!na(pmod->rsq)) {
	    tex_dcolumn_double(pmod->rsq, r2);
	    pprintf(prn, "%s & %s \\\\\n", I_("Unadjusted $R^2$"), r2);
	}
	if (!NO_RBAR_SQ(pmod->aux) && !na(pmod->adjrsq)) {
	    tex_dcolumn_double(pmod->adjrsq, r2);
	    pprintf(prn, "%s & %s \\\\\n", I_("Adjusted $\\bar{R}^2$"), r2);
	}
    }
}

/* ......................................................... */ 

static void info_stats_lines (const MODEL *pmod, PRN *prn)
{
    const double *crit = pmod->criterion;
    const char *info_str[] = {
	N_("Akaike information criterion"),
	N_("Schwarz Bayesian criterion")
    };
    const char *info_abbrev[] = {
	N_("AIC"),
	N_("BIC")
    };    

    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG ||
	pmod->aux == AUX_COINT || pmod->aux == AUX_WHITE ||
	pmod->aux == AUX_AR) {
	return;
    }

    if (na(crit[C_AIC]) || na(crit[C_BIC])) {
	return;
    }

    if (PLAIN_FORMAT(prn->format)) { 
	pprintf(prn, "  %s (%s) = %.*g\n", _(info_str[0]), _(info_abbrev[0]),
		GRETL_DIGITS, crit[C_AIC]);
	pprintf(prn, "  %s (%s) = %.*g\n", _(info_str[1]), _(info_abbrev[1]),
		GRETL_DIGITS, crit[C_BIC]);
    }

    else if (RTF_FORMAT(prn->format)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_(info_str[0]), crit[C_AIC]);
	pprintf(prn, RTFTAB "%s = %g\n", I_(info_str[1]), crit[C_BIC]);
    }

    else if (TEX_FORMAT(prn->format)) {  
	char cval[32];

	tex_dcolumn_double(crit[C_AIC], cval);
	pprintf(prn, "%s & %s \\\\\n", I_(info_str[0]), cval);
	tex_dcolumn_double(crit[C_BIC], cval);
	pprintf(prn, "%s & %s \\\\\n", I_(info_str[1]), cval);
    }
}

/* ......................................................... */ 

static void ladstats (const MODEL *pmod, PRN *prn)
{
    int utf = PLAIN_FORMAT(prn->format);

    if (TEX_FORMAT(prn->format)) {  
	char x1str[32], x2str[32];

	tex_dcolumn_double(pmod->rho, x1str); /* "rho" is abused here! */
	tex_dcolumn_double(pmod->ess, x2str);
	pprintf(prn, "%s & %s \\\\\n",
		I_("Sum of absolute residuals"), x1str); 
	pprintf(prn, "%s & %s \\\\\n",
		I_("Sum of squared residuals"), x2str); 
    }
    
    else {
	pprintf(prn, "  %s = %.*g\n", 
		(utf)? _("Sum of absolute residuals") :
		I_("Sum of absolute residuals"),
		GRETL_DIGITS, pmod->rho);
	pprintf(prn, "  %s = %.*g\n", 
		(utf)? _("Sum of squared residuals") :
		I_("Sum of squared residuals"),
		GRETL_DIGITS, pmod->ess);

	if (utf) {
	    int ladcode = gretl_model_get_int(pmod, "ladcode");

	    if (ladcode == 0) {
		pputs(prn, _("\nWarning: solution is probably not unique\n"));
	    }
	}
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
    if (pmod->ifc && pmod->ncoeff <= 2) {
	if (PLAIN_FORMAT(prn->format)) {
	    pprintf(prn, "  %s = %d\n", _("Degrees of freedom"), pmod->dfd);
	}
	else if (TEX_FORMAT(prn->format)) {
	    pprintf(prn, "%s & %d \\\\\n", I_("Degrees of freedom"), pmod->dfd);
	}
	else if (RTF_FORMAT(prn->format)) {
	    pprintf(prn, RTFTAB "%s = %d\n", I_("Degrees of freedom"), pmod->dfd);
	}
	return;
    }

    if (PLAIN_FORMAT(prn->format)) {
	char tmp[32];

	sprintf(tmp, "%s (%d, %d)", _("F-statistic"), pmod->dfn, pmod->dfd);
	if (na(pmod->fstt)) {
	    pprintf(prn, "  %s %s\n", tmp, _("undefined"));
	} else {
	    pprintf(prn, "  %s = %.*g", tmp, GRETL_DIGITS, pmod->fstt);
	    print_f_pval_str(fdist(pmod->fstt, pmod->dfn, pmod->dfd), prn);
	}
    }

    else if (TEX_FORMAT(prn->format)) {
	if (na(pmod->fstt)) {
	    pprintf(prn, "$F(%d, %d)$ & \\multicolumn{1}{c}{\\rm %s} \\\\\n", 
		    pmod->dfn, pmod->dfd, I_("undefined"));
	} else {
	    char x1str[32];

	    tex_dcolumn_double(pmod->fstt, x1str);
	    pprintf(prn, "$F(%d, %d)$ & %s \\\\\n", pmod->dfn, pmod->dfd, x1str);
	    print_f_pval_str(fdist(pmod->fstt, pmod->dfn, pmod->dfd), prn);
	}
    }

    else if (RTF_FORMAT(prn->format)) {
	char tmp[32];

	sprintf(tmp, "%s (%d, %d)", I_("F-statistic"), pmod->dfn, pmod->dfd);
	if (na(pmod->fstt)) {
	    pprintf(prn, RTFTAB "%s %s\n", tmp, I_("undefined"));
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
	    pprintf(prn, "  %s = %.*g\n", _("Durbin-Watson statistic"), 
		    GRETL_DIGITS, pmod->dw);
	    pprintf(prn, "  %s = %.*g\n", _("First-order autocorrelation coeff."), 
		    GRETL_DIGITS, pmod->rho);
	} 
    }

    else if (TEX_FORMAT(prn->format)) {
	char x1str[32], x2str[32];

	tex_dcolumn_double(pmod->dw, x1str);
	tex_dcolumn_double(pmod->rho, x2str);
	pprintf(prn, "%s & %s \\\\\n%s & %s \\\\\n",
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
    int ldepvar = gretl_model_get_int(pmod, "ldepvar");
    int T = pmod->nobs - 1;

    sderr = pmod->sderr[ldepvar - 2];

    if (pmod->ess <= 0.0 || (T * sderr * sderr) >= 1.0) return;

    h = pmod->rho * sqrt(T/(1 - T * sderr * sderr));

    if (PLAIN_FORMAT(prn->format)) {
        char tmp[128];

        sprintf(tmp, _("Durbin's h stat. %g"), h);
        pprintf(prn, "  %s\n", tmp);

        sprintf(tmp, _("(Using variable %d for h stat, with T' = %d)"), 
                pmod->list[ldepvar], T);
        pprintf(prn, "  %s\n", tmp);
    }

    else if (RTF_FORMAT(prn->format)) {
        char tmp[128];

        sprintf(tmp, I_("Durbin's h stat. %g"), h);
        pprintf(prn, RTFTAB "%s\n", tmp);

        sprintf(tmp, I_("(Using variable %d for h stat, with T' = %d)"), 
                pmod->list[ldepvar], T);
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

static int least_signif_coeff (const MODEL *pmod)
{
    int i, k = 0;
    double tstat, tmin = 4.0;
    
    for (i=pmod->ifc; i<pmod->ncoeff; i++) {
	tstat = fabs(pmod->coeff[i] / pmod->sderr[i]);
	if (tstat < tmin) {
	    tmin = tstat;
	    k = i;
	}
    }

    if (tprob(tmin, pmod->dfd) > .10) {
	return pmod->list[k+2];
    }

    return 0;
}

/* ......................................................... */ 

static void pval_max_line (const MODEL *pmod, const DATAINFO *pdinfo, 
			   PRN *prn)
{
    int k = pmod->ncoeff - pmod->ifc;

    if (k < 2) return;

    if ((k = least_signif_coeff(pmod))) {
	char tmp[128];

	sprintf(tmp, _("Excluding the constant, p-value was highest "
		       "for variable %d (%s)"), k, pdinfo->varname[k]);
	pprintf(prn, "%s\n\n", tmp);
    }
}

/* ......................................................... */

static const char *aux_string (int aux, int format)
{
    if (aux == AUX_SQ) {
	return N_("Auxiliary regression for non-linearity test "
		 "(squared terms)");
    } else if (aux == AUX_LOG) {
	return N_("Auxiliary regression for non-linearity test "
		 "(log terms)");
    } else if (aux == AUX_WHITE) {
	return N_("White's test for heteroskedasticity");
    } else if (aux == AUX_CHOW) {
	return N_("Augmented regression for Chow test");
    } else if (aux == AUX_COINT) {
	if (TEX_FORMAT(format)) return N_("Cointegrating regression -- ");
	else return N_("Cointegrating regression - ");
    } else if (aux == AUX_ADF) {
	if (TEX_FORMAT(format)) return N_("Augmented Dickey--Fuller regression");
	else return N_("Augmented Dickey-Fuller regression");
    } else if (aux == AUX_DF) {
	if (TEX_FORMAT(format)) return N_("Dickey--Fuller regression");
	else return N_("Dickey-Fuller regression");
    } else if (aux == AUX_KPSS) {
	return N_("KPSS regression");
    } else if (aux == AUX_RESET) {
	return N_("Auxiliary regression for RESET specification test");
    } else if (aux == AUX_GROUPWISE) {
	return N_("Groupwise heteroskedasticity");
    }

    else return "";
}

/* ......................................................... */

const char *estimator_string (int ci, int format)
{
    if (ci == OLS || ci == VAR) return N_("OLS");
    else if (ci == WLS) return N_("WLS"); 
    else if (ci == ARCH) return N_("WLS (ARCH)");
    else if (ci == TSLS) return N_("TSLS");
    else if (ci == HSK) return N_("Heteroskedasticity-corrected");
    else if (ci == AR) return N_("AR");
    else if (ci == LAD) return N_("LAD");
    else if (ci == HCCM) return N_("HCCM");
    else if (ci == PROBIT) return N_("Probit");
    else if (ci == LOGIT) return N_("Logit");
    else if (ci == TOBIT) return N_("Tobit");
    else if (ci == POOLED) return N_("Pooled OLS");
    else if (ci == NLS) return N_("NLS");
    else if (ci == LOGISTIC) return N_("Logistic");
    else if (ci == GARCH) return N_("GARCH");
    else if (ci == CORC) {
	if (TEX_FORMAT(format)) return N_("Cochrane--Orcutt");
	else return N_("Cochrane-Orcutt");
    }
    else if (ci == HILU) {
	if (TEX_FORMAT(format)) return N_("Hildreth--Lu");
	else return N_("Hildreth-Lu");
    }
    else if (ci == PWE) {
	if (TEX_FORMAT(format)) return N_("Prais--Winsten");
	else return N_("Prais-Winsten");
    }

    else return "";
}

static const char *my_estimator_string (const MODEL *pmod,
					int format)
{
    if (pmod->ci == ARMA) {
	return (pmod->list[0] > 4)? N_("ARMAX") : N_("ARMA");
    } else if (pmod->ci == WLS) {
	if (gretl_model_get_int(pmod, "iters")) {
	    return N_("Maximum Likelihood");
	} else {
	    return N_("WLS");
	}
    } else {
	return estimator_string(pmod->ci, format);
    } 
}

/* ......................................................... */

void get_test_type_string (const GRETLTEST *test, char *str, int format)
{
    if (PLAIN_FORMAT(format)) {
	if (test->param[0]) {
	    sprintf(str, _(test->type), test->param);
	} else {
	    strcpy(str, _(test->type));
	} 
    } else {
	if (test->param[0]) {
	    sprintf(str, I_(test->type), test->param);
	} else {
	    strcpy(str, I_(test->type));
	}
    }
}

void get_test_stat_string (const GRETLTEST *test, char *str, int format)
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
    case GRETL_TEST_RESET:
	if (tex) {
	    sprintf(str, "$F(%d, %d)$ = %g", test->dfn, test->dfd, test->value);
	} else {
	    sprintf(str, "F(%d, %d) = %g", test->dfn, test->dfd, test->value);
	}
	break;
    case GRETL_TEST_LMF:
	sprintf(str, "LMF = %g", test->value);
	break;
    case GRETL_TEST_HARVEY_COLLIER:
	if (tex) {
	    sprintf(str, "Harvey--Collier $t(%d)$ = %g", test->dfn, test->value);
	} else {
	    sprintf(str, "Harvey-Collier t(%d) = %g", test->dfn, test->value);
	}
	break;
    case GRETL_TEST_NORMAL_CHISQ:
	if (tex) {
	    sprintf(str, "$\\chi^2_2$ = %g", test->value); 
	} else {
	    sprintf(str, "%s(2) = %g", _("Chi-square"), test->value);
	}
	break;
    case GRETL_TEST_LR:
	if (tex) {
	    sprintf(str, "$\\chi^2_%d$ = %g", test->dfn, test->value); 
	} else {
	    sprintf(str, "%s(%d) = %g", _("Chi-square"), test->dfn, test->value);
	}
	break;
    default:
	*str = 0;
    }
}

void get_test_pval_string (const GRETLTEST *test, char *str, int format)
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
    case GRETL_TEST_RESET:
	if (tex) {
	    sprintf(str, "$P$($F(%d, %d) >$ %g) = %g", 
		    test->dfn, test->dfd, test->value, test->pvalue);
	} else {
	    sprintf(str, "P(F(%d, %d) > %g) = %g", 
		    test->dfn, test->dfd, test->value, test->pvalue);
	}
	break;
    case GRETL_TEST_LMF:
	if (tex) {
	    sprintf(str, "$P$($F(%d, %d) >$ %g) = %g", 
		    test->dfn, test->dfd, test->value, test->pvalue);
	} else {
	    sprintf(str, "P(F(%d,%d) > %g) = %g", 
		    test->dfn, test->dfd, test->value, test->pvalue);
	}
	break;
    case GRETL_TEST_HARVEY_COLLIER:
	if (tex) {
	    sprintf(str, "$P$($t_{%d} >$ %g)  = %g", 
		    test->dfn, test->value, test->pvalue);
	} else {
	    sprintf(str, "P(t(%d) > %g) = %g", 
		    test->dfn, test->value, test->pvalue);
	}
	break;
    case GRETL_TEST_NORMAL_CHISQ:
    case GRETL_TEST_LR:
	sprintf(str, "%g", test->pvalue);
	break;
    default:
	*str = 0;
    }
}

/* ......................................................... */

static void print_model_tests (const MODEL *pmod, PRN *prn)
{
    int i;
    char test_str[64], pval_str[64], type_str[96];

    if (PLAIN_FORMAT(prn->format)) {
	for (i=0; i<pmod->ntests; i++) {
	    get_test_type_string(&pmod->tests[i], type_str, prn->format);
	    get_test_stat_string(&pmod->tests[i], test_str, prn->format);
	    get_test_pval_string(&pmod->tests[i], pval_str, prn->format);
	    pprintf(prn, "%s -\n"
		    "  %s: %s\n"
		    "  %s: %s\n"
		    "  %s = %s\n\n",
		    type_str, 
		    _("Null hypothesis"), _((pmod->tests[i]).h_0), 
		    _("Test statistic"), test_str, 
		    _("with p-value"), pval_str);
	}
    }

    else if (TEX_FORMAT(prn->format)) {
	if (pmod->ntests > 0) {
	    pputs(prn, "\\vspace{1em}\n\\begin{raggedright}\n");
	    for (i=0; i<pmod->ntests; i++) {
		if (i > 0) {
		    pputs(prn, "\\vspace{1ex}\n");;
		}
		get_test_type_string(&pmod->tests[i], type_str, prn->format);
		get_test_stat_string(&pmod->tests[i], test_str, prn->format);
		get_test_pval_string(&pmod->tests[i], pval_str, prn->format);
		pprintf(prn, "%s --\\\\\n"
			"\\quad %s: %s\\\\\n"
			"\\quad %s: %s\\\\\n"
			"\\quad %s = %s\\\\\n",
			type_str, 
			I_("Null hypothesis"), I_((pmod->tests[i]).h_0), 
			I_("Test statistic"), test_str, 
			I_("with p-value"), pval_str);
	    }
	    pputs(prn, "\\end{raggedright}\n");
	}
    }

    else if (RTF_FORMAT(prn->format)) {
	for (i=0; i<pmod->ntests; i++) {
	    get_test_type_string(&pmod->tests[i], type_str, prn->format);
	    get_test_stat_string(&pmod->tests[i], test_str, prn->format);
	    get_test_pval_string(&pmod->tests[i], pval_str, prn->format);
	    pprintf(prn, "\\par \\ql %s -\\par\n"
		    " %s: %s\\par\n"
		    " %s: %s\\par\n"
		    " %s = %s\\par\n\n",
		    type_str, 
		    I_("Null hypothesis"), I_((pmod->tests[i]).h_0),
		    I_("Test statistic"), test_str, 
		    I_("with p-value"), pval_str);
	}
    }
}

/* ......................................................... */

static int 
print_tsls_instruments (const int *list, const DATAINFO *pdinfo, PRN *prn)
{
    int i, j, pos = 0;
    int ccount = 0;
    char vname[16];
    int tex = TEX_FORMAT(prn->format);
    int utf = PLAIN_FORMAT(prn->format);

    if (utf) {
	pprintf(prn, "%s: ", _("Instruments"));
    } else {
	pprintf(prn, "%s: ", I_("Instruments"));
    }

    ccount += strlen(_("Instruments") + 2);

    for (i=2; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    pos = i;
	    continue;
	}
	if (pos && list[i] > 0) {
	    int dup = 0;

	    for (j=2; j<pos; j++) {
		if (list[i] == list[j]) {
		    dup = 1;
		    break;
		}
	    }
	    if (!dup) {
		if (tex) {
		    tex_escape(vname, pdinfo->varname[list[i]]);
		} else {
		    strcpy(vname, pdinfo->varname[list[i]]);
		}
		pprintf(prn, "%s ", vname);
		ccount += strlen(vname) + 1;
		if (ccount >= 76) {
		    if (tex) {
			pputs(prn, "\\\\\n");
		    } else if (RTF_FORMAT(prn->format)) {
			pputs(prn, "\\par\n");
		    } else {
			pputs(prn, "\n  "); 
		    }
		    ccount = 0;
		}
	    }
	}
    }

    if (ccount > 0) {
	if (tex) {
	    pputs(prn, "\\\\\n");
	} else if (RTF_FORMAT(prn->format)) {
	    pputs(prn, "\\par\n");
	} else {
	    pputs(prn, "\n");
	}
    }

    return 0;
}

/* ......................................................... */

static void hac_vcv_line (const MODEL *pmod, PRN *prn)
{
    int lag;

    if (pmod->aux == AUX_SCR) {
	lag = pmod->order;
    } else {
	lag = gretl_model_get_int(pmod, "hac_lag");
    }

    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, _("Serial correlation-robust standard errors, "
		       "lag order %d\n"), lag);
    } else {
	pprintf(prn, I_("Serial correlation-robust standard errors, "
		       "lag order %d\n"), lag);
    }	
}

static void hc_vcv_line (const MODEL *pmod, PRN *prn)
{
    int v = gretl_model_get_int(pmod, "hc_version");

    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, _("Heteroskedasticity-robust standard errors, "
		       "variant HC%d\n"), v);
    } else {
	pprintf(prn, I_("Heteroskedasticity-robust standard errors, "
			"variant HC%d\n"), v);
    }	
}

static void garch_vcv_line (const MODEL *pmod, PRN *prn)
{
    int v = gretl_model_get_int(pmod, "garch_vcv");
    int tex = TEX_FORMAT(prn->format);
    int utf = PLAIN_FORMAT(prn->format);
    const char *vcvstr = NULL;

    switch (v) {
    case VCV_HESSIAN:
	vcvstr = N_("Standard errors based on Hessian");
	break;
    case VCV_IM:
	vcvstr = N_("Standard errors based on Information Matrix");
	break;
    case VCV_OP:
	vcvstr = N_("Standard errors based on Outer Products matrix");
	break;
    case VCV_QML:
	vcvstr = N_("QML standard errors");
	break;
    case VCV_BW:
	if (tex) {
	    vcvstr = N_("Bollerslev--Wooldridge standard errors");
	} else {
	    vcvstr = N_("Bollerslev-Wooldridge standard errors");
	}
	break;
    default:
	break;
    }

    if (vcvstr != NULL) {
	pprintf(prn, "%s\n", (utf)? _(vcvstr) : I_(vcvstr));
    }
}

/* ......................................................... */

static void model_print_newline (PRN *prn)
{
    if (TEX_FORMAT(prn->format)) {
	pputs(prn, "\\\\\n");
    } else if (RTF_FORMAT(prn->format)) {
	pputs(prn, "\\par\n");
    } else {
	pputc(prn, '\n');
    }
}

/* ......................................................... */

static void print_model_heading (const MODEL *pmod, 
				 const DATAINFO *pdinfo, 
				 PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN], vname[16];
    int t1 = pmod->t1, t2 = pmod->t2;
    int tex = TEX_FORMAT(prn->format);
    int utf = PLAIN_FORMAT(prn->format);

    ntodate(startdate, t1, pdinfo);
    ntodate(enddate, t2, pdinfo);

    switch (pmod->aux) {
    case AUX_SQ:
    case AUX_LOG:
    case AUX_WHITE:
    case AUX_CHOW:
    case AUX_COINT:
    case AUX_ADF:
    case AUX_DF:
    case AUX_KPSS:
    case AUX_RESET:
    case AUX_GROUPWISE:
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
	    pprintf(prn, "\n%s ", _("Breusch-Godfrey test for"));
	} else {
	    pprintf(prn, "\n%s ", I_("Breusch-Godfrey test for"));
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
    case AUX_SYS:
	pprintf(prn, "\n%s %d: ", 
		(utf)? _("Equation") : I_("Equation"), pmod->ID + 1);
	break;	
    case VAR:
	pprintf(prn, "\n%s %d: ", 
		(utf)? _("Equation") : I_("Equation"), pmod->ID);
	break;
    case AUX_ADD:
    default:
	if (pmod->ID < 0) pputs(prn, "\n");
	else if (pmod->name) {
	    pprintf(prn, "\n%s:\n", pmod->name);
	} else {
	    pprintf(prn, "\n%s %d: ", (utf)? _("Model") : I_("Model"), pmod->ID);
	}
	break;
    }

    if (pmod->aux == AUX_SYS) {
	pprintf(prn, (utf)?
		_("%s estimates using the %d observations %s%s%s") :
		I_("%s estimates using the %d observations %s%s%s"),
		_(gretl_system_short_string(pmod)),
		pmod->nobs, startdate, (tex)? "--" : "-", enddate);
    } else if (!dataset_is_panel(pdinfo)) {
	if (pmod->missmask != NULL) {
	    int mc = model_missval_count(pmod);

	    pprintf(prn, (utf)?
		    _("%s estimates using %d observations from %s%s%s") :
		    I_("%s estimates using %d observations from %s%s%s"),
		    _(my_estimator_string(pmod, prn->format)), 
		    pmod->nobs, startdate, (tex)? "--" : "-", enddate);
	    model_print_newline(prn);
	    pprintf(prn, "%s: %d",
		    (utf)? _("Missing or incomplete observations dropped") :
		    I_("Missing or incomplete observations dropped"), mc);
	} else {
	    pprintf(prn, (utf)?
		    _("%s estimates using the %d observations %s%s%s") :
		    I_("%s estimates using the %d observations %s%s%s"),
		    _(my_estimator_string(pmod, prn->format)), 
		    pmod->nobs, startdate, (tex)? "--" : "-", enddate);
	}
    } else {
	int effn = gretl_model_get_int(pmod, "n_included_units");

	pprintf(prn, (utf)?
		_("%s estimates using %d observations") :
		I_("%s estimates using %d observations"),
		_(my_estimator_string(pmod, prn->format)), 
		pmod->nobs);
	if (effn > 0) {
	    model_print_newline(prn);
	    pprintf(prn, (utf)? _("Included %d cross-sectional units") :
		    I_("Included %d cross-sectional units"), effn);
	}
    }

    model_print_newline(prn);

    /* special names for dependent variable in cases of certain sorts
       of auxiliary regressions */
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG) {
	pprintf(prn, "%s: %s", 
		(utf)? _("Dependent variable") : I_("Dependent variable"),
		(tex)? "$\\hat{u}$" : "uhat");
    }
    else if (pmod->aux == AUX_WHITE) {
	pprintf(prn, "%s: %s", 
		(utf)? _("Dependent variable") : I_("Dependent variable"),
		(tex)? "$\\hat{u}^2$" : "uhat^2");
    }
    else if (pmod->aux == AUX_ARCH) {
	pprintf(prn, "%s: %s", 
		(utf)? _("Dependent variable") : I_("Dependent variable"),
		(tex)? "$u_t^2$" : "ut^2");
    }
    else if (pmod->ci == NLS) {
	if (tex) tex_escape(vname, pmod->params[0]);
	pprintf(prn, "%s: %s", 
		(utf)? _("Dependent variable") : I_("Dependent variable"),
		(tex)? vname : pmod->params[0]);
    }
    else { 
	int v = (pmod->ci == ARMA || pmod->ci == GARCH)? 
	    pmod->list[4] : pmod->list[1];

	if (tex) tex_escape(vname, pdinfo->varname[v]);
	pprintf(prn, "%s: %s", 
		(utf)? _("Dependent variable") : I_("Dependent variable"),
		(tex)? vname : pdinfo->varname[v]);
    }

    model_print_newline(prn);

    /* supplementary strings below the estimator and sample info */

    /* list of instruments for TSLS */
    if (pmod->ci == TSLS) {
	print_tsls_instruments (pmod->list, pdinfo, prn);
    }

    /* VCV variants */
    if (pmod->aux == AUX_SCR || gretl_model_get_int(pmod, "hac_lag")) {
	hac_vcv_line(pmod, prn);
    }
    else if (gretl_model_get_int(pmod, "hc")) {
	hc_vcv_line(pmod, prn);
    }
    else if (gretl_model_get_int(pmod, "garch_vcv")) {
	garch_vcv_line(pmod, prn);
    }

    /* WLS on panel data */
    else if (gretl_model_get_int(pmod, "unit_weights") && !pmod->aux) {
	if (tex) {
	    pputs(prn, "\\\\\n");
	}
	if (gretl_model_get_int(pmod, "iters")) {
	    pprintf(prn, (utf)? _("Allowing for groupwise heteroskedasticity") : 
		    I_("Allowing for groupwise heteroskedasticity"));
	} else {
	    pprintf(prn, (utf)? _("Weights based on per-unit error variances") : 
		    I_("Weights based on per-unit error variances"));
	}
	pputc(prn, '\n');
    }

    /* weight variable for WLS and ARCH */
    else if ((pmod->ci == WLS && !pmod->aux) || pmod->ci == ARCH) {
	if (tex) {
	    tex_escape(vname, pdinfo->varname[pmod->nwt]);
	    pputs(prn, "\\\\\n");
	}
	pprintf(prn, "%s: %s\n", 
		(utf)? _("Variable used as weight") : I_("Variable used as weight"), 
		(tex)? vname : pdinfo->varname[pmod->nwt]);
    }

    /* rhohat for CORC and HILU (TeX) */
    else if (pmod->ci == CORC || pmod->ci == HILU || pmod->ci == PWE) {
	if (tex) {
	    pprintf(prn, "$\\hat{\\rho}$ = %g\n", 
		    gretl_model_get_double(pmod, "rho_in"));
	}
    } 

    /* message about new variable created */
    if (PLAIN_FORMAT(prn->format) && gretl_msg[0] != '\0' &&
	strstr(gretl_msg, _("Replaced")) == NULL &&
	strstr(gretl_msg, _("Generated")) == NULL) {
	pprintf(prn, "%s\n", gretl_msg);
    }

    if (pmod->missmask == NULL && gretl_model_get_int(pmod, "wt_dummy")) { 
	/* FIXME alt formats */
	pprintf(prn, "%s %d\n", 
		(utf)? _("Weight var is a dummy variable, effective obs =") :
		I_("Weight var is a dummy variable, effective obs ="),
		pmod->nobs);
    } 

    if (RTF_FORMAT(prn->format)) {
	pputs(prn, "\\par\n");
    } else {
	pputc(prn, '\n');
    }
}

static void model_format_start (PRN *prn)
{
    if (TEX_FORMAT(prn->format)) {
	if (STANDALONE(prn->format)) {
	    gretl_tex_preamble(prn, 0);
	}
	pputs(prn, "\\begin{center}\n");
    }

    else if (RTF_FORMAT(prn->format)) {
	pputs(prn, "{\\rtf1\\par\n\\qc ");
	return;
    }
}

#define RTF_COEFF_ROW  "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                       "\\cellx500\\cellx1900\\cellx3300\\cellx4700\\cellx6100" \
                       "\\cellx7500\\cellx8000\n\\intbl"

#define RTF_DISCRETE_ROW  "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                       "\\cellx500\\cellx1900\\cellx3300\\cellx4700\\cellx6100" \
                       "\\cellx8000\n\\intbl"

#define RTF_ROOT_ROW   "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                       "\\cellx500\\cellx1500\\cellx2900\\cellx4300" \
                       "\\cellx5700\\cellx7100\n\\intbl"


static void print_coeff_table_start (const MODEL *pmod, PRN *prn, int discrete)
{
    if (PLAIN_FORMAT(prn->format)) {
	if (discrete) {
	    pputs(prn, _("      VARIABLE      COEFFICIENT        STDERROR"
			   "       T STAT       SLOPE\n"));
	    pprintf(prn, "                                                 "
		    "                 %s\n", _("(at mean)"));
	} else if (pmod->ci == NLS) {
	    pputs(prn, _("      PARAMETER      ESTIMATE          STDERROR"
			   "       T STAT   2Prob(t > |T|)\n\n"));
	} else {
	    pputs(prn, _("      VARIABLE      COEFFICIENT        STDERROR"
			   "       T STAT   2Prob(t > |T|)\n\n"));
	}
	return;
    } 
    else {
	char col1[16], col2[16];

	if (pmod->ci == NLS) {
	    strcpy(col1, N_("Parameter"));
	    strcpy(col2, N_("Estimate"));
	} else {
	    strcpy(col1, N_("Variable"));
	    strcpy(col2, N_("Coefficient"));
	}	    

	if (TEX_FORMAT(prn->format)) {
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
		    "        \\multicolumn{1}{c}{%s%s} \\\\[1ex]\n",
		    pt, pt, pt, pt, pt, pt, pt, pt, I_(col1),
		    I_(col2), I_("Std.\\ Error"), 
		    I_("$t$-statistic"), 
		    (discrete)? I_("Slope"): I_("p-value"),
		    (discrete)? "$^*$" : "");
	    return;
	}   

	if (RTF_FORMAT(prn->format)) {
	    if (discrete) {
		pprintf(prn, "{" RTF_DISCRETE_ROW
			" \\qr \\cell"
			" \\qc {\\i %s}\\cell"
			" \\qc {\\i %s}\\cell"
			" \\qc {\\i %s}\\cell"
			" \\qc {\\i %s}\\cell"
			" \\qc {\\i %s{\\super *}}\\cell"
			" \\intbl \\row\n",
			I_(col1), I_(col2), I_("Std. Error"), 
			I_("t-statistic"), I_("Slope"));
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
			I_(col1), I_(col2), I_("Std. Error"), 
			I_("t-statistic"), I_("p-value"));
	    }
	    return;
	} 
    }
}

static void print_coeff_table_end (PRN *prn)
{
    if (PLAIN_FORMAT(prn->format)) {
	pputc(prn, '\n');
    } else if (TEX_FORMAT(prn->format)) {
	pputs(prn, "\\end{tabular*}\n\n");
    } else if (RTF_FORMAT(prn->format)) {
	pputs(prn, "}\n\n");
    }
}

static void model_format_end (PRN *prn)
{
    if (TEX_FORMAT(prn->format)) {
	pputs(prn, "\n\\end{center}\n");
	if (STANDALONE(prn->format)) {
	    pputs(prn, "\n\\end{document}\n"); 
	}
    }

    else if (RTF_FORMAT(prn->format)) {
	pputs(prn, "\n}\n");
    }
} 

static int 
print_coefficients (const MODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    int i, err = 0, gotnan = 0;
    int n = pmod->ncoeff;
    int longnames = 0;
    int gn = -1;

    if (pmod->ci == GARCH) {
	gn = pmod->list[0] - 4;
    }

    if (PLAIN_FORMAT(prn->format) && 
	pmod->aux != AUX_ARCH &&
	pmod->ci != NLS && 
	pmod->ci != ARMA && 
	pmod->ci != GARCH) {
	for (i=0; i<n; i++) {
	    int len = strlen(pdinfo->varname[pmod->list[i+2]]);

	    if (len > 8) {
		longnames = 1;
	    }
	}
    }

    for (i=0; i<n; i++) {
	if (PLAIN_FORMAT(prn->format)) {
	    if (i == gn) {
		pputc(prn, '\n');
	    }
	    err = print_coeff(pdinfo, pmod, i + 2, longnames, prn);
	} else if (TEX_FORMAT(prn->format)) {
	    if (i == gn) {
		pputs(prn, "\\\\ \n");
	    }
	    err = tex_print_coeff(pdinfo, pmod, i + 2, prn);
	} else if (RTF_FORMAT(prn->format)) {
	    if (i == gn) {
		pputc(prn, '\n');
	    }
	    err = rtf_print_coeff(pdinfo, pmod, i + 2, prn);
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
		pt, pt);
    }
}

static void print_middle_table_end (PRN *prn)
{
    if (TEX_FORMAT(prn->format)) {
	pputs(prn, "\\end{tabular}\n\n");
    } else if (RTF_FORMAT(prn->format)) {
	pputs(prn, "\\par\n");
    } else {
	pputc(prn, '\n');
    }
}

static void r_squared_message (PRN *prn)
{
    pprintf(prn, "%s.\n\n",    
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

    else if (RTF_FORMAT(prn->format)) { /* FIXME */
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

static void print_ll (const MODEL *pmod, PRN *prn)
{
    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "  %s = %.3f\n", _("Log-likelihood"), pmod->lnL);
    } else if (RTF_FORMAT(prn->format)) {
	pprintf(prn, RTFTAB "%s = %.3f\n", I_("Log-likelihood"), pmod->lnL);
    } else if (TEX_FORMAT(prn->format)) {
	char xstr[32];

	tex_dcolumn_double(pmod->lnL, xstr);
	pprintf(prn, "%s & %s \\\\\n", I_("Log-likelihood"), xstr);
    }
}

/**
 * printmodel:
 * @pmod: pointer to gretl model.
 * @pdinfo: data information struct.
 * @opt: may contain OPT_O to print covariance matrix.
 * @prn: gretl printing struct.
 *
 * Print to @prn the estimates in @pmod plus associated statistics.
 * 
 * Returns: 0 on success, 1 if some of the values to print were %NaN.
 */

int printmodel (MODEL *pmod, const DATAINFO *pdinfo, gretlopt opt, 
		PRN *prn)
{
    int gotnan = 0;
    int is_discrete = (pmod->ci == PROBIT || pmod->ci == LOGIT);

    if (prn == NULL || (opt & OPT_Q)) {
	return 0;
    }

    if (prn->format != GRETL_PRINT_FORMAT_PLAIN) {
	model_format_start(prn);
    } else if (pmod->ci == ARMA || pmod->ci == GARCH || 
	       pmod->ci == TOBIT || pmod->ci == WLS) {
	int iters = gretl_model_get_int(pmod, "iters");

	if (iters > 0) {
	    pprintf(prn, _("Convergence achieved after %d iterations\n"), iters);
	}
    }

    print_model_heading(pmod, pdinfo, prn);

    print_coeff_table_start(pmod, prn, is_discrete);

    gotnan = print_coefficients(pmod, pdinfo, prn);

    if (pmod->ci == AR) {
	print_rho_terms(pmod, prn); 
    }

    print_coeff_table_end(prn);

    if (pmod->aux == AUX_ARCH || pmod->aux == AUX_ADF || 
	pmod->aux == AUX_RESET || pmod->aux == AUX_SCR ||
	pmod->aux == AUX_DF || pmod->aux == AUX_KPSS) {
	goto close_format;
    }

    if (is_discrete) {
	print_discrete_statistics(pmod, pdinfo, prn);
	/* AIC, BIC? */
	goto close_format;
    }

    if (pmod->ci == TOBIT) {
	print_middle_table_start(prn);
	depvarstats(pmod, prn);
	print_tobit_stats(pmod, prn);
	info_stats_lines(pmod, prn);
	print_middle_table_end(prn);
	goto close_format;
    }	

    if (pmod->ci == LAD) {
	print_middle_table_start(prn);
	depvarstats(pmod, prn);
	ladstats(pmod, prn);
	print_middle_table_end(prn);
	goto close_format;
    }

    if (pmod->ci == GARCH) {
	print_middle_table_start(prn);
	depvarstats(pmod, prn);
	print_arma_stats(pmod, prn);
	print_middle_table_end(prn);
	goto close_format;
    }    

    if (pmod->aux == AUX_SYS) {
	print_middle_table_start(prn);
	depvarstats(pmod, prn);
	essline(pmod, prn, 0);
	print_middle_table_end(prn);
	goto close_format;
    }    

    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG || 
	pmod->aux == AUX_AR) {
	print_middle_table_start(prn);
	rsqline(pmod, prn);
	print_middle_table_end(prn);
	goto close_format;
    }

    if (pmod->aux == AUX_COINT) {
	rsqline(pmod, prn);
	dwline(pmod, prn);
	if (!PLAIN_FORMAT(prn->format)) {
	    print_middle_table_end(prn);
	}
	goto close_format;
    }

    if (!pmod->ifc && pmod->ci != NLS && PLAIN_FORMAT(prn->format)) {
	noconst(pmod, prn);
    }
    
    if (pmod->aux == AUX_WHITE) { 
	rsqline(pmod, prn);
	print_whites_results(pmod, prn);
	goto close_format;
    }

    if (pmod->ci == LOGISTIC) {
	original_stats_message(prn);
    }

    if (pmod->ci == OLS || pmod->ci == VAR || pmod->ci == TSLS 
	|| pmod->ci == HCCM || pmod->ci == POOLED || pmod->ci == NLS
	|| (pmod->ci == AR && pmod->arinfo->arlist[0] == 1)
	|| pmod->ci == ARMA || pmod->ci == LOGISTIC || pmod->ci == TOBIT
	|| (pmod->ci == WLS && gretl_model_get_int(pmod, "wt_dummy"))) {
	print_middle_table_start(prn);
	if (pmod->ci != VAR) {
	    depvarstats(pmod, prn);
	}
	if (essline(pmod, prn, 0)) {
	    print_middle_table_end(prn);
	    goto close_format;
	}

	rsqline(pmod, prn);

	if (pmod->ci != NLS) {
	    Fline(pmod, prn);
	}

	if (dataset_is_time_series(pdinfo)) {
	    if (pmod->ci == OLS || pmod->ci == VAR ||
		(pmod->ci == WLS && gretl_model_get_int(pmod, "wt_dummy"))) {
		dwline(pmod, prn);
		if (pmod->ci != VAR && gretl_model_get_int(pmod, "ldepvar")) {
		    dhline(pmod, prn);
		} 
	    }
	    /* FIXME -- check output below */
	    if (pmod->ci == HCCM || pmod->ci == TSLS) {
		dwline(pmod, prn);
	    }
	}

	if (pmod->ci == ARMA) {
	    print_arma_stats(pmod, prn);
	} else {
	    info_stats_lines(pmod, prn);
	}

	print_middle_table_end(prn);

	if (pmod->ci == ARMA) {
	    print_arma_roots(pmod, prn);
	}

	if (pmod->ci == TSLS && PLAIN_FORMAT(prn->format)) {
	    r_squared_message(prn);
	}
    }

    else if (pmod->ci == WLS && gretl_model_get_int(pmod, "iters")) {
	/* panel ML estimation */
	print_middle_table_start(prn);
	depvarstats(pmod, prn);
	print_ll(pmod, prn);
	info_stats_lines(pmod, prn);
	print_middle_table_end(prn);
    }

    else if (pmod->ci == HSK || pmod->ci == ARCH ||
	     (pmod->ci == WLS && !gretl_model_get_int(pmod, "wt_dummy"))) {

	weighted_stats_message(prn);
	print_middle_table_start(prn);

	if (essline(pmod, prn, 1)) {
	    print_middle_table_end(prn);
	    goto close_format;
	}

	rsqline(pmod, prn);
	Fline(pmod, prn);

	if (dataset_is_time_series(pdinfo)) {
	    dwline(pmod, prn);
	}

	info_stats_lines(pmod, prn);

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

    else if (pmod->ci == AR || pmod->ci == CORC || 
	     pmod->ci == HILU || pmod->ci == PWE) {

	rho_differenced_stats_message(prn);

	print_middle_table_start(prn);
	if (essline(pmod, prn, 0)) {
	    print_middle_table_end(prn);
	    goto close_format;
	}
	rsqline(pmod, prn);
	Fline(pmod, prn);
	dwline(pmod, prn);
	info_stats_lines(pmod, prn);
	print_middle_table_end(prn);
    }

    if (PLAIN_FORMAT(prn->format) && 
	pmod->ci != ARMA && pmod->ci != NLS) {
	pval_max_line(pmod, pdinfo, prn);
    }

 close_format:

    if (pmod->ci != LAD && (opt & OPT_O)) {
	outcovmx(pmod, pdinfo, prn);
    }

    if (pmod->ntests > 0) {
	print_model_tests(pmod, prn);
    }

#if 0
    if (PLAIN_FORMAT(prn->format) && pmod->aux == AUX_ADF) {
	info_stats_lines(pmod, prn);
    }
#endif

    if (!PLAIN_FORMAT(prn->format)) {
	model_format_end(prn);
    }
    
    if (gotnan) {
	pmod->errcode = E_NAN;
    }

    return gotnan;
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
			int c, int longnames, PRN *prn)
{
    double t, pvalue = 999.0;
    int gotnan = 0;
    int do_pval = (pmod->ci != LOGIT && pmod->ci != PROBIT);
    char varname[24];

    /* special treatment for ARCH model coefficients, NLS, ARMA, VAR */
    if (pmod->aux == AUX_ARCH) {
	make_cname(pdinfo->varname[pmod->list[c]], varname);
    } else if (pmod->ci == NLS || pmod->ci == ARMA || pmod->ci == GARCH) {
	strcpy(varname, pmod->params[c-1]);
    } else {
	strcpy(varname, pdinfo->varname[pmod->list[c]]);
    }

    if (longnames) {
	pprintf(prn, " %13s ", varname);
    } else if (pmod->ci == GARCH) {
	pprintf(prn, "      %8s ", varname);
    } else if (pmod->ci == ARMA) {
	pprintf(prn, " %3d) %8s ", c - 1, varname);
    } else {
	pprintf(prn, " %3d) %8s ", pmod->list[c], varname);
    }

    bufspace((strlen(varname) > 12)? 1 : 2, prn);

    /* print coeff value if well-defined */
    if (isnan(pmod->coeff[c-2]) || na(pmod->coeff[c-2])) {
	pprintf(prn, "%*s", UTF_WIDTH(_("undefined"), 17), _("undefined"));
	gotnan = 1;
    } else {
	gretl_print_value(pmod->coeff[c-2], prn);
    }

    bufspace(2, prn);

    /* get out if std error is undefined */
    if (isnan(pmod->sderr[c-2]) || na(pmod->sderr[c-2])) {
	pprintf(prn, "%*s\n", UTF_WIDTH(_("undefined"), 16), _("undefined"));
	return 1;
    }

    gretl_print_value (pmod->sderr[c-2], prn); 

    /* std error is well-defined, but is it positive? */
    if (pmod->sderr[c-2] > 0.) {
	t = pmod->coeff[c-2] / pmod->sderr[c-2];
	if (fabs(t) >= 1000.0) { /* || t < .001 ? */
	    char numstr[9];

	    sprintf(numstr, "%#8.2G", t);
	    pprintf(prn, " %8s", numstr);
	} else {
	    pprintf(prn, " %7.3f", t);
	}
	if (pmod->aux == AUX_ADF || pmod->aux == AUX_DF) {
	    if (c == gretl_model_get_int(pmod, "dfnum")) {
		char pvalstr[16];

		pvalue = gretl_model_get_double(pmod, "dfpval");
		print_pval_str(pvalue, pvalstr);
		pprintf(prn, "%*s", UTF_WIDTH(pvalstr, 12), pvalstr);
	    } 
	    do_pval = 0;
	}
	if (do_pval) {
	    char pvalstr[16];

	    pvalue = tprob(t, pmod->dfd);
	    print_pval_str(pvalue, pvalstr);
	    pprintf(prn, "%*s", UTF_WIDTH(pvalstr, 12), pvalstr);
	}
    } else if (do_pval) { /* zero standard error */
	do_pval = 0;
	pprintf(prn, "     %*s", UTF_WIDTH(_("undefined"), 12), _("undefined"));
    }

    if (do_pval) {
	if (pvalue < 0.01) pputs(prn, " ***");
	else if (pvalue < 0.05) pputs(prn, " **");
	else if (pvalue < 0.10) pputs(prn, " *");
    } 
    else if (pmod->list[c] != 0 && 
	     (pmod->ci == LOGIT || pmod->ci == PROBIT)) { 
	double *slopes = gretl_model_get_data(pmod, "slopes");

	if (slopes != NULL) {
	    gretl_print_value(slopes[c-2], prn);
	}
    }

    pputc(prn, '\n');

    return gotnan;
}

/* ......................................................... */ 

static void rtf_print_double (double xx, PRN *prn)
{
    char numstr[32];

    xx = screen_zero(xx);
    sprintf(numstr, "%.*g", GRETL_DIGITS, xx);
    gretl_fix_exponent(numstr);
    pprintf(prn, " \\qc %s\\cell", numstr);
}

/* ......................................................... */ 

static int rtf_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			    int c, PRN *prn)
{
    double t, pvalue = 999.0;
    int gotnan = 0;
    int do_pval = (pmod->ci != LOGIT && pmod->ci != PROBIT);
    char varname[12];


    /* special treatment for ARCH model coefficients, NLS, ARMA, VAR */
    if (pmod->aux == AUX_ARCH) {
	make_cname(pdinfo->varname[pmod->list[c]], varname);
    } else if (pmod->ci == NLS || pmod->ci == ARMA || pmod->ci == GARCH) {
	strcpy(varname, pmod->params[c-1]);
    } else {
	strcpy(varname, pdinfo->varname[pmod->list[c]]);
    }

    pputs(prn, RTF_COEFF_ROW);

    pprintf(prn, " \\qr %d\\cell \\ql %s\\cell", 
	    ((pmod->ci == ARMA || pmod->ci == GARCH)? c - 1 : pmod->list[c]), 
	    varname);

    if (isnan(pmod->coeff[c-2]) || na(pmod->coeff[c-2])) {
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	gotnan = 1;
    } else {
	rtf_print_double(pmod->coeff[c-2], prn);
    }

    if (isnan(pmod->sderr[c-2]) || na(pmod->sderr[c-2])) {
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	goto rtf_coeff_getout;
    } 

    rtf_print_double(pmod->sderr[c-2], prn); 

    if (pmod->sderr[c-2] > 0.) {
	t = pmod->coeff[c-2] / pmod->sderr[c-2];
	pprintf(prn, " \\qc %.4f\\cell", t);
	if (pmod->aux == AUX_ADF || pmod->aux == AUX_DF) {
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
	    pputs(prn, " \\ql ***\\cell");
	else if (pvalue < 0.05) 
	    pputs(prn, " \\ql **\\cell");
	else if (pvalue < 0.10) 
	    pputs(prn, " \\ql *\\cell");
	else 
	    pputs(prn, " \\ql \\cell");
    } else if (pmod->list[c] != 0 && 
	       (pmod->ci == LOGIT || pmod->ci == PROBIT)) { 
	double *slopes = gretl_model_get_data(pmod, "slopes");

	if (slopes != NULL) {
	    rtf_print_double(slopes[c-2], prn);
	}	
    }

 rtf_coeff_getout:
    pputs(prn, " \\intbl \\row\n");

    return gotnan;
}

/* ......................................................... */ 

static void print_rho (const ARINFO *arinfo, int c, int dfd, PRN *prn)
{
    char ustr[32];
    double xx = arinfo->rho[c] / arinfo->sderr[c];

    if (PLAIN_FORMAT(prn->format)) {
	sprintf(ustr, "u_%d", arinfo->arlist[c]);
	pprintf(prn, "%14s", ustr); 
	bufspace(3, prn);
	gretl_print_value (arinfo->rho[c], prn);
	bufspace(2, prn);
	gretl_print_value (arinfo->sderr[c], prn); 
	pprintf(prn, " %7.3f ", xx, tprob(xx, dfd));
	if (1) {
	    char pvalstr[16];
	    double pval;

	    pval = tprob(xx, dfd);
	    print_pval_str(pval, pvalstr);
	    pprintf(prn, "%*s\n", UTF_WIDTH(pvalstr, 12), pvalstr);
	}
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

    if (pmod->arinfo->arlist[0] > 1) {
	dfd = pmod->dfd + (pmod->ncoeff - pmod->arinfo->arlist[0]);
    } else {
	dfd = pmod->dfd;
    }

    for (i=1; i<=pmod->arinfo->arlist[0]; i++) {
	print_rho(pmod->arinfo, i, dfd, prn);
	xx += pmod->arinfo->rho[i]; 
    }

    if (pmod->arinfo->arlist[0] > 1 && PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "\n%s = %#g\n\n", _("Sum of AR coefficients"), xx);
    }
}

static void tex_float_str (double x, char *str)
{
    if (x < 0.) {
	strcpy(str, "$-$");
	sprintf(str + 3, "%.*g", GRETL_DIGITS, -x);
    } else {
	sprintf(str, "%.*g", GRETL_DIGITS, x);
    }
}

const char *roots_hdr = N_("                        Real  Imaginary"
			   "    Modulus  Frequency");
const char *root_fmt = "%8s%3d%17.4f%11.4f%11.4f%11.4f\n";
const char *roots_sep = "  -----------------------------------------"
                        "------------------";

static void print_arma_roots (const MODEL *pmod, PRN *prn)
{
    cmplx *roots = gretl_model_get_data(pmod, "roots");

    if (roots != NULL) {
	int p = pmod->list[1];
	int q = pmod->list[2];
	int i;
	double mod, fr;

	if (PLAIN_FORMAT(prn->format)) {
	    pprintf(prn, "\n%s\n%s\n", _(roots_hdr), roots_sep);
	} else if (TEX_FORMAT(prn->format)) {
	    pputs(prn, "\n\\vspace{1em}\n\n");
	    pputs(prn, "\\begin{tabular}{llrrrrr}\n");
	    pprintf(prn, "& & & %s & %s & %s & %s \\\\ \\hline\n", 
		    I_("Real"), I_("Imaginary"), I_("Modulus"), I_("Frequency"));
	} else if (RTF_FORMAT(prn->format)) {
	    pputs(prn, "\n\\par\n{" RTF_ROOT_ROW);
	    pprintf(prn, "\\qr \\cell \\qc \\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell \\intbl \\row\n",
		    I_("Real"), I_("Imaginary"), I_("Modulus"), I_("Frequency"));
	}

	if (p > 0) {
	    if (PLAIN_FORMAT(prn->format)) {
		pprintf(prn, "  %s\n", _("AR"));
	    } else if (TEX_FORMAT(prn->format)) {
		pprintf(prn, "%s \\\\ \n", I_("AR"));
	    } else if (RTF_FORMAT(prn->format)) {
		pputs(prn, RTF_ROOT_ROW);
		pprintf(prn, "\\ql %s\\cell\\ql \\cell\\ql \\cell\\ql \\cell\\ql \\cell"
			"\\ql\\cell \\intbl \\row\n", I_("AR"));
	    }
	    for (i=0; i<p; i++) {
		if (roots[i].i != 0) {
		    mod = roots[i].r * roots[i].r + roots[i].i * roots[i].i;
		    mod = sqrt(mod);
		} else {
		    mod = fabs(roots[i].r);
		}
		fr = atan2(roots[i].i, roots[i].r) / (2.0 * M_PI);
		if (PLAIN_FORMAT(prn->format)) {
		    pprintf(prn, root_fmt, _("Root"), i + 1, 
			    roots[i].r, roots[i].i, mod, fr);
		} else if (TEX_FORMAT(prn->format)) {
		    pprintf(prn, "& %s & %d & $%.4f$ & $%.4f$ & $%.4f$ & $%.4f$ \\\\ ",
			    I_("Root"), i + 1, roots[i].r, roots[i].i, mod, fr);
		    if (i == p - 1 && q == 0) pputs(prn, "\\hline\n");
		    else pputc(prn, '\n');
		} else if (RTF_FORMAT(prn->format)) {
		    pputs(prn, RTF_ROOT_ROW);
		    pprintf(prn, "\\ql \\cell \\ql %s %d \\cell"
			    " \\qr %.4f\\cell"
			    " \\qr %.4f\\cell"
			    " \\qr %.4f\\cell"
			    " \\qr %.4f\\cell \\intbl \\row\n",
			    I_("Root"), i + 1, roots[i].r, roots[i].i, mod, fr);
		}
	    }
	} 

	if (q > 0) {
	    if (PLAIN_FORMAT(prn->format)) {
		pprintf(prn, "  %s\n", _("MA"));
	    } else if (TEX_FORMAT(prn->format)) {
		pprintf(prn, "%s \\\\ \n", I_("MA"));
	    } else if (RTF_FORMAT(prn->format)) {
		pputs(prn, RTF_ROOT_ROW);
		pprintf(prn, "\\ql %s\\cell\\ql \\cell\\ql \\cell\\ql \\cell\\ql \\cell"
			"\\ql\\cell \\intbl \\row\n", I_("MA"));
	    }
	    for (i=p; i<p+q; i++) {
		if (roots[i].i != 0) {
		    mod = roots[i].r * roots[i].r + roots[i].i * roots[i].i;
		    mod = sqrt(mod);
		} else {
		    mod = fabs(roots[i].r);
		}
		fr = atan2(roots[i].i, roots[i].r) / (2.0 * M_PI);
		if (PLAIN_FORMAT(prn->format)) {
		    pprintf(prn, root_fmt, _("Root"), i - p + 1, 
			    roots[i].r, roots[i].i, mod, fr);
		} else if (TEX_FORMAT(prn->format)) {
		    pprintf(prn, "& %s & %d & $%.4f$ & $%.4f$ & $%.4f$ & $%.4f$ \\\\ ",
			    I_("Root"), i - p + 1, roots[i].r, roots[i].i, mod, fr);
		    if (i == p + q - 1) pputs(prn, "\\hline\n");
		    else pputc(prn, '\n');
		} else if (RTF_FORMAT(prn->format)) {
		    pputs(prn, RTF_ROOT_ROW);
		    pprintf(prn, "\\ql \\cell \\ql %s %d \\cell"
			    " \\qr %.4f\\cell"
			    " \\qr %.4f\\cell"
			    " \\qr %.4f\\cell"
			    " \\qr %.4f\\cell \\intbl \\row\n",
			    I_("Root"), i - p + 1, roots[i].r, roots[i].i, mod, fr);
		}
	    }
	}

	if (PLAIN_FORMAT(prn->format)) {
	    pprintf(prn, "%s\n\n", roots_sep);
	} else if (TEX_FORMAT(prn->format)) {
	    pputs(prn, "\\end{tabular}\n");
	} else if (RTF_FORMAT(prn->format)) {
	    pputs(prn, "}\n");
	}
    }
}

static void print_tobit_stats (const MODEL *pmod, PRN *prn)
{
    int cenc = gretl_model_get_int(pmod, "censobs");
    double cenpc = 100.0 * cenc / pmod->nobs;

    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "  %s: %d (%.1f%%)\n", _("Censored observations"), cenc, cenpc);
	pprintf(prn, "  %s = %.*g\n", _("sigma"), GRETL_DIGITS, pmod->sigma);
	pprintf(prn, "  %s = %.3f\n", _("Log-likelihood"), pmod->lnL);
    }
    else if (RTF_FORMAT(prn->format)) {
	pprintf(prn, RTFTAB "%s: %d (%.1f%%)\n", I_("Censored observations"), cenc, cenpc);
	pprintf(prn, RTFTAB "%s = %g\n", I_("sigma"), pmod->sigma);
	pprintf(prn, RTFTAB "%s = %.3f\n", I_("Log-likelihood"), pmod->lnL);
    }
    else if (TEX_FORMAT(prn->format)) {
	char xstr[32];

	pprintf(prn, "%s & \\multicolumn{1}{r}{%.1f\\%%} \\\\\n", 
		I_("Censored observations"), cenpc);
	tex_dcolumn_double(pmod->sigma, xstr);
	pprintf(prn, "$\\hat{\\sigma}$ & %s \\\\\n", xstr);
	tex_dcolumn_double(pmod->lnL, xstr);
	pprintf(prn, "%s & %s \\\\\n", I_("Log-likelihood"), xstr);
    }
}

static void print_arma_stats (const MODEL *pmod, PRN *prn)
{
    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "  %s = %.3f\n", _("Log-likelihood"), pmod->lnL);
	pprintf(prn, "  %s = %.3f\n", _("AIC"), pmod->criterion[C_AIC]);
	pprintf(prn, "  %s = %.3f\n", _("BIC"), pmod->criterion[C_BIC]);
    }
    else if (RTF_FORMAT(prn->format)) {
	pprintf(prn, RTFTAB "%s = %.3f\n", I_("Log-likelihood"), pmod->lnL);
	pprintf(prn, RTFTAB "%s = %.3f\n", I_("AIC"), pmod->criterion[C_AIC]);
	pprintf(prn, RTFTAB "%s = %.3f\n", I_("BIC"), pmod->criterion[C_BIC]);
    }
    else if (TEX_FORMAT(prn->format)) {
	char xstr[32];

	tex_dcolumn_double(pmod->lnL, xstr);
	pprintf(prn, "%s & %s \\\\\n", I_("Log-likelihood"), xstr);
	tex_dcolumn_double(pmod->criterion[C_AIC], xstr);
	pprintf(prn, "%s & %s \\\\\n", I_("AIC"), xstr);
	tex_dcolumn_double(pmod->criterion[C_BIC], xstr);
	pprintf(prn, "%s & %s \\\\\n", I_("BIC"), xstr);
    }
}

static void print_discrete_statistics (const MODEL *pmod, 
				       const DATAINFO *pdinfo,
				       PRN *prn)
{
    int i;
    int correct = gretl_model_get_int(pmod, "correct");
    double pc_correct = 100 * (double) correct / pmod->nobs;
    double model_chisq = gretl_model_get_double(pmod, "chisq");

    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, "  %s %s = %.3f\n", _("Mean of"), 
		pdinfo->varname[pmod->list[1]], pmod->ybar);
	pprintf(prn, "  %s = %d (%.1f%%)\n", _("Number of cases 'correctly predicted'"), 
		correct, pc_correct);
	pprintf(prn, "  f(beta'x) %s = %.3f\n", _("at mean of independent vars"), 
		pmod->sdy);
	pprintf(prn, "  %s = %.3f\n", _("Log-likelihood"), pmod->lnL);
	if (pmod->aux != AUX_OMIT && pmod->aux != AUX_ADD) {
	    i = pmod->ncoeff - 1;
	    pprintf(prn, "  %s: %s(%d) = %.3f (%s %f)\n",
		    _("Likelihood ratio test"), _("Chi-square"), 
		    i, model_chisq, _("p-value"), chisq(model_chisq, i));
	}
	pputc(prn, '\n');
    }

    else if (RTF_FORMAT(prn->format)) {
	pputc(prn, '\n');
	pprintf(prn, "\\par {\\super *}%s\n", I_("Evaluated at the mean"));
	pprintf(prn, "\\par %s %s = %.3f\n", I_("Mean of"), 
		pdinfo->varname[pmod->list[1]], pmod->ybar);
	pprintf(prn, "\\par %s = %d (%.1f%%)\n", I_("Number of cases 'correctly predicted'"), 
		correct, pc_correct);
	pprintf(prn, "\\par f(beta'x) %s = %.3f\n", I_("at mean of independent vars"), 
		pmod->sdy);
	pprintf(prn, "\\par %s = %.3f\n", I_("Log-likelihood"), pmod->lnL);
	if (pmod->aux != AUX_OMIT && pmod->aux != AUX_ADD) {
	    i = pmod->ncoeff - 1;
	    pprintf(prn, "\\par %s: %s(%d) = %.3f (%s %f)\\par\n",
		    I_("Likelihood ratio test"), I_("Chi-square"), 
		    i, model_chisq, I_("p-value"), chisq(model_chisq, i));
	}
	pputc(prn, '\n');
    }

    else if (TEX_FORMAT(prn->format)) {
	char lnlstr[16];

	pprintf(prn, "\\begin{center}\n$^*$%s\n\\end{center}\n", 
		I_("Evaluated at the mean"));

	pputs(prn, "\\vspace{1em}\n\\begin{raggedright}\n");
	pprintf(prn, "%s %s = %.3f\\\\\n", I_("Mean of"), 
		pdinfo->varname[pmod->list[1]], pmod->ybar);
	pprintf(prn, "%s = %d (%.1f percent)\\\\\n", 
		I_("Number of cases `correctly predicted'"), 
		correct, pc_correct);
	pprintf(prn, "$f(\\beta'x)$ %s = %.3f\\\\\n", I_("at mean of independent vars"), 
		pmod->sdy);
	tex_float_str(pmod->lnL, lnlstr);
	pprintf(prn, "%s = %s\\\\\n", I_("Log-likelihood"), lnlstr);
	if (pmod->aux != AUX_OMIT && pmod->aux != AUX_ADD) {
	    i = pmod->ncoeff - 1;
	    pprintf(prn, "%s: $\\chi^2_{%d}$ = %.3f (%s %f)\\\\\n",
		    I_("Likelihood ratio test"), 
		    i, model_chisq, I_("p-value"), chisq(model_chisq, i));
	}
	pputs(prn, "\\end{raggedright}\n");
    }
}

/* ....................................................... */

static void mp_other_stats (const mp_results *mpvals, PRN *prn)
{
    char fstr[16];
    int len = 24;

    if (doing_nls()) len = 36;
    
    pprintf(prn, "%-*s", len, _("Standard error"));
    gretl_print_fullwidth_double(mpvals->sigma, GRETL_MP_DIGITS, prn);
    pputc(prn, '\n');

    pprintf(prn, "%-*s", len, _("Error Sum of Squares"));
    gretl_print_fullwidth_double(mpvals->ess, GRETL_MP_DIGITS, prn);
    pputc(prn, '\n');

    pprintf(prn, "%-*s", len, _("Unadjusted R-squared"));
    gretl_print_fullwidth_double(mpvals->rsq, GRETL_MP_DIGITS, prn);
    pputc(prn, '\n');

    pprintf(prn, "%-*s", len, _("Adjusted R-squared"));
    gretl_print_fullwidth_double(mpvals->adjrsq, GRETL_MP_DIGITS, prn);
    pputc(prn, '\n');

    sprintf(fstr, "F(%d, %d)", mpvals->dfn, mpvals->dfd);
    pprintf(prn, "%-*s", len, fstr);

    if (na(mpvals->fstt)) {
	pprintf(prn, "            %s", _("undefined"));
    } else {
	gretl_print_fullwidth_double(mpvals->fstt, GRETL_MP_DIGITS, prn);
    }

    pputs(prn, "\n\n");
}

/* ....................................................... */

static void print_mpvals_coeff (const mp_results *mpvals, 
				int c, PRN *prn)
{
    pprintf(prn, " %3d) %8s ", mpvals->varlist[c+2], mpvals->varnames[c+1]);

    gretl_print_fullwidth_double(mpvals->coeff[c], 
				 GRETL_MP_DIGITS, prn);
    gretl_print_fullwidth_double(mpvals->sderr[c], 
				 GRETL_MP_DIGITS, prn); 

    pputc(prn, '\n');
}

/* ....................................................... */

void print_mpols_results (const mp_results *mpvals, DATAINFO *pdinfo,
			  PRN *prn)
{
    int i;
    char startdate[OBSLEN], enddate[OBSLEN];
    int nobs = mpvals->t2 - mpvals->t1 + 1;

    ntodate(startdate, mpvals->t1, pdinfo);
    ntodate(enddate, mpvals->t2, pdinfo);

    pputc(prn, '\n');

    if (!PLAIN_FORMAT(prn->format)) {
	pputs(prn, "FIXME: this is still to be implemented!\n\n");
    }

    if (PLAIN_FORMAT(prn->format)) {
	pprintf(prn, _("Multiple-precision OLS estimates using "
		       "the %d observations %s-%s\n"),
		nobs, startdate, enddate);
	pprintf(prn, "%s: %s\n\n", _("Dependent variable"),
		mpvals->varnames[0]);

	pputs(prn, _("      VARIABLE         COEFFICIENT          "
		       "        STD. ERROR\n"));
    }

    for (i=0; i<mpvals->ncoeff; i++) {
	print_mpvals_coeff(mpvals, i, prn);
    }
    pputc(prn, '\n');

    mp_other_stats (mpvals, prn);
}
