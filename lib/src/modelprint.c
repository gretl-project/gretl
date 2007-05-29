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

/* modelprint.c */ 

#include "libgretl.h"
#include "libset.h"
#include "system.h"
#include "texprint.h"

#define NO_RBAR_SQ(a) (a == AUX_SQ || a == AUX_LOG || a == AUX_WHITE || a == AUX_AR)

static int 
print_coefficients (const MODEL *pmod, const DATAINFO *pdinfo, PRN *prn);
static void depvarstats (const MODEL *pmod, PRN *prn);
static void print_rho_terms (const MODEL *pmod, PRN *prn);
static void print_binary_statistics (const MODEL *pmod, 
				     const DATAINFO *pdinfo,
				     PRN *prn);
static void print_arma_stats (const MODEL *pmod, PRN *prn);
static void print_arma_roots (const MODEL *pmod, PRN *prn);
static void print_tobit_stats (const MODEL *pmod, PRN *prn);
static void print_heckit_stats (const MODEL *pmod, PRN *prn);
static void print_poisson_offset (const MODEL *pmod, const DATAINFO *pdinfo, 
				  PRN *prn);
static void print_ll (const MODEL *pmod, PRN *prn);

#define RTFTAB "\\par \\ql \\tab "

#define XDIGITS(m) (((m)->ci == MPOLS)? GRETL_MP_DIGITS : GRETL_DIGITS)

#define ordered_model(m) ((m->ci == LOGIT || m->ci == PROBIT) && \
                           gretl_model_get_int(m, "ordered"))

#define binary_model(m) ((m->ci == LOGIT || m->ci == PROBIT) && \
                         !gretl_model_get_int(m, "ordered"))

#define pooled_model(m) (m->ci == OLS && \
                         gretl_model_get_int(m, "pooled")) 

void model_coeff_init (model_coeff *mc)
{
    mc->b = NADBL;
    mc->se = NADBL;
    mc->tval = NADBL;
    mc->pval = NADBL;
    mc->slope = NADBL;  
    mc->show_pval = 1;
    mc->df_pval = 0;
    mc->name[0] = '\0';
}

static void print_y_median (const MODEL *pmod, PRN *prn)
{
    double m = gretl_model_get_double(pmod, "ymedian");

    if (na(m)) {
	return;
    }

    if (plain_format(prn)) {
	pprintf(prn, "  %s = %.*g\n", _("Median of dependent variable"), 
		XDIGITS(pmod), m);
    } else if (tex_format(prn)) {
	char s[32];

	tex_dcolumn_double(m, s);
	pprintf(prn, "%s & %s \\\\\n", I_("Median of dependent variable"), s);
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_("Median of dependent variable"), 
		m);
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", I_("Median of dependent variable"), 
		prn_delim(prn), m);
    }	
}

static void depvarstats (const MODEL *pmod, PRN *prn)
{
    if (pmod->ci == LAD) {
	print_y_median(pmod, prn);
    }

    if (na(pmod->ybar) || na(pmod->sdy)) {
	return;
    }

    if (plain_format(prn)) {
	pprintf(prn, "  %s = %.*g\n", _("Mean of dependent variable"), 
		XDIGITS(pmod), pmod->ybar);
	pprintf(prn, "  %s = %.*g\n", _("Standard deviation of dep. var."), 
		XDIGITS(pmod), pmod->sdy);
    } else if (tex_format(prn)) {
	char x1str[32], x2str[32];

	tex_dcolumn_double(pmod->ybar, x1str);
	tex_dcolumn_double(pmod->sdy, x2str);
	pprintf(prn, "%s & %s \\\\\n %s & %s \\\\\n",
		I_("Mean of dependent variable"), x1str,
		I_("S.D. of dependent variable"), x2str);
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_("Mean of dependent variable"), 
		pmod->ybar);
	pprintf(prn, RTFTAB "%s = %g\n", I_("Standard deviation of dep. var."), 
		pmod->sdy);
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", I_("Mean of dependent variable"), 
		prn_delim(prn), pmod->ybar);
	pprintf(prn, "\"%s\"%c%.15g\n", I_("Standard deviation of dep. var."), 
		prn_delim(prn), pmod->sdy);
    }	
}

static void garch_variance_line (const MODEL *pmod, PRN *prn)
{
    const char *varstr = N_("Unconditional error variance");
    double v = pmod->sigma * pmod->sigma;

    if (plain_format(prn)) {  
	pprintf(prn, "  %s = %.*g\n", _(varstr), GRETL_DIGITS, v);
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_(varstr), v);
    } else if (tex_format(prn)) {
	char xstr[32];

	tex_dcolumn_double(v, xstr);
	pprintf(prn, "%s & %s \\\\\n", I_(varstr), xstr);
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", I_(varstr), prn_delim(prn), v);
    }
}

static int 
real_essline (const MODEL *pmod, double ess, double sigma, PRN *prn)
{
    if (ess < 0) {
	if (plain_format(prn)) {
	    pprintf(prn, _("Error sum of squares (%g) is not > 0"), ess);
	    pputs(prn, "\n\n");
	}
	return 1;
    }

    if (plain_format(prn)) {    
	pprintf(prn, "  %s = %.*g\n", _("Sum of squared residuals"), 
		XDIGITS(pmod), ess);
	pprintf(prn, "  %s = %.*g\n", _("Standard error of residuals"), 
		XDIGITS(pmod), sigma);
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_("Sum of squared residuals"), 
		ess);
	pprintf(prn, RTFTAB "%s = %g\n", I_("Standard error of residuals"), 
		sigma);
    } else if (tex_format(prn)) {
	char x1str[32], x2str[32];

	tex_dcolumn_double(ess, x1str);
	tex_dcolumn_double(sigma, x2str);
	pprintf(prn, "%s & %s \\\\\n%s ($\\hat{\\sigma}$) & %s \\\\\n",
		I_("Sum of squared residuals"), x1str,
		I_("Standard error of residuals"), x2str);
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", I_("Sum of squared residuals"), 
		prn_delim(prn), ess);
	pprintf(prn, "\"%s\"%c%.15g\n", I_("Standard error of residuals"), 
		prn_delim(prn), sigma);
    }	

    return 0;
}

static int GMM_crit_line (const MODEL *pmod, PRN *prn)
{
    if (plain_format(prn)) {    
	pprintf(prn, "  %s = %.*g\n", _("GMM criterion"), 
		XDIGITS(pmod), pmod->ess);
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_("GMM criterion"), 
		pmod->ess);
    } else if (tex_format(prn)) {
	char xstr[32];

	tex_dcolumn_double(pmod->ess, xstr);
	pprintf(prn, "%s & %s \\\\\n",
		I_("GMM criterion"), xstr);
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", I_("GMM criterion"), 
		prn_delim(prn), pmod->ess);
    }	

    return 0;
}

static int essline (const MODEL *pmod, PRN *prn)
{
    if (pmod->ci == GMM) {
	return GMM_crit_line(pmod, prn);
    } else {
	return real_essline(pmod, pmod->ess, pmod->sigma, prn);
    }
}

static int essline_original (const MODEL *pmod, PRN *prn)
{
    double ess = gretl_model_get_double(pmod, "ess_orig");
    double sigma = gretl_model_get_double(pmod, "sigma_orig");

    if (na(ess) || na(sigma)) {
	return 1;
    }
    
    return real_essline(pmod, ess, sigma, prn);
}

static void rsqline (const MODEL *pmod, PRN *prn)
{
    const char *plainrsq[] = {
	N_("Unadjusted R-squared"),
	N_("Uncentered R-squared")
    };
    const char *texrsq[] = {
	N_("Unadjusted $R^2$"),
	N_("Uncentered $R^2$")
    };
    const char *rtfrsq[] = {
	N_("Unadjusted R{\\super 2}"),
	N_("Uncentered R{\\super 2}")
    };
    int ridx = 0;
    int adjr2 = 1;

    if (na(pmod->rsq)) {
	return;
    }

    if (!plain_format(prn) && pmod->ci == TSLS) {
	return;
    }

    if (NO_RBAR_SQ(pmod->aux) || na(pmod->adjrsq)) {
	adjr2 = 0;
    }

    if (gretl_model_get_int(pmod, "uncentered")) {
	ridx = 1;
    }

    if (plain_format(prn)) { 
	pprintf(prn, "  %s = %.*g\n", _(plainrsq[ridx]), 
		XDIGITS(pmod), pmod->rsq);
	if (adjr2) {
	    pprintf(prn, "  %s = %.*g\n", _("Adjusted R-squared"),  
		    XDIGITS(pmod), pmod->adjrsq);
	}
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_(rtfrsq[ridx]), pmod->rsq);
	if (adjr2) {
	    pprintf(prn, RTFTAB "%s = %g\n", I_("Adjusted R{\\super 2}"),  
		    pmod->adjrsq);
	}	
    } else if (tex_format(prn)) {  
	char r2[32];

	tex_dcolumn_double(pmod->rsq, r2);
	pprintf(prn, "%s & %s \\\\\n", I_(texrsq[ridx]), r2);
	if (adjr2) {
	    tex_dcolumn_double(pmod->adjrsq, r2);
	    pprintf(prn, "%s & %s \\\\\n", I_("Adjusted $\\bar{R}^2$"), r2);
	}
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", I_(plainrsq[ridx]), 
		prn_delim(prn), pmod->rsq);
	if (adjr2) {
	    pprintf(prn, "\"%s\"%c%.15g\n", I_("Adjusted R-squared"),  
		    prn_delim(prn), pmod->adjrsq);
	}
    }	
}

static void pseudorsqline (const MODEL *pmod, PRN *prn)
{
    if (na(pmod->rsq)) {
	return;
    }

    if (plain_format(prn)) { 
	pprintf(prn, "  %s = %.*g\n", _("McFadden's pseudo-R-squared"), 
		XDIGITS(pmod), pmod->rsq);
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_("McFadden's pseudo-R{\\super 2}"), 
		pmod->rsq);
    } else if (tex_format(prn)) {  
	char r2[32];

	tex_dcolumn_double(pmod->rsq, r2);
	pprintf(prn, "%s = %s \\\\\n", I_("McFadden's pseudo-$R^2$"), r2);
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", I_("McFadden's pseudo-R-squared"), 
		prn_delim(prn), pmod->rsq);
    }	
}

static const char *aic_str = N_("Akaike information criterion");
static const char *bic_str = N_("Schwarz Bayesian criterion");
static const char *hqc_str = N_("Hannan-Quinn criterion");
static const char *tex_hqc_str = N_("Hannan--Quinn criterion");

static const char *aic_abbrev = N_("AIC");
static const char *bic_abbrev = N_("BIC");
static const char *hqc_abbrev = N_("HQC");

static void info_stats_lines (const MODEL *pmod, PRN *prn)
{
    const double *crit = pmod->criterion;

    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG ||
	pmod->aux == AUX_WHITE || pmod->aux == AUX_AR) {
	return;
    }

    if (na(crit[C_AIC]) || na(crit[C_BIC]) || na(crit[C_HQC])) {
	return;
    }

    if (plain_format(prn)) { 
	pprintf(prn, "  %s (%s) = %.*g\n", _(aic_str), _(aic_abbrev),
		XDIGITS(pmod), crit[C_AIC]);
	pprintf(prn, "  %s (%s) = %.*g\n", _(bic_str), _(bic_abbrev),
		XDIGITS(pmod), crit[C_BIC]);
	pprintf(prn, "  %s (%s) = %.*g\n", _(hqc_str), _(hqc_abbrev),
		XDIGITS(pmod), crit[C_HQC]);
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_(aic_str), crit[C_AIC]);
	pprintf(prn, RTFTAB "%s = %g\n", I_(bic_str), crit[C_BIC]);
	pprintf(prn, RTFTAB "%s = %g\n", I_(hqc_str), crit[C_HQC]);
    } else if (tex_format(prn)) {  
	char cval[32];

	tex_dcolumn_double(crit[C_AIC], cval);
	pprintf(prn, "%s & %s \\\\\n", I_(aic_str), cval);
	tex_dcolumn_double(crit[C_BIC], cval);
	pprintf(prn, "%s & %s \\\\\n", I_(bic_str), cval);
	tex_dcolumn_double(crit[C_HQC], cval);
	pprintf(prn, "%s & %s \\\\\n", I_(tex_hqc_str), cval);
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s (%s)\"%c%.15g\n", I_(aic_str), I_(aic_abbrev),
		prn_delim(prn), crit[C_AIC]);
	pprintf(prn, "\"%s (%s)\"%c%.15g\n", I_(bic_str), I_(bic_abbrev),
		prn_delim(prn), crit[C_BIC]);
	pprintf(prn, "\"%s (%s)\"%c%.15g\n", I_(hqc_str), I_(hqc_abbrev),
		prn_delim(prn), crit[C_HQC]);
    }	
}

static void print_liml_equation_data (const MODEL *pmod, PRN *prn)
{
    double lmin = gretl_model_get_double(pmod, "lmin");
    int idf = gretl_model_get_int(pmod, "idf");

    if (!gretl_model_get_int(pmod, "restricted")) {
	print_ll(pmod, prn);
    }

#if 0
    info_stats_lines(pmod, prn);
#endif

    if (!na(lmin)) {
	pprintf(prn, "  Smallest eigenvalue = %g\n", lmin);

	if (idf > 0) {
	    double X2 = pmod->nobs * log(lmin);

	    pprintf(prn, "  %s:\n", _("LR over-identification test"));
	    pprintf(prn, "    %s(%d) = %g %s %g\n", _("Chi-square"),
		    idf, X2, _("with p-value"), chisq_cdf_comp(X2, idf));
	} else if (idf == 0) {
	    pprintf(prn, "  %s\n", _("Equation is just identified"));
	}
    }
}

static void print_arbond_AR_test (double z, int order, PRN *prn)
{
    double pv = normal_pvalue_2(z);

    pputs(prn, "  ");

    if (plain_format(prn)) {
	pprintf(prn, _("Test for AR(%d) errors:"), order);
    } else {
	pprintf(prn, I_("Test for AR(%d) errors:"), order);
    }

    if (tex_format(prn)) {
	char numstr[32];

	tex_float_string(z, 4, numstr);
	pprintf(prn, " $z$ = %s & [%.4f]", numstr, pv);
    } else {
	pprintf(prn, " z = %g (%s %.4f)", z, 
		(plain_format(prn))? _("p-value") : I_("p-value"), pv);
    }

    gretl_prn_newline(prn);
}

enum {
    AB_SARGAN,
    AB_WALD,
    J_TEST
};

static void 
print_GMM_chi2_test (const MODEL *pmod, double x, int j, PRN *prn)
{
    const char *strs[] = {
	N_("Sargan over-identification test"),
	N_("Wald (joint) test"),
	N_("J test")
    };
    const char *texstrs[] = {
	N_("Sargan test"),
	N_("Wald (joint) test"),
	N_("J test")
    };
    double pv;
    int df;

    if (j == AB_SARGAN) {
	df = gretl_model_get_int(pmod, "sargan_df");
    } else if (j == AB_WALD) {
	df = gretl_model_get_int(pmod, "wald_df");
    } else {
	df = gretl_model_get_int(pmod, "J_df");
    }

    pv = chisq_cdf_comp(x, df);

    if (tex_format(prn)) {
	pprintf(prn, "%s: ", I_(texstrs[j]));
	pprintf(prn, "$\\chi^2(%d)$ = %g & [%.4f]", df, x, pv);
    } else if (plain_format(prn)) {
	pprintf(prn, "  %s:", _(strs[j]));
	if (j == J_TEST) {
	    pputc(prn, ' ');
	} else {
	    pputs(prn, "\n   ");
	}
	pprintf(prn, "%s(%d) = %g (%s %.4f)", _("Chi-square"),
		df, x, _("p-value"), pv);
    } else {
	pprintf(prn, "  %s:", I_(strs[j]));
	if (j == J_TEST) {
	    pputc(prn, ' ');
	} else {
	    gretl_prn_newline(prn);
	    pputs(prn, "   ");
	}
	pprintf(prn, "%s(%d) = %g (%s %.4f)", I_("Chi-square"),
		df, x, I_("p-value"), pv);
    }

    gretl_prn_newline(prn);
}

static void print_GMM_test_data (const MODEL *pmod, PRN *prn)
{
    double x;

    if (tex_format(prn)) {
	pputs(prn, "\\end{tabular}\n\n\\vspace{1ex}\n");
	pputs(prn, "\\begin{tabular}{ll}\n");
    }

    if (pmod->ci == GMM) {
	x = gretl_model_get_double(pmod, "J_test");
	if (!na(x)) {
	    print_GMM_chi2_test(pmod, x, J_TEST, prn);
	}
	return;
    }

    x = gretl_model_get_double(pmod, "AR1");
    if (!na(x)) {
	print_arbond_AR_test(x, 1, prn);
    }

    x = gretl_model_get_double(pmod, "AR2");
    if (!na(x)) {
	print_arbond_AR_test(x, 2, prn);
    }

    x = gretl_model_get_double(pmod, "sargan");
    if (!na(x)) {
	print_GMM_chi2_test(pmod, x, AB_SARGAN, prn);
    }

    x = gretl_model_get_double(pmod, "wald");
    if (!na(x)) {
	print_GMM_chi2_test(pmod, x, AB_WALD, prn);
    }
}

static void ladstats (const MODEL *pmod, PRN *prn)
{
    double ladsum = gretl_model_get_double(pmod, "ladsum");
    int utf = plain_format(prn);

    if (tex_format(prn)) {  
	char x1str[32], x2str[32];

	tex_dcolumn_double(ladsum, x1str);
	tex_dcolumn_double(pmod->ess, x2str);
	pprintf(prn, "%s & %s \\\\\n",
		I_("Sum of absolute residuals"), x1str); 
	pprintf(prn, "%s & %s \\\\\n",
		I_("Sum of squared residuals"), x2str); 
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", I_("Sum of absolute residuals"),
		prn_delim(prn), ladsum);
	pprintf(prn, "\"%s\"%c%.15g\n", I_("Sum of squared residuals"),
		prn_delim(prn), pmod->ess);
    } else {
	pprintf(prn, "  %s = %.*g\n", 
		(utf)? _("Sum of absolute residuals") :
		I_("Sum of absolute residuals"),
		GRETL_DIGITS, ladsum);
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

static void print_f_pval_str (double pval, PRN *prn)
{
    int utf = plain_format(prn);

    if (utf || rtf_format(prn)) {
	if (pval < .00001) {
	    pprintf(prn, " (%s < %.5f)\n", 
		    (utf)? _("p-value") : I_("p-value"), 0.00001);
	} else {
	    pprintf(prn, " (%s = %.3g)\n", 
		    (utf)? _("p-value") : I_("p-value"), pval);
	}
    } else if (tex_format(prn)) {
	if (pval < .00001) {
	    return;
	} else {
	    char pstr[32];

	    tex_dcolumn_double(pval, pstr);
	    pprintf(prn, "%s $F()$ & %s \\\\\n", I_("p-value for"), pstr);
	}
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%g\n", I_("p-value"), prn_delim(prn), pval);
    }
}

static void panel_variance_lines (const MODEL *pmod, PRN *prn)
{
    double ws2 = gretl_model_get_double(pmod, "within-variance");
    double bs2 = gretl_model_get_double(pmod, "between-variance");
    double theta = gretl_model_get_double(pmod, "gls-theta");

    if (na(ws2) || na(bs2)) {
	return;
    }

    if (plain_format(prn)) {
	pprintf(prn, "  %s = %g\n", _("'Within' variance"), ws2);
	pprintf(prn, "  %s = %g\n", _("'Between' variance"), bs2);
	if (!na(theta)) {
	    pprintf(prn, "  %s = %g\n", _("theta used for quasi-demeaning"), theta);
	}
    } else if (tex_format(prn)) {
	char xstr[32];

	tex_dcolumn_double(ws2, xstr);
	pprintf(prn, "$\\hat{\\sigma}^2_{\\varepsilon}$ & %s \\\\\n", xstr);
	tex_dcolumn_double(bs2, xstr);
	pprintf(prn, "$\\hat{\\sigma}^2_u$ & %s \\\\\n", xstr);
	if (!na(theta)) {
	    tex_dcolumn_double(theta, xstr);
	    pprintf(prn, "$\\theta$ & %s \\\\\n", xstr);
	}
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g", I_("'Within' variance"), ws2);
	pprintf(prn, RTFTAB "%s = %g", I_("'Between' variance"), bs2);
	if (!na(theta)) {
	    pprintf(prn, RTFTAB "%s = %g", I_("theta used for quasi-demeaning"), theta);
	}
    } else if (csv_format(prn)) {
	char d = prn_delim(prn);

	pprintf(prn, "\"%s\"%c%.15g\n", I_("'Within' variance"), d, ws2);
	pprintf(prn, "\"%s\"%c%.15g\n", I_("'Between' variance"), d, bs2);
	if (!na(theta)) {
	    pprintf(prn, "\"%s\"%c%.15g\n", I_("theta used for quasi-demeaning"), 
		    d, theta);
	}
    }	
}

static void Fline (const MODEL *pmod, PRN *prn)
{
    if (pmod->ifc && pmod->ncoeff <= 2) {
	if (plain_format(prn)) {
	    pprintf(prn, "  %s = %d\n", _("Degrees of freedom"), pmod->dfd);
	} else if (tex_format(prn)) {
	    pprintf(prn, "%s & %d \\\\\n", I_("Degrees of freedom"), pmod->dfd);
	} else if (rtf_format(prn)) {
	    pprintf(prn, RTFTAB "%s = %d\n", I_("Degrees of freedom"), pmod->dfd);
	} else if (csv_format(prn)) {
	    pprintf(prn, "\"%s\"%c%d\n", I_("Degrees of freedom"), 
		    prn_delim(prn), pmod->dfd);
	}
	return;
    }

    if (plain_format(prn)) {
	char tmp[32];

	sprintf(tmp, "%s (%d, %d)", _("F-statistic"), pmod->dfn, pmod->dfd);
	if (na(pmod->fstt)) {
	    pprintf(prn, "  %s %s\n", tmp, _("undefined"));
	} else {
	    pprintf(prn, "  %s = %.*g", tmp, XDIGITS(pmod), pmod->fstt);
	    print_f_pval_str(f_cdf_comp(pmod->fstt, pmod->dfn, pmod->dfd), prn);
	}
    } else if (tex_format(prn)) {
	if (na(pmod->fstt)) {
	    pprintf(prn, "$F(%d, %d)$ & \\multicolumn{1}{c}{\\rm %s} \\\\\n", 
		    pmod->dfn, pmod->dfd, I_("undefined"));
	} else {
	    char x1str[32];

	    tex_dcolumn_double(pmod->fstt, x1str);
	    pprintf(prn, "$F(%d, %d)$ & %s \\\\\n", pmod->dfn, pmod->dfd, x1str);
	    print_f_pval_str(f_cdf_comp(pmod->fstt, pmod->dfn, pmod->dfd), prn);
	}
    } else if (rtf_format(prn)) {
	char tmp[32];

	sprintf(tmp, "%s (%d, %d)", I_("F-statistic"), pmod->dfn, pmod->dfd);
	if (na(pmod->fstt)) {
	    pprintf(prn, RTFTAB "%s %s\n", tmp, I_("undefined"));
	} else {
	    pprintf(prn, RTFTAB "%s = %g", tmp, pmod->fstt);
	    print_f_pval_str(f_cdf_comp(pmod->fstt, pmod->dfn, pmod->dfd), prn);
	}
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s (%d, %d)\"%c", I_("F-statistic"), pmod->dfn, pmod->dfd,
		prn_delim(prn));
	if (na(pmod->fstt)) {
	    pprintf(prn, "\"%s\"\n", I_("undefined"));
	} else {
	    pprintf(prn, "%.15g%c", pmod->fstt, prn_delim(prn));
	    print_f_pval_str(f_cdf_comp(pmod->fstt, pmod->dfn, pmod->dfd), prn);
	}
    }	
}

static void dwline (const MODEL *pmod, PRN *prn)
{
    if (na(pmod->dw)) {
	return;
    }

    if (plain_format(prn)) {
	pprintf(prn, "  %s = %.*g\n", _("Durbin-Watson statistic"), 
		XDIGITS(pmod), pmod->dw);
	if (!na(pmod->rho)) {
	    pprintf(prn, "  %s = %.*g\n", _("First-order autocorrelation coeff."), 
		    XDIGITS(pmod), pmod->rho);
	}
    } else if (tex_format(prn)) {
	char xstr[32];

	tex_dcolumn_double(pmod->dw, xstr);
	pprintf(prn, "%s & %s \\\\\n",
		I_("Durbin--Watson statistic"), xstr); 
	if (!na(pmod->rho)) {
	    tex_dcolumn_double(pmod->rho, xstr);
	    pprintf(prn, "%s & %s \\\\\n",
		    I_("First-order autocorrelation coeff."), xstr);
	}
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_("Durbin-Watson statistic"), 
		pmod->dw);
	if (!na(pmod->rho)) {
	    pprintf(prn, RTFTAB "%s = %g\n", I_("First-order autocorrelation coeff."), 
		    pmod->rho);
	} 
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", I_("Durbin-Watson statistic"), 
		prn_delim(prn), pmod->dw);
	if (!na(pmod->rho)) {
	    pprintf(prn, "\"%s\"%c%.15g\n", I_("First-order autocorrelation coeff."), 
		    prn_delim(prn), pmod->rho);
	}
    }
}

static void dhline (const MODEL *pmod, PRN *prn)
{
    double sderr, h = 0.0;
    int ldepvar = gretl_model_get_int(pmod, "ldepvar");
    int T = pmod->nobs - 1;

    sderr = pmod->sderr[ldepvar - 2];

    if (pmod->ess <= 0.0 || (T * sderr * sderr) >= 1.0) return;

    h = pmod->rho * sqrt(T/(1 - T * sderr * sderr));

    if (plain_format(prn)) {
        char tmp[128];

        sprintf(tmp, _("Durbin's h stat. %g"), h);
        pprintf(prn, "  %s\n", tmp);

        sprintf(tmp, _("(Using variable %d for h stat, with T' = %d)"), 
                pmod->list[ldepvar], T);
        pprintf(prn, "  %s\n", tmp);
    } else if (rtf_format(prn)) {
        char tmp[128];

        sprintf(tmp, I_("Durbin's h stat. %g"), h);
        pprintf(prn, RTFTAB "%s\n", tmp);

        sprintf(tmp, I_("(Using variable %d for h stat, with T' = %d)"), 
                pmod->list[ldepvar], T);
        pprintf(prn, RTFTAB "%s\n", tmp);
    } else if (tex_format(prn)) {
	char x1str[32];

	tex_dcolumn_double(h, x1str);
	pprintf(prn, "%s & %s \\\\\n",
		I_("Durbin's $h$ statistic"), x1str);
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", I_("Durbin's h"), prn_delim(prn), h);
    }
}

static int least_signif_coeff (const MODEL *pmod)
{
    double tstat, tmin = 4.0;
    int i, k = 0;
    
    for (i=pmod->ifc; i<pmod->ncoeff; i++) {
	tstat = fabs(pmod->coeff[i] / pmod->sderr[i]);
	if (tstat < tmin) {
	    tmin = tstat;
	    k = i;
	}
    }

    if (coeff_pval(pmod->ci, tmin, pmod->dfd) > .10) {
	return pmod->list[k+2];
    }

    return 0;
}

static void pval_max_line (const MODEL *pmod, const DATAINFO *pdinfo, 
			   PRN *prn)
{
    int k = pmod->ncoeff - pmod->ifc;

    if (k < 3) return;

    if ((k = least_signif_coeff(pmod))) {
	char tmp[128];

	if (pmod->ifc) {
	    sprintf(tmp, _("Excluding the constant, p-value was highest "
			   "for variable %d (%s)"), k, pdinfo->varname[k]);
	} else {
	    sprintf(tmp, _("P-value was highest for variable %d (%s)"), 
		    k, pdinfo->varname[k]);
	}	    
	pprintf(prn, "%s\n\n", tmp);
    }
}

static const char *aux_string (int aux, PRN *prn)
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
	if (tex_format(prn)) return N_("Cointegrating regression -- ");
	else return N_("Cointegrating regression - ");
    } else if (aux == AUX_ADF) {
	if (tex_format(prn)) return N_("Augmented Dickey--Fuller regression");
	else return N_("Augmented Dickey-Fuller regression");
    } else if (aux == AUX_DF) {
	if (tex_format(prn)) return N_("Dickey--Fuller regression");
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

static const char *simple_estimator_string (int ci, PRN *prn)
{
    if (ci == OLS || ci == VAR) return N_("OLS");
    else if (ci == WLS)  return N_("WLS"); 
    else if (ci == ARCH) return N_("WLS (ARCH)");
    else if (ci == TSLS) return N_("TSLS");
    else if (ci == HSK)  return N_("Heteroskedasticity-corrected");
    else if (ci == AR)   return N_("AR");
    else if (ci == LAD)  return N_("LAD");
    else if (ci == MPOLS) return N_("High-Precision OLS");
    else if (ci == PROBIT) return N_("Probit");
    else if (ci == LOGIT)  return N_("Logit");
    else if (ci == TOBIT)  return N_("Tobit");
    else if (ci == HECKIT) return N_("Heckit");
    else if (ci == POISSON) return N_("Poisson");
    else if (ci == NLS) return N_("NLS");
    else if (ci == MLE) return N_("ML");
    else if (ci == GMM) return N_("GMM");
    else if (ci == LOGISTIC) return N_("Logistic");
    else if (ci == GARCH) return N_("GARCH");
    else if (ci == CORC) {
	if (tex_format(prn)) return N_("Cochrane--Orcutt");
	else return N_("Cochrane-Orcutt");
    } else if (ci == HILU) {
	if (tex_format(prn)) return N_("Hildreth--Lu");
	else return N_("Hildreth-Lu");
    } else if (ci == PWE) {
	if (tex_format(prn)) return N_("Prais--Winsten");
	else return N_("Prais-Winsten");
    } else if (ci == ARBOND) {
	if (tex_format(prn)) return N_("Arellano--Bond");
	else return N_("Arellano-Bond");
    } else {
	return "";
    }
}

const char *estimator_string (const MODEL *pmod, PRN *prn)
{
    if (pmod->ci == ARMA) {
	if (gretl_model_get_int(pmod, "armax")) {
	    return N_("ARMAX");
	} else if (gretl_model_get_int(pmod, "arima_d") ||
		   gretl_model_get_int(pmod, "arima_D")) {
	    return N_("ARIMA");
	} else {
	    return N_("ARMA");
	}
    } else if (pmod->ci == WLS) {
	if (gretl_model_get_int(pmod, "iters")) {
	    return N_("Maximum Likelihood");
	} else {
	    return N_("WLS");
	}
    } else if (pmod->ci == PANEL) {
	if (gretl_model_get_int(pmod, "fixed-effects")) {
	    return N_("Fixed-effects");
	} else if (gretl_model_get_int(pmod, "random-effects")) {
	    return N_("Random-effects (GLS)");
	} else {
	    return N_("Between-groups");
	}
    } else if (pooled_model(pmod)) {
	return N_("Pooled OLS");
    } else if (pmod->ci == ARBOND) {
	if (gretl_model_get_int(pmod, "step") == 2) {
	    return N_("2-step Arellano-Bond");
	} else {
	    return N_("1-step Arellano-Bond");
	}
    } else if (pmod->ci == GMM) {
	int s = gretl_model_get_int(pmod, "step");

	if (s > 2) {
	    return N_("Iterated GMM");
	} else if (s == 2) {
	    return N_("2-step GMM");
	} else {
	    return N_("1-step GMM");
	}	
    } else if (pmod->ci == LOGIT) {
	if (gretl_model_get_int(pmod, "ordered")) {
	    return N_("Ordered Logit");
	} else {
	    return N_("Logit");
	}
    } else if (pmod->ci == PROBIT) {
	if (gretl_model_get_int(pmod, "ordered")) {
	    return N_("Ordered Probit");
	} else {
	    return N_("Probit");
	}
    } else {
	return simple_estimator_string(pmod->ci, prn);
    }
}

static int any_tests (const MODEL *pmod)
{
    if (pmod->ntests > 0) {
	return 1;
    }

    if (pmod->ci == TSLS && gretl_model_get_int(pmod, "stage1-dfn")) {
	return 1;
    }

    return 0;
}

static void maybe_print_first_stage_F (const MODEL *pmod, PRN *prn)
{
    double F = gretl_model_get_double(pmod, "stage1-F");
    int dfn, dfd;

    if (na(F)) {
	return;
    }

    dfn = gretl_model_get_int(pmod, "stage1-dfn");
    dfd = gretl_model_get_int(pmod, "stage1-dfd");

    if (dfn <= 0 || dfd <= 0) {
	return;
    }

    if (plain_format(prn)) {
	pprintf(prn, "%s (%d, %d) = %.*g\n", _("First-stage F-statistic"),
		dfn, dfd, GRETL_DIGITS, F);
	pprintf(prn, "  %s\n\n", _("A value < 10 may indicate weak instruments"));
    } else if (tex_format(prn)) {
	char x1str[32];

	tex_dcolumn_double(F, x1str);
	pprintf(prn, "First-stage $F(%d, %d)$ = %s \\\\\n", dfn, dfd, x1str);
    } else if (rtf_format(prn)) {
	pprintf(prn, "%s (%d, %d) = %g\n", I_("First-stage F-statistic"), 
		dfn, dfd, F);
    }
}

static void print_model_tests (const MODEL *pmod, PRN *prn)
{
    int i;

    if (tex_format(prn)) {
	pputs(prn, "\\vspace{1em}\n\\begin{raggedright}\n");
	for (i=0; i<pmod->ntests; i++) {
	    if (i > 0) {
		pputs(prn, "\\vspace{1ex}\n");
	    }
	    gretl_model_test_print(pmod, i, prn);
	}
	if (pmod->ntests > 0) {
	    pputs(prn, "\\vspace{1ex}\n");
	}
	maybe_print_first_stage_F(pmod, prn);
	pputs(prn, "\\end{raggedright}\n");
    } else {
	for (i=0; i<pmod->ntests; i++) {
	    gretl_model_test_print(pmod, i, prn);
	}
	maybe_print_first_stage_F(pmod, prn);
    }
}

static int 
print_tsls_instruments (const int *list, const DATAINFO *pdinfo, PRN *prn)
{
    int i, j, pos = 0;
    int ccount = 0;
    int ninst = 0;
    char vname[16];
    int tex = tex_format(prn);
    int utf = plain_format(prn);

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
		if (ccount >= 64) {
		    if (tex) {
			pputs(prn, "\\\\\n");
		    } else if (rtf_format(prn)) {
			pputs(prn, "\\par\n");
		    } else {
			pputs(prn, "\n  "); 
		    }
		    ccount = 0;
		}
		ninst++;
	    }
	}
    }

    if (ninst == 0) {
	pputs(prn, _("none"));
    }

    if (ccount > 0) {
	if (tex) {
	    pputs(prn, "\\\\\n");
	} else if (rtf_format(prn)) {
	    pputs(prn, "\\par\n");
	} else {
	    pputs(prn, "\n");
	}
    }

    return 0;
}

static void arbond_asy_vcv_line (const MODEL *pmod, PRN *prn)
{
    if (csv_format(prn)) {
	pprintf(prn, "\"%s\"", I_("Asymptotic standard errors (unreliable)"));
    } else if (plain_format(prn)) {
	pputs(prn, _("Asymptotic standard errors (unreliable)"));
    } else {
	pputs(prn, I_("Asymptotic standard errors (unreliable)"));
    } 

    pputc(prn, '\n');
}


static void panel_robust_vcv_line (PRN *prn)
{
    if (csv_format(prn)) {
	pprintf(prn, "\"%s\"", I_("Robust (HAC) standard errors"));
    } else if (plain_format(prn)) {
	pputs(prn, _("Robust (HAC) standard errors"));
    } else {
	pputs(prn, I_("Robust (HAC) standard errors"));
    } 

    pputc(prn, '\n');
}

static void beck_katz_vcv_line (PRN *prn)
{
    if (csv_format(prn)) {
	pprintf(prn, "\"%s\"", I_("Beck-Katz standard errors"));
    } else if (plain_format(prn)) {
	pputs(prn, _("Beck-Katz standard errors"));
    } else if (tex_format(prn)) {
	pputs(prn, I_("Beck--Katz standard errors"));
    } else {
	pputs(prn, I_("Beck-Katz standard errors"));
    } 

    pputc(prn, '\n');
}

static void beck_katz_failed_line (PRN *prn)
{
    if (plain_format(prn)) {
	pputs(prn, _("Could not compute Beck-Katz standard errors"));
	pputc(prn, '\n');
    }
}

static void hac_vcv_line (const MODEL *pmod, PRN *prn)
{
    const char *kstrs[] = {
	N_("Bartlett kernel"),
	N_("Parzen kernel"),
	N_("QS kernel")
    };
    int kern = gretl_model_get_int(pmod, "hac_kernel");
    int h = gretl_model_get_int(pmod, "hac_lag");
    int white = gretl_model_get_int(pmod, "hac_prewhiten");
    double bt;

    if (kern == KERNEL_QS) {
	bt = gretl_model_get_double(pmod, "qs_bandwidth");
	if (plain_format(prn)) {
	    pprintf(prn, _("HAC standard errors, "
			   "bandwidth %.2f"), bt);
	} else {
	    pprintf(prn, I_("HAC standard errors, "
			    "bandwidth %.2f"), bt);
	}
    } else {
	if (plain_format(prn)) {
	    pprintf(prn, _("HAC standard errors, "
			   "bandwidth %d"), h);
	} else {
	    pprintf(prn, I_("HAC standard errors, "
			    "bandwidth %d"), h);
	}
    }

    pputc(prn, ' ');
    if (plain_format(prn)) {
	pprintf(prn, "(%s", _(kstrs[kern - 1]));
    } else {
	pprintf(prn, "(%s", I_(kstrs[kern - 1]));
		    
    }
    if (white) {
	pputs(prn, ", ");
	pputs(prn, (plain_format(prn))? _("prewhitened") : 
	      I_("prewhitened"));
    }

    pputs(prn, ")\n");
}

static void hc_vcv_line (const MODEL *pmod, PRN *prn)
{
    int hcv = gretl_model_get_int(pmod, "hc_version");
    int jack = 0;

    if (hcv == 4) {
	jack = 1;
	hcv--;
    }

    if (plain_format(prn)) {
	pprintf(prn, "%s, %s%sHC%d%s", 
		_("Heteroskedasticity-robust standard errors"),
		(jack)? "" : _("variant"),
		(jack)? "" : " ",
		hcv, (jack)? " (jackknife)" : "");
    } else {
	pprintf(prn, "%s, %s%sHC%d%s", 
		I_("Heteroskedasticity-robust standard errors"),
		(jack)? "" : I_("variant"),
		(jack)? "" : " ",
		hcv, (jack)? " (jackknife)" : "");
    }

    if (rtf_format(prn)) {
	pputs(prn, "\\par\n");
    } else {
	pputc(prn, '\n');
    }
}

static void ml_vcv_line (const MODEL *pmod, PRN *prn)
{
    int v = gretl_model_get_int(pmod, "ml_vcv");
    int tex = tex_format(prn);
    int utf = plain_format(prn);
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
	if (csv_format(prn)) {
	    pprintf(prn, "\"%s\"\n", I_(vcvstr));
	} else {
	    pprintf(prn, "%s\n", (utf)? _(vcvstr) : I_(vcvstr));
	}
    }
}

static void tex_vecm_depvar_name (char *s, const char *vname)
{
    char tmp[14];
    int gotit = 0;

    if (sscanf(vname, "d_%13s", tmp)) {
	char myvar[24];

	tex_escape(myvar, tmp);
	sprintf(s, "$\\Delta$%s", myvar);
	gotit = 1;
    }

    if (!gotit) {
	tex_escape(s, vname); 
    }    
}

static void tex_arbond_depvar_name (char *s, const char *vname)
{
    char vnesc[32];

    tex_escape(vnesc, vname);
    sprintf(s, "$\\Delta$%s", vnesc);
}

void print_model_vcv_info (const MODEL *pmod, PRN *prn)
{
    if (gretl_model_get_int(pmod, "hac_kernel")) {
	hac_vcv_line(pmod, prn);
    } else if (gretl_model_get_int(pmod, "hc")) {
	hc_vcv_line(pmod, prn);
    } else if (gretl_model_get_int(pmod, "ml_vcv")) {
	ml_vcv_line(pmod, prn);
    } else if (gretl_model_get_int(pmod, "panel_hac")) {
	panel_robust_vcv_line(prn);
    } else if (gretl_model_get_int(pmod, "panel_bk")) {
	beck_katz_vcv_line(prn);
    } else if (gretl_model_get_int(pmod, "panel_bk_failed")) {
	beck_katz_failed_line(prn);
    } else if (pmod->ci == ARBOND && gretl_model_get_int(pmod, "asy")) {
	arbond_asy_vcv_line(pmod, prn);
    }
}

static void print_extra_list (const char *tag, const int *list, 
			      const DATAINFO *pdinfo, PRN *prn)
{
    int i, v, len;

    len = pputs(prn, _(tag));

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v < pdinfo->v) {
	    len += pprintf(prn, " %s", pdinfo->varname[v]);
	} else {
	    len += pprintf(prn, " %d", v);
	}
	if (len > 68 && i < list[0]) {
	    pputc(prn, '\n');
	    len = 0;
	}
    }

    pputc(prn, '\n');
}

static void print_model_zerolist (const MODEL *pmod, 
				  const DATAINFO *pdinfo,
				  PRN *prn)
{
    const int *zlist = gretl_model_get_data(pmod, "zerolist");
    const char *tag = N_("Omitted because all values were zero:");

    if (pmod->ci == PANEL && gretl_model_get_int(pmod, "between")) {
	return;
    }

    print_extra_list(tag, zlist, pdinfo, prn);
}

static void print_model_droplist (const MODEL *pmod, 
				  const DATAINFO *pdinfo,
				  PRN *prn)
{
    const int *dlist = gretl_model_get_data(pmod, "droplist");
    const char *tag = N_("Omitted due to exact collinearity:");

    if (pmod->ci == PANEL && gretl_model_get_int(pmod, "between")) {
	return;
    }

    print_extra_list(tag, dlist, pdinfo, prn);
}

static void print_tsls_droplist (const MODEL *pmod, 
				 const DATAINFO *pdinfo,
				 PRN *prn)
{
    const int *dlist = gretl_model_get_data(pmod, "inst_droplist");
    int i, v;

    pputs(prn, _("Redundant instruments:"));
    for (i=1; i<=dlist[0]; i++) {
	v = dlist[i];
	if (v < pdinfo->v) {
	    pprintf(prn, " %s", pdinfo->varname[v]);
	} else {
	    pprintf(prn, " %d", v);
	}
    }
    pputc(prn, '\n');
}

static void print_arma_depvar (const MODEL *pmod,
			       const DATAINFO *pdinfo,
			       PRN *prn)
{
    int tex = tex_format(prn);
    int utf = plain_format(prn);
    int yno = gretl_model_get_depvar(pmod);
    int d = gretl_model_get_int(pmod, "arima_d");
    int D = gretl_model_get_int(pmod, "arima_D");
    char vname[64];

    *vname = 0;

    if (tex) {
	char tmp[32];

	if (d > 0 || D > 0) {
	    strcat(vname, "$");
	}
	if (d == 1) {
	    strcat(vname, "(1-L)");
	} else if (d == 2) {
	    strcat(vname, "(1-L)^2");
	}
	if (D == 1) {
	    strcat(vname, "(1-L^s)");
	} else if (D == 2) {
	    strcat(vname, "(1-L^s)^2");
	}
	if (d > 0 || D > 0) {
	    strcat(vname, "$");
	}
	tex_escape(tmp, pdinfo->varname[yno]);
	strcat(vname, tmp);
    } else {
	if (d == 1) {
	    strcat(vname, "(1-L)");
	} else if (d == 2) {
	    strcat(vname, "(1-L)^2");
	}
	if (D == 1) {
	    strcat(vname, "(1-Ls)");
	} else if (D == 2) {
	    strcat(vname, "(1-Ls)^2");
	}
	if (d > 0 || D > 0) {
	    strcat(vname, " ");
	}
	strcat(vname, pdinfo->varname[yno]);
    }

    pprintf(prn, "%s: %s", 
	    (utf)? _("Dependent variable") : I_("Dependent variable"),
	    vname);
}

static void arma_extra_info (const MODEL *pmod, PRN *prn)
{
    int acode = gretl_model_get_int(pmod, "arma_flags");

    if (acode & ARMA_X12A) {
	pputs(prn, _("Estimated using X-12-ARIMA"));
	pputs(prn, " (");
	pputs(prn, (acode & ARMA_EXACT)? _("exact ML") : _("conditional ML"));
	pputs(prn, ")\n");
    } else if (acode & ARMA_EXACT) {
	pputs(prn, _("Estimated using Kalman filter"));
	pputs(prn, " (");
	pputs(prn, _("exact ML"));
	pputs(prn, ")\n");
    } else {
	pputs(prn, _("Estimated using BHHH method"));
	pputs(prn, " (");
	pputs(prn, _("conditional ML"));
	pputs(prn, ")\n");
    }	
}

static void print_model_heading (const MODEL *pmod, 
				 const DATAINFO *pdinfo, 
				 gretlopt opt, 
				 PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN], vname[32];
    int t1 = pmod->t1, t2 = pmod->t2;
    int tex = tex_format(prn);
    int utf = plain_format(prn);
    int csv = csv_format(prn);
    int dvnl = 1;
    int order = 0;

    if (pmod->aux != AUX_VAR && pmod->aux != AUX_VECM) {
	ntodate(startdate, t1, pdinfo);
	ntodate(enddate, t2, pdinfo);
    }

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
	    pprintf(prn, "\n%s\n", _(aux_string(pmod->aux, prn)));
	} else if (tex) {
	    pprintf(prn, "\n%s\n", I_(aux_string(pmod->aux, prn)));
	} else if (csv) {
	    pprintf(prn, "\"%s\"\n", I_(aux_string(pmod->aux, prn)));
	} else { /* RTF */
	    pprintf(prn, "%s\\par\n", I_(aux_string(pmod->aux, prn)));
	}
	break;
    case AUX_AR:
	order = gretl_model_get_int(pmod, "BG_order");
	if (utf) { 	
	    pprintf(prn, "\n%s ", _("Breusch-Godfrey test for"));
	} else {
	    pprintf(prn, "\n%s ", I_("Breusch-Godfrey test for"));
	} 
	if (order > 1) {
	    pprintf(prn, "%s %d\n", (utf)? _("autocorrelation up to order") :
		    I_("autocorrelation up to order"), 
		    order);
	} else {
	    pprintf(prn, "%s\n", (utf)? _("first-order autocorrelation") :
		    I_("first-order autocorrelation"));
	}
	break;	
    case AUX_ARCH:
	order = gretl_model_get_int(pmod, "arch_order");
	pprintf(prn, "\n%s %d\n", 
		(utf)? _("Test for ARCH of order") : 
		I_("Test for ARCH of order"), 
		order);
	break;	
    case AUX_SYS:
	pprintf(prn, "%s %d: ", 
		(utf)? _("Equation") : I_("Equation"), pmod->ID + 1);
	break;	
    case AUX_VAR:
	pprintf(prn, "\n%s %d: ", 
		(utf)? _("Equation") : I_("Equation"), pmod->ID);
	break;
    case AUX_VECM:
	pprintf(prn, "%s %d: ", 
		(utf)? _("Equation") : I_("Equation"), pmod->ID);
	break;
    case AUX_ADD:
    default:
	if (pmod->ID < 0 || (opt & OPT_S)) {
	    if (!csv) {
		pputc(prn, '\n');
	    }
	} else if (pmod->name) {
	    if (csv) {
		pprintf(prn, "\"%s:\"\n", pmod->name);
	    } else {
		pprintf(prn, "\n%s:\n", pmod->name);
	    }
	} else {
	    if (csv) {
		pprintf(prn, "\"%s %d: ", I_("Model"), pmod->ID);
	    } else {
		pprintf(prn, "\n%s %d: ", (utf)? _("Model") : I_("Model"), pmod->ID);
	    }
	}
	break;
    }

    if (pmod->aux == AUX_VAR || pmod->aux == AUX_VECM) {
	;
    } else if (pmod->aux == AUX_SYS) {
	pprintf(prn, (utf)?
		_("%s estimates using the %d observations %s%s%s") :
		I_("%s estimates using the %d observations %s%s%s"),
		_(gretl_system_short_string(pmod)),
		pmod->nobs, startdate, (tex)? "--" : "-", enddate);
    } else if (!dataset_is_panel(pdinfo)) {
	if (pmod->missmask != NULL) {
	    int mc = model_missval_count(pmod);

	    if (pmod->ci == HECKIT) {
		int Tmax = pmod->t2 - pmod->t1 + 1;
		
		mc = Tmax - gretl_model_get_int(pmod, "totobs");
	    }

	    pprintf(prn, (utf)?
		    _("%s estimates using %d observations from %s%s%s") :
		    I_("%s estimates using %d observations from %s%s%s"),
		    _(estimator_string(pmod, prn)), 
		    pmod->nobs, startdate, (tex)? "--" : "-", enddate);
	    if (mc > 0) {
		gretl_prn_newline(prn);
		pprintf(prn, "%s: %d",
			(utf)? _("Missing or incomplete observations dropped") :
			I_("Missing or incomplete observations dropped"), mc);
	    }
	} else {
	    pprintf(prn, (utf)?
		    _("%s estimates using the %d observations %s%s%s") :
		    I_("%s estimates using the %d observations %s%s%s"),
		    _(estimator_string(pmod, prn)), 
		    pmod->nobs, startdate, (tex)? "--" : "-", enddate);
	}
    } else {
	int effn = gretl_model_get_int(pmod, "n_included_units");
	int Tmin = gretl_model_get_int(pmod, "Tmin");
	int Tmax = gretl_model_get_int(pmod, "Tmax");

	pprintf(prn, (utf)?
		_("%s estimates using %d observations") :
		I_("%s estimates using %d observations"),
		_(estimator_string(pmod, prn)), 
		pmod->nobs);
	if (effn > 0) {
	    gretl_prn_newline(prn);
	    pprintf(prn, (utf)? _("Included %d cross-sectional units") :
		    I_("Included %d cross-sectional units"), effn);
	}
	if (Tmin > 0 && Tmax > 0) {
	    gretl_prn_newline(prn);
	    if (Tmin == Tmax) {
		pprintf(prn, (utf)? _("Time-series length = %d") :
			I_("Time-series length = %d"), Tmin);
	    } else {
		pprintf(prn, (utf)? _("Time-series length: minimum %d, maximum %d") :
			I_("Time-series length: minimum %d, maximum %d"), 
			Tmin, Tmax);
	    }
	}
    }

    if (csv) pputc(prn, '"');

    if (pmod->aux != AUX_VAR && pmod->aux != AUX_VECM) {
	gretl_prn_newline(prn);
    }

    if (pmod->ci == ARMA && plain_format(prn)) {
	arma_extra_info(pmod, prn);
    }

    if (csv) pputc(prn, '"');

    /* special formulations for dependent variable in various cases */
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG) {
	pprintf(prn, "%s: %s", 
		(utf)? _("Dependent variable") : I_("Dependent variable"),
		(tex)? "$\\hat{u}$" : "uhat");
    } else if (pmod->aux == AUX_WHITE) {
	pprintf(prn, "%s: %s", 
		(utf)? _("Dependent variable") : I_("Dependent variable"),
		(tex)? "$\\hat{u}^2$" : "uhat^2");
    } else if (pmod->aux == AUX_ARCH) {
	pprintf(prn, "%s: %s", 
		(utf)? _("Dependent variable") : I_("Dependent variable"),
		(tex)? "$u_t^2$" : "ut^2");
    } else if (pmod->ci == NLS) {
	if (tex) tex_escape(vname, pmod->depvar);
	pprintf(prn, "%s: %s", 
		(utf)? _("Dependent variable") : I_("Dependent variable"),
		(tex)? vname : pmod->depvar);
    } else if (pmod->ci == MLE || pmod->ci == GMM) {
	if (pmod->depvar != NULL) {
	    if (tex) {
		pprintf(prn, "\\verb!%s!", pmod->depvar);
	    } else {
		pputs(prn, pmod->depvar);
	    }
	} else {
	    dvnl = 0;
	}
    } else if (pmod->ci == ARMA) {
	print_arma_depvar(pmod, pdinfo, prn);
    } else { 
	int v = gretl_model_get_depvar(pmod);

	if (tex) {
	    if (pmod->aux == AUX_VECM) {
		tex_vecm_depvar_name(vname, pdinfo->varname[v]);
	    } else if (pmod->ci == ARBOND) {
		tex_arbond_depvar_name(vname, pdinfo->varname[v]);
	    } else {
		tex_escape(vname, pdinfo->varname[v]);
	    }
	}

	if (pmod->aux == AUX_VAR || pmod->aux == AUX_VECM) {
	    pputs(prn, (tex)? vname : pdinfo->varname[v]);
	} else {
	    pprintf(prn, "%s: %s", 
		    (utf)? _("Dependent variable") : I_("Dependent variable"),
		    (tex)? vname : pdinfo->varname[v]);
	}
    }

    if (csv) pputc(prn, '"');

    if (dvnl) {
	gretl_prn_newline(prn);
    }

    /* supplementary strings below the estimator and sample info */

    /* list of instruments for TSLS */
    if (pmod->ci == TSLS) {
	int method = gretl_model_get_int(pmod, "method");

	if (method != SYS_FIML && method != SYS_LIML) {
	    print_tsls_instruments(pmod->list, pdinfo, prn);
	}
    }

    /* VCV variants */
    print_model_vcv_info(pmod, prn);

    /* WLS on panel data */
    if (gretl_model_get_int(pmod, "unit_weights") && !pmod->aux) {
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

    /* weight variable for WLS */
    else if ((pmod->ci == WLS && !pmod->aux)) {
	if (tex) {
	    tex_escape(vname, pdinfo->varname[pmod->nwt]);
	}
	if (csv) pputc(prn, '"');
	pprintf(prn, "%s: %s", 
		(utf)? _("Variable used as weight") : I_("Variable used as weight"), 
		(tex)? vname : pdinfo->varname[pmod->nwt]);
	if (csv) pputc(prn, '"');
	pputc(prn, '\n');
    }

    /* weight variable for ARCH */
    else if (pmod->ci == ARCH) {
	if (csv) pputc(prn, '"');
	pprintf(prn, "%s: %s", 
		(utf)? _("Variable used as weight") : I_("Variable used as weight"), 
		(tex)? "$1/\\hat{\\sigma}_t$" : "1/sigma");
	if (csv) pputc(prn, '"');
	pputc(prn, '\n');
    }    

    /* rhohat for CORC and HILU (TeX) */
    else if (pmod->ci == CORC || pmod->ci == HILU || pmod->ci == PWE) {
	if (tex) {
	    pprintf(prn, "$\\hat{\\rho}$ = %g\n", 
		    gretl_model_get_double(pmod, "rho_in"));
	}
    } 

    /* y-hat formula for logistic regression */
    else if (pmod->ci == LOGISTIC) {
	if (tex) {
	    pprintf(prn, "$\\hat{y} = %g / (1 + e^{-X\\hat{\\beta}})$\n", 
		    gretl_model_get_double(pmod, "lmax"));  
	} else {
	    pprintf(prn, "yhat = %g / (1 + exp(-X*b))\n",  
		    gretl_model_get_double(pmod, "lmax"));
	}
    }

    /* TSLS: message about redundant instruments */
    if (plain_format(prn) && pmod->ci == TSLS &&
	gretl_model_get_data(pmod, "inst_droplist") != NULL) {
	print_tsls_droplist(pmod, pdinfo, prn);
    }  

    /* messages about collinear and/or zero regressors */
    if (plain_format(prn)) {
	if (gretl_model_get_data(pmod, "zerolist") != NULL) {
	    print_model_zerolist(pmod, pdinfo, prn);
	}
	if (gretl_model_get_data(pmod, "droplist") != NULL) {
	    print_model_droplist(pmod, pdinfo, prn);
	}
    }    

    if (pmod->missmask == NULL && gretl_model_get_int(pmod, "wt_dummy")) { 
	/* FIXME alt formats */
	pprintf(prn, "%s %d\n", 
		(utf)? _("Weight var is a dummy variable, effective obs =") :
		I_("Weight var is a dummy variable, effective obs ="),
		pmod->nobs);
    } 

    if (rtf_format(prn)) {
	pputs(prn, "\\par\n");
    } else {
	pputc(prn, '\n');
    }
}

static void model_format_start (PRN *prn)
{
    if (tex_format(prn)) {
	if (tex_doc_format(prn)) {
	    gretl_tex_preamble(prn, 0);
	} else {
	    pputs(prn, "%% You'll need to \\usepackage{dcolumn}\n\n");
	}
	pputs(prn, "\\begin{center}\n");
    } else if (rtf_format(prn)) {
	if (rtf_doc_format(prn)) {
	    pputs(prn, "{\\rtf1\\par\n\\qc ");
	} else {
	    pputs(prn, "\\par\n\\qc ");
	}
    }
}

#define RTF_COEFF_ROW  "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                       "\\cellx1900\\cellx3300\\cellx4700\\cellx6100" \
                       "\\cellx7500\\cellx8000\n\\intbl"

#define RTF_BINARY_ROW "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                       "\\cellx1900\\cellx3300\\cellx4700\\cellx6100" \
                       "\\cellx8000\n\\intbl"

#define RTF_ROOT_ROW   "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                       "\\cellx500\\cellx1500\\cellx2900\\cellx4300" \
                       "\\cellx5700\\cellx7100\n\\intbl"

void print_coeff_heading (int mode, PRN *prn)
{
    if (mode == COEFF_HEADING_VARNAME) {
	pputs(prn, _("      VARIABLE       COEFFICIENT        STDERROR"
		     "      T STAT   P-VALUE\n\n"));
    } else {
	pputs(prn, _("      PARAMETER       ESTIMATE          STDERROR"
		     "      T STAT   P-VALUE\n\n"));
    } 
}

static void print_coeff_table_start (const MODEL *pmod, PRN *prn)
{
    int slopes = binary_model(pmod) && !gretl_model_get_int(pmod, "show-pvals");
    int use_param = pmod->ci == NLS || pmod->ci == MLE || pmod->ci == GMM;

    if (plain_format(prn)) {
	if (pmod->ci == MPOLS) {
	    pputs(prn, _("      VARIABLE            COEFFICIENT          "
		       "        STDERROR\n"));
	    pputc(prn, '\n');
	} else if (slopes) {
	    pputs(prn, _("      VARIABLE       COEFFICIENT        STDERROR"
			 "      T STAT       SLOPE\n"));
	    pprintf(prn, "                                                 "
		    "                 %s\n", _("(at mean)"));
	} else {
	    print_coeff_heading(use_param, prn);
	}
	return;
    } else if (csv_format(prn)) {
	char d = prn_delim(prn);

	if (pmod->ci == MPOLS) {
	    pprintf(prn, "\"%s\"%c\"%s\"%c\"%s\"\n",
		    I_("VARIABLE"), d, I_("COEFFICIENT"), d, I_("STDERROR"));
	} else if (slopes) {
	    pprintf(prn, "\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    I_("VARIABLE"), d, I_("COEFFICIENT"), d, I_("STDERROR"),
		    d, I_("T STAT"), d, I_("SLOPE at mean"));
	} else if (use_param) {
	    pprintf(prn, "\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    I_("PARAMETER"), d, I_("ESTIMATE"), d, I_("STDERROR"),
		    d, I_("T STAT"), d, I_("P-VALUE"));
	} else {
	    pprintf(prn, "\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    I_("VARIABLE"), d, I_("COEFFICIENT"), d, I_("STDERROR"),
		    d, I_("T STAT"), d, I_("P-VALUE"));
	}	
    } else {
	char col1[16], col2[16];

	if (use_param) {
	    strcpy(col1, N_("Parameter"));
	    strcpy(col2, N_("Estimate"));
	} else {
	    strcpy(col1, N_("Variable"));
	    strcpy(col2, N_("Coefficient"));
	}	    

	if (tex_format(prn)) {
	    tex_coeff_table_start(col1, col2, slopes, prn);
	    return;
	}   

	if (rtf_format(prn)) {
	    if (slopes) {
		pprintf(prn, "{" RTF_BINARY_ROW
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

static void print_coeff_table_end (const MODEL *pmod, PRN *prn)
{
    if (plain_format(prn) || csv_format(prn)) {
	pputc(prn, '\n');
    } else if (tex_format(prn)) {
	tex_coeff_table_end(prn);
    } else if (rtf_format(prn)) {
	pputs(prn, "}\n\n");
    }

    if (plain_format(prn) && gretl_model_get_int(pmod, "near-singular")) {
	const char *msg = N_("Warning: data matrix close to singularity!");

	pprintf(prn, "%s\n\n", _(msg));
    }
}

static void model_format_end (PRN *prn)
{
    if (tex_format(prn)) {
	pputs(prn, "\n\\end{center}\n");
	if (tex_doc_format(prn)) {
	    pputs(prn, "\n\\end{document}\n"); 
	}
    } else if (rtf_doc_format(prn)) {
	pputs(prn, "\n}\n");
    }
} 


static void print_middle_table_start (PRN *prn)
{
    if (tex_format(prn)) {
	char pt = get_local_decpoint();

	pprintf(prn, 
		"\\vspace{1em}\n\n"
		"\\begin{tabular}{lD{%c}{%c}{-1}}\n",
		pt, pt);
    }
}

static void print_middle_table_end (PRN *prn)
{
    if (tex_format(prn)) {
	pputs(prn, "\\end{tabular}\n\n");
    } else if (rtf_format(prn)) {
	pputs(prn, "\\par\n");
    } else {
	pputc(prn, '\n');
    }
}

static void r_squared_message (const MODEL *pmod, PRN *prn)
{
    if (na(pmod->rsq)) {
	return;
    }

    pprintf(prn, "%s.\n\n",    
	    _("R-squared is computed as the square of the correlation "
	      "between observed and\nfitted values of the dependent variable"));
}

static void weighted_stats_message (PRN *prn)
{
    const char *msg = N_("Statistics based on the weighted data");

    if (plain_format(prn)) {
	pprintf(prn, "%s:\n\n", _(msg));
    } else if (tex_format(prn)) {
	pprintf(prn, "\\vspace{1em}%s:\n\n", I_(msg));
    } else { /* RTF */
	pprintf(prn, "\\par \\qc\n%s:\n\n", I_(msg));	
    }
}

static void original_stats_message (PRN *prn)
{
    const char *msg = N_("Statistics based on the original data");

    if (plain_format(prn)) {
	pprintf(prn, "%s:\n\n", _(msg));
    } else if (tex_format(prn)) {
	pprintf(prn, "\\vspace{1em}\n%s:\n\n", I_(msg));
    } else { /* RTF */
	pprintf(prn, "\\par \\qc\n%s:\n\n", I_(msg));
    }
}

static void rho_differenced_stats_message (PRN *prn)
{
    const char *msg = N_("Statistics based on the rho-differenced data");

    if (plain_format(prn)) {    
	pprintf(prn, "%s:\n\n", _(msg));
    } else if (tex_format(prn)) {
	pprintf(prn, "\\vspace{1em}\n%s:\n\n", I_(msg));
    } else { /* RTF */
	pprintf(prn, "\\par \\qc\n%s:\n\n", I_(msg));
    }	
}

static void print_whites_results (const MODEL *pmod, PRN *prn)
{
    if (plain_format(prn)) {
	pprintf(prn, "\n%s: TR^2 = %f,\n", _("Test statistic"), 
		pmod->rsq * pmod->nobs);
	pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n", 
		_("with p-value"), _("Chi-square"), 
		pmod->ncoeff - 1, pmod->rsq * pmod->nobs,
		chisq_cdf_comp(pmod->rsq * pmod->nobs, pmod->ncoeff - 1)); 
    } else if (rtf_format(prn)) { /* FIXME */
	pprintf(prn, "\\par \\ql\n%s: TR{\\super 2} = %f,\n", I_("Test statistic"), 
		pmod->rsq * pmod->nobs);
	pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n", 
		I_("with p-value"), I_("Chi-square"), 
		pmod->ncoeff - 1, pmod->rsq * pmod->nobs,
		chisq_cdf_comp(pmod->rsq * pmod->nobs, pmod->ncoeff - 1)); 
    } else if (tex_format(prn)) {
	pprintf(prn, "\n%s: $TR^2$ = %f,\n", I_("Test statistic"), 
		pmod->rsq * pmod->nobs);
	pprintf(prn, "%s = $P$($\\chi^2(%d)$ > %f) = %f\n\n",
		I_("with p-value"), 
		pmod->ncoeff - 1, pmod->rsq * pmod->nobs,
		chisq_cdf_comp(pmod->rsq * pmod->nobs, pmod->ncoeff - 1)); 
    }
}

static void print_ll (const MODEL *pmod, PRN *prn)
{
    int lldig;

    if (na(pmod->lnL)) {
	return;
    }

    if (pmod->ci == ARMA) {
	lldig = 8;
    } else {
	lldig = XDIGITS(pmod);
    }

    if (plain_format(prn)) {
	double jll;

	pprintf(prn, "  %s = %.*g\n", _("Log-likelihood"), lldig, pmod->lnL);
	jll = gretl_model_get_double(pmod, "jll");
	if (!na(jll)) {
	    char jllstr[64];

	    sprintf(jllstr, _("Log-likelihood for %s"), 
		    (const char *) gretl_model_get_data(pmod, "log-parent"));
	    pprintf(prn, "  (%s = %.*g)\n", jllstr, lldig, jll);
	}
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %.*g\n", I_("Log-likelihood"), GRETL_DIGITS,
		pmod->lnL);
    } else if (tex_format(prn)) {
	char xstr[32];

	tex_dcolumn_double(pmod->lnL, xstr);
	pprintf(prn, "%s & %s \\\\\n", I_("Log-likelihood"), xstr);
    }
}



static char active_decpoint (void)
{
    char test[4];

    sprintf(test, "%.1f", 1.0);
    return test[1];
}

#define fixed_effects_model(m) (m->ci == PANEL && \
                                gretl_model_get_int(m, "fixed-effects"))

#define random_effects_model(m) (m->ci == PANEL && \
                                 gretl_model_get_int(m, "random-effects"))

#define between_model(m) (m->ci == PANEL && \
                          gretl_model_get_int(m, "between"))

/**
 * printmodel:
 * @pmod: pointer to gretl model.
 * @pdinfo: data information struct.
 * @opt: may contain %OPT_O to print covariance matrix, %OPT_S
 * to get a "simple" print (just coefficients and standard
 * errors).
 * @prn: gretl printing struct.
 *
 * Print to @prn the estimates in @pmod plus associated statistics.
 * 
 * Returns: 0 on success, 1 if some of the values to print were %NaN.
 */

int printmodel (MODEL *pmod, const DATAINFO *pdinfo, gretlopt opt, 
		PRN *prn)
{
    int binary = binary_model(pmod);
    int gotnan = 0;

    if (prn == NULL || (opt & OPT_Q)) {
	return 0;
    }

    if (csv_format(prn)) {
	if (active_decpoint() == ',') {
	    gretl_print_set_delim(prn, '\t');
	} else {
	    gretl_print_set_delim(prn, ',');
	}
    }

    if (!plain_format(prn)) {
	model_format_start(prn);
    } else if (pmod->ci != NLS) {
	int iters = gretl_model_get_int(pmod, "iters");

	if (iters > 0) {
	    pprintf(prn, _("Convergence achieved after %d iterations\n"), iters);
	} else {
	    int fncount = gretl_model_get_int(pmod, "fncount");
	    int grcount = gretl_model_get_int(pmod, "grcount");

	    if (fncount > 0) {
		pprintf(prn, _("Function evaluations: %d\n"), fncount);
		pprintf(prn, _("Evaluations of gradient: %d\n"), grcount);
	    }
	}
    }

    print_model_heading(pmod, pdinfo, opt, prn);

    print_coeff_table_start(pmod, prn);
    gotnan = print_coefficients(pmod, pdinfo, prn);

    if (pmod->ci == AR) {
	print_rho_terms(pmod, prn); 
    } else if (pmod->ci == POISSON) {
	print_poisson_offset(pmod, pdinfo, prn);
    }

    print_coeff_table_end(pmod, prn);

    if (pmod->aux == AUX_ARCH || pmod->aux == AUX_ADF || 
	pmod->aux == AUX_RESET || pmod->aux == AUX_DF || 
	pmod->aux == AUX_KPSS) {
	goto close_format;
    }

    if (binary) {
	print_binary_statistics(pmod, pdinfo, prn);
	goto close_format;
    }

    if (pmod->ci == POISSON) {
	print_middle_table_start(prn);
	depvarstats(pmod, prn);
	pseudorsqline(pmod, prn);
	print_ll(pmod, prn);
	info_stats_lines(pmod, prn);
	print_middle_table_end(prn);
	goto close_format;
    }	

    if (pmod->ci == TOBIT) {
	print_middle_table_start(prn);
	depvarstats(pmod, prn);
	print_tobit_stats(pmod, prn);
	print_ll(pmod, prn);
	info_stats_lines(pmod, prn);
	print_middle_table_end(prn);
	goto close_format;
    }
	
    if (ordered_model(pmod)) {
	print_middle_table_start(prn);
	print_ll(pmod, prn);
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
	garch_variance_line(pmod, prn);
	print_ll(pmod, prn);
	info_stats_lines(pmod, prn);
	print_middle_table_end(prn);
	goto close_format;
    } 

    if (pmod->ci == ARBOND || pmod->ci == GMM) {
	print_middle_table_start(prn);
	essline(pmod, prn);
	print_GMM_test_data(pmod, prn);
	print_middle_table_end(prn);
	goto close_format;
    }   

    if (pmod->ci == HECKIT) {
	print_middle_table_start(prn);
	depvarstats(pmod, prn);
	print_heckit_stats(pmod, prn);
	print_middle_table_end(prn);
	goto close_format;
    }
	
    if (pmod->aux == AUX_SYS) {
	print_middle_table_start(prn);
	depvarstats(pmod, prn);
	essline(pmod, prn);
	if (gretl_model_get_int(pmod, "method") == SYS_LIML) {
	    print_liml_equation_data(pmod, prn);
	}
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
	info_stats_lines(pmod, prn);
	if (!plain_format(prn)) {
	    print_middle_table_end(prn);
	}
	goto close_format;
    }

    if (opt & OPT_S) {
	/* --simple-print */
	goto close_format;
    }

    if (pmod->aux == AUX_WHITE) { 
	rsqline(pmod, prn);
	print_whites_results(pmod, prn);
	goto close_format;
    }

    if (pmod->ci == ARMA) {
	print_middle_table_start(prn);
	depvarstats(pmod, prn);
	print_arma_stats(pmod, prn);
	print_ll(pmod, prn);
	info_stats_lines(pmod, prn);
	print_middle_table_end(prn);
	print_arma_roots(pmod, prn);
	goto close_format;
    }

    if (pmod->ci == LOGISTIC) {
	original_stats_message(prn);
    }    

    if (pmod->ci == OLS || pmod->ci == VAR || pmod->ci == TSLS 
	|| pmod->ci == NLS || pmod->ci == MPOLS
	|| (pmod->ci == AR && pmod->arinfo->arlist[0] == 1)
	|| pmod->ci == LOGISTIC || pmod->ci == TOBIT
	|| pmod->ci == PANEL 
	|| (pmod->ci == WLS && gretl_model_get_int(pmod, "wt_dummy"))) {
	print_middle_table_start(prn);
	if (pmod->ci != VAR) {
	    depvarstats(pmod, prn);
	}
	if (essline(pmod, prn)) {
	    print_middle_table_end(prn);
	    goto close_format;
	}

	rsqline(pmod, prn);

	if (random_effects_model(pmod)) {
	    panel_variance_lines(pmod, prn);
	} else if (pmod->ci != NLS && pmod->aux != AUX_VECM) {
	    Fline(pmod, prn);
	}

	if (pmod->ci == PANEL || 
	    (pmod->ci == OLS && dataset_is_panel(pdinfo))) {
	    dwline(pmod, prn);
	} 

	if (dataset_is_time_series(pdinfo)) {
	    if (pmod->ci == OLS || pmod->ci == MPOLS || pmod->ci == VAR ||
		(pmod->ci == WLS && gretl_model_get_int(pmod, "wt_dummy"))) {
		dwline(pmod, prn);
		if (pmod->ci != VAR && pmod->aux != AUX_VECM &&
		    gretl_model_get_int(pmod, "ldepvar")) {
		    dhline(pmod, prn);
		} 
	    }
	    /* FIXME -- check output below */
	    if (pmod->ci == TSLS) {
		dwline(pmod, prn);
	    }
	}

	if (pmod->aux != AUX_VECM) {
	    if (pmod->ci == OLS || pmod->ci == MPOLS ||
		fixed_effects_model(pmod)) {
		print_ll(pmod, prn);
	    }
	    info_stats_lines(pmod, prn);
	}

	print_middle_table_end(prn);

	if (pmod->ci == TSLS && plain_format(prn)) {
	    r_squared_message(pmod, prn);
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

    else if (pmod->ci == MLE) {
	print_middle_table_start(prn);
	print_ll(pmod, prn);
	info_stats_lines(pmod, prn);
	print_middle_table_end(prn);
    }	

    else if (pmod->ci == HSK || pmod->ci == ARCH ||
	     (pmod->ci == WLS && !gretl_model_get_int(pmod, "wt_dummy"))) {

	weighted_stats_message(prn);
	print_middle_table_start(prn);

	if (essline(pmod, prn)) {
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

	if (essline_original(pmod, prn)) {
	    print_middle_table_end(prn);
	    goto close_format;
	}

	print_middle_table_end(prn);
    }

    else if (pmod->ci == AR || pmod->ci == CORC || 
	     pmod->ci == HILU || pmod->ci == PWE) {

	rho_differenced_stats_message(prn);

	print_middle_table_start(prn);
	if (essline(pmod, prn)) {
	    print_middle_table_end(prn);
	    goto close_format;
	}
	rsqline(pmod, prn);
	Fline(pmod, prn);
	dwline(pmod, prn);
	info_stats_lines(pmod, prn);
	print_middle_table_end(prn);
    }

    if (plain_format(prn) && pmod->ci != MLE && pmod->ci != PANEL &&
	pmod->ci != ARMA && pmod->ci != NLS && pmod->ci != GMM &&
	!pmod->aux) {
	pval_max_line(pmod, pdinfo, prn);
    }

 close_format:

    if (opt & OPT_O) {
	outcovmx(pmod, pdinfo, prn);
    }

    if (any_tests(pmod)) {
	print_model_tests(pmod, prn);
    }

    if (!plain_format(prn)) {
	model_format_end(prn);
    }
    
    if (gotnan) {
	pmod->errcode = E_NAN;
    }

    return gotnan;
}

static void print_pval_str (double pval, char *str)
{
    if (pval < .00001) {
	sprintf(str, "<%.5f", 0.00001);
    } else {
	sprintf(str, "%.5f", pval);
    }
}

static void plain_print_pval (double pval, PRN *prn)
{
    char pvstr[16];

    if (pval < .00001) {
	sprintf(pvstr, "<%.5f", 0.00001);
    } else {
	sprintf(pvstr, "%.5f", pval);
    }

    pprintf(prn, "%*s", UTF_WIDTH(pvstr, 10), pvstr);

    if (pval < 0.01) {
	pputs(prn, " ***");
    } else if (pval < 0.05) {
	pputs(prn, " **");
    } else if (pval < 0.10) {
	pputs(prn, " *");
    }	
}

static void plain_print_tval (double tval, PRN *prn)
{
    if (fabs(tval) >= 1000) { 
	/* || t < .001 ? */
	char numstr[9];

	sprintf(numstr, "%#8.2G", tval);
	pprintf(prn, " %8s", numstr);
    } else if (tval <= -100) {
	pprintf(prn, " %7.2f", tval);
    } else {
	pprintf(prn, " %7.3f", tval);
    }
}

static int 
prepare_model_coeff (const MODEL *pmod, const DATAINFO *pdinfo,
		     int i, model_coeff *mc, PRN *prn)
{
    int gotnan = 0;

    model_coeff_init(mc);

    mc->show_pval = !binary_model(pmod) || gretl_model_get_int(pmod, "show-pvals");

    if (tex_format(prn)) {
	make_tex_coeff_name(pmod, pdinfo, i, mc->name);
    } else {
	gretl_model_get_param_name(pmod, pdinfo, i, mc->name);
    }

    if (xna(pmod->coeff[i])) {
	gotnan = 1;
    } else {
	mc->b = pmod->coeff[i];
    }

    if (!xna(pmod->sderr[i])) {
	mc->se = pmod->sderr[i];
    }

    if (!na(mc->b) && !na(mc->se) && mc->se > 0.0) {
	mc->tval = mc->b / mc->se;
	if (xna(mc->tval)) {
	    mc->tval = NADBL;
	}
    }

    if (mc->show_pval && !na(mc->tval)) {
	if (pmod->aux == AUX_ADF || pmod->aux == AUX_DF) {
	    if (i + 2 == gretl_model_get_int(pmod, "dfnum")) {
		mc->pval = gretl_model_get_double(pmod, "dfpval");
	    } 
	    mc->df_pval = 1;
	} else {	
	    mc->pval = coeff_pval(pmod->ci, mc->tval, pmod->dfd);
	}
    }

    if (!gotnan && !mc->show_pval && 
	binary_model(pmod) && 
	pmod->list[i+2] != 0) { 
	double *slopes = gretl_model_get_data(pmod, "slopes");

	if (slopes != NULL) {
	    mc->slope = slopes[i];
	}
    }

    return gotnan;
}

static void plain_print_coeff (const model_coeff *mc, PRN *prn)
{
    pprintf(prn, "  %-15s ", mc->name);

    if (na(mc->b)) {
	pprintf(prn, "%*s", UTF_WIDTH(_("undefined"), 17), _("undefined"));
    } else {
	gretl_print_value(mc->b, prn);
    }

    /* get out if std error is undefined */
    if (na(mc->se)) {
	pprintf(prn, "%*s\n", UTF_WIDTH(_("undefined"), 16), _("undefined"));
	return;
    }

    gretl_print_value(mc->se, prn); 

    if (!na(mc->tval)) {
	plain_print_tval(mc->tval, prn);
    }

    if (!na(mc->slope)) {
	/* slope for binary models */
	gretl_print_value(mc->slope, prn);
    } else if (mc->show_pval) {
	if (na(mc->pval)) {
	    if (!mc->df_pval) {
		pprintf(prn, "     %*s", UTF_WIDTH(_("undefined"), 10), _("undefined"));
	    }
	} else {
	    plain_print_pval(mc->pval, prn);
	} 
    }

    pputc(prn, '\n');
}

static void rtf_print_double (double xx, PRN *prn)
{
    char numstr[32];

    xx = screen_zero(xx);
    sprintf(numstr, "%.*g", GRETL_DIGITS, xx);
    gretl_fix_exponent(numstr);
    pprintf(prn, " \\qc %s\\cell", numstr);
}

static void rtf_print_coeff (const model_coeff *mc, PRN *prn)
{
    pputs(prn, RTF_COEFF_ROW);
    pprintf(prn, "\\ql %s\\cell", mc->name);

    if (na(mc->b)) {
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
    } else {
	rtf_print_double(mc->b, prn);
    }

    if (na(mc->se)) {
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	goto rtf_coeff_getout;
    } 

    rtf_print_double(mc->se, prn); 

    if (!na(mc->tval)) {
	pprintf(prn, " \\qc %.4f\\cell", mc->tval);
    } else {
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
    }

    if (!na(mc->slope)) {
	rtf_print_double(mc->slope, prn);
    } else if (mc->show_pval) {
	if (na(mc->pval)) {
	    if (mc->df_pval) {
		pprintf(prn, " \\qc %s\\cell", I_("unknown"));
	    } else {
		pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	    }
	} else {
	    char pvalstr[16];

	    print_pval_str(mc->pval, pvalstr);
	    pprintf(prn, " \\qc %s\\cell", pvalstr);

	    if (mc->pval < 0.01) {
		pputs(prn, " \\ql ***\\cell");
	    } else if (mc->pval < 0.05) { 
		pputs(prn, " \\ql **\\cell");
	    } else if (mc->pval < 0.10) {
		pputs(prn, " \\ql *\\cell");
	    } else {
		pputs(prn, " \\ql \\cell");
	    }
	} 
    }

 rtf_coeff_getout:

    pputs(prn, " \\intbl \\row\n");
}

static void csv_print_coeff (const model_coeff *mc, PRN *prn)
{
    char d = prn_delim(prn);

    pprintf(prn, "\"%s\"", mc->name);

    if (na(mc->b)) {
	pprintf(prn, "%c\"%s\"", d, I_("undefined"));
    } else {
	pprintf(prn, "%c%.15g", d, mc->b);
    }
    
    /* get out if std error is undefined */
    if (na(mc->se)) {
	pprintf(prn, "%c\"%s\"\n", d, I_("undefined"));
	return;
    }

    pprintf(prn, "%c%.15g", d, mc->se);

    if (!na(mc->tval)) {
	pprintf(prn, "%c%.15g", d, mc->tval);
    } else {
	pprintf(prn, "%c\"%s\"\n", d, I_("undefined"));
    }

    if (!na(mc->slope)) {
	/* slope for binary models */
	pprintf(prn, "%c%.15g", d, mc->slope);
    } else if (mc->show_pval) {
	if (na(mc->pval)) {
	    if (mc->df_pval) {
		pprintf(prn, "%c\"%s\"\n", d, I_("unknown"));
	    } else {
		pprintf(prn, "%c\"%s\"\n", d, I_("undefined"));
	    }
	} else {
	    pprintf(prn, "%c%.15g", d, mc->pval);
	} 
    }

    pputc(prn, '\n');
}

static void print_mp_coeff (const model_coeff *mc, PRN *prn)
{
    pprintf(prn, "  %-*s", VNAMELEN - 1, mc->name);

    gretl_print_fullwidth_double(mc->b, GRETL_MP_DIGITS, prn);
    gretl_print_fullwidth_double(mc->se, GRETL_MP_DIGITS, prn); 

    pputc(prn, '\n');
}

void print_coeff (const model_coeff *mc, PRN *prn)
{
    if (plain_format(prn)) {
	plain_print_coeff(mc, prn);
    } else if (rtf_format(prn)) {
	rtf_print_coeff(mc, prn);
    } else if (tex_format(prn)) {
	tex_print_coeff(mc, prn);
    } else if (csv_format(prn)) {
	csv_print_coeff(mc, prn);
    }
}

void print_arch_coeffs (const double *a, const double *se,
			int T, int order, PRN *prn,
			int aux)
{
    model_coeff mc;
    int i;

    if (aux) {
	pprintf(prn, "\n%s %d\n\n", _("Test for ARCH of order"),
		order);
	pputs(prn, _("      PARAMETER       ESTIMATE          STDERROR"
		     "      T STAT   P-VALUE\n\n"));
    } else {
	gretl_prn_newline(prn);
    }

    for (i=0; i<=order; i++) {
	model_coeff_init(&mc);
	mc.b = a[i];
	mc.se = se[i];
	mc.tval = a[i] / se[i];
	mc.pval = t_pvalue_2(mc.tval, T - (order + 1));

	if (tex_format(prn)) {
	    sprintf(mc.name, "$\\alpha_%d$", i);
	} else {
	    sprintf(mc.name, "alpha(%d)", i);
	}

	print_coeff(&mc, prn);
    }
}

static void print_coeff_separator (const char *s, PRN *prn)
{
    if (tex_format(prn)) {
	/* FIXME */
	pputs(prn, "\\\\ \n");
    } else {
	if (s != NULL) {
	    pputc(prn, '\n');
	    print_centered((rtf_format(prn))? I_(s) : _(s), 78, prn);
	    pputc(prn, '\n');
	}
	pputc(prn, '\n');
    }
}

static int 
print_coefficients (const MODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    const char *sepstr = NULL;
    int seppos = -1;
    model_coeff mc;
    int nc = pmod->ncoeff;
    int i, err = 0, gotnan = 0;

    if (pmod->ci == PANEL) {
	nc = pmod->list[0] - 1;
    }

    gretl_model_get_coeff_separator(pmod, &sepstr, &seppos);
    if (seppos == -1 && pmod->ci == GARCH) {
	seppos = pmod->list[0] - 4;
    }

    for (i=0; i<nc; i++) {

	err = prepare_model_coeff(pmod, pdinfo, i, &mc, prn);
	if (err) gotnan = 1;

	if (i == seppos) {
	    print_coeff_separator(sepstr, prn);
	}

	if (plain_format(prn) && pmod->ci == MPOLS) {
	    print_mp_coeff(&mc, prn);
	    continue;
	}	    

	print_coeff(&mc, prn);
    }

    if (pmod->ci == ARCH) {
	double *a = gretl_model_get_data(pmod, "arch_coeff");
	double *se = gretl_model_get_data(pmod, "arch_sderr");
	int order = gretl_model_get_int(pmod, "arch_order");

	if (a != NULL && se != NULL && order > 0) {
	    print_arch_coeffs(a, se, pmod->nobs, order, prn, 0);
	}
    }

    return gotnan;
} 

static void print_rho (const ARINFO *arinfo, int c, int dfd, PRN *prn)
{
    char ustr[32];
    double xx = arinfo->rho[c] / arinfo->sderr[c];

    if (plain_format(prn)) {
	char pvalstr[16];
	double pval;

	sprintf(ustr, "u_%d", arinfo->arlist[c+1]);
	pprintf(prn, "%14s", ustr); 
	bufspace(3, prn);
	gretl_print_value (arinfo->rho[c], prn);
	bufspace(2, prn);
	gretl_print_value (arinfo->sderr[c], prn); 
	pval = t_pvalue_2(xx, dfd);
	pprintf(prn, " %7.3f ", xx, pval);
	print_pval_str(pval, pvalstr);
	pprintf(prn, "%*s\n", UTF_WIDTH(pvalstr, 12), pvalstr);
    } else if (tex_format(prn)) {
	char coeff[32], sderr[32];

	tex_dcolumn_double(arinfo->rho[c], coeff);
	tex_dcolumn_double(arinfo->sderr[c], sderr);

	sprintf(ustr, "$\\hat{u}_{t-%d}$", arinfo->arlist[c+1]);

	pprintf(prn, "%s &\n"
		"  %s &\n"
		"    %s &\n"
		"      %.4f &\n"
		"        %.4f \\\\\n",  
		ustr,
		coeff,
		sderr,
		arinfo->rho[c] / arinfo->sderr[c],
		t_pvalue_2(xx, dfd));
    } else if (rtf_format(prn)) {
	char pvalstr[16];
	double pval;

	pputs(prn, RTF_COEFF_ROW);
	pprintf(prn, "\\ql u(-%d)\\cell", arinfo->arlist[c+1]);
	rtf_print_double(arinfo->rho[c], prn);
	rtf_print_double(arinfo->sderr[c], prn);
	pprintf(prn, " \\qc %.4f\\cell", xx);
	pval = t_pvalue_2(xx, dfd);
	print_pval_str(pval, pvalstr);
	pprintf(prn, " \\qc %s\\cell", pvalstr);
	if (pval < 0.01) {
	    pputs(prn, " \\ql ***\\cell");
	} else if (pval < 0.05) { 
	    pputs(prn, " \\ql **\\cell");
	} else if (pval < 0.10) {
	    pputs(prn, " \\ql *\\cell");
	} else {
	    pputs(prn, " \\ql \\cell");
	}
	pputs(prn, " \\intbl \\row\n");
    }
}

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

    if (plain_format(prn)) {
	pprintf(prn, "\n%s:\n\n", _("Estimates of the AR coefficients"));
    }

    if (pmod->arinfo->arlist[0] > 1) {
	dfd = pmod->dfd + (pmod->ncoeff - pmod->arinfo->arlist[0]);
    } else {
	dfd = pmod->dfd;
    }

    for (i=1; i<=pmod->arinfo->arlist[0]; i++) {
	print_rho(pmod->arinfo, i - 1, dfd, prn);
	xx += pmod->arinfo->rho[i-1]; 
    }

    if (pmod->arinfo->arlist[0] > 1 && plain_format(prn)) {
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

static void print_root (double rx, double ix, double mod, double fr,
			int i, int hline, PRN *prn)
{
     if (plain_format(prn)) {
	 pprintf(prn, root_fmt, _("Root"), i, rx, ix, mod, fr);
     } else if (tex_format(prn)) {
	 pprintf(prn, "& %s & %d & $%.4f$ & $%.4f$ & $%.4f$ & $%.4f$ \\\\ ",
		 I_("Root"), i, rx, ix, mod, fr);
	 if (hline) {
	     pputs(prn, "\\hline\n");
	 } else {
	     pputc(prn, '\n');
	 }
     } else if (rtf_format(prn)) {
	 pputs(prn, RTF_ROOT_ROW);
	 pprintf(prn, "\\ql \\cell \\ql %s %d \\cell"
		 " \\qr %.4f\\cell"
		 " \\qr %.4f\\cell"
		 " \\qr %.4f\\cell"
		 " \\qr %.4f\\cell \\intbl \\row\n",
		 I_("Root"), i, rx, ix, mod, fr);
     }
}

static void root_start (const char *tag, PRN *prn)
{
    if (plain_format(prn)) {
	pprintf(prn, "  %s\n", _(tag));
    } else if (tex_format(prn)) {
	pprintf(prn, "%s \\\\ \n", I_(tag));
    } else if (rtf_format(prn)) {
	pputs(prn, RTF_ROOT_ROW);
	pprintf(prn, "\\ql %s\\cell\\ql \\cell\\ql \\cell\\ql \\cell\\ql \\cell"
		"\\ql\\cell \\intbl \\row\n", I_(tag));
    }
}

static void print_arma_roots (const MODEL *pmod, PRN *prn)
{
    cmplx *roots = gretl_model_get_data(pmod, "roots");

    if (roots != NULL) {
	int p = arma_model_nonseasonal_AR_order(pmod);
	int q = arma_model_nonseasonal_MA_order(pmod);
	int P = gretl_model_get_int(pmod, "arma_P");
	int Q = gretl_model_get_int(pmod, "arma_Q");
	int i, k, hline;
	double mod, fr;

	if (plain_format(prn)) {
	    pprintf(prn, "\n%s\n%s\n", _(roots_hdr), roots_sep);
	} else if (tex_format(prn)) {
	    pputs(prn, "\n\\vspace{1em}\n\n");
	    pputs(prn, "\\begin{tabular}{llrrrrr}\n");
	    pprintf(prn, "& & & %s & %s & %s & %s \\\\ \\hline\n", 
		    I_("Real"), I_("Imaginary"), I_("Modulus"), I_("Frequency"));
	} else if (rtf_format(prn)) {
	    pputs(prn, "\n\\par\n{" RTF_ROOT_ROW);
	    pprintf(prn, "\\qr \\cell \\qc \\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell \\intbl \\row\n",
		    I_("Real"), I_("Imaginary"), I_("Modulus"), I_("Frequency"));
	}

	if (p > 0) {
	    k = 1;
	    root_start(N_("AR"), prn);
	    for (i=0; i<p; i++) {
		if (roots[i].i != 0) {
		    mod = roots[i].r * roots[i].r + roots[i].i * roots[i].i;
		    mod = sqrt(mod);
		} else {
		    mod = fabs(roots[i].r);
		}
		fr = atan2(roots[i].i, roots[i].r) / M_2PI;
		if (i == p - 1 && q == 0 && P == 0 && Q == 0) {
		    hline = 1;
		} else {
		    hline = 0;
		}
		print_root(roots[i].r, roots[i].i, mod, fr, k++, hline, prn);
	    }
	}

	if (P > 0) {
	    k = 1;
	    root_start(N_("AR (seasonal)"), prn);
	    for (i=p; i<p+P; i++) {
		if (roots[i].i != 0) {
		    mod = roots[i].r * roots[i].r + roots[i].i * roots[i].i;
		    mod = sqrt(mod);
		} else {
		    mod = fabs(roots[i].r);
		}
		fr = atan2(roots[i].i, roots[i].r) / M_2PI;
		if (i == p + P - 1 && q == 0 && Q == 0) {
		    hline = 1;
		} else {
		    hline = 0;
		}
		print_root(roots[i].r, roots[i].i, mod, fr, k++, hline, prn);
	    }
	}

	if (q > 0) {
	    k = 1;
	    root_start(N_("MA"), prn);
	    for (i=p+P; i<p+P+q; i++) {
		if (roots[i].i != 0) {
		    mod = roots[i].r * roots[i].r + roots[i].i * roots[i].i;
		    mod = sqrt(mod);
		} else {
		    mod = fabs(roots[i].r);
		}
		fr = atan2(roots[i].i, roots[i].r) / M_2PI;
		if (i == p + P + q - 1 && Q == 0) {
		    hline = 1;
		} else {
		    hline = 0;
		}
		print_root(roots[i].r, roots[i].i, mod, fr, k++, hline, prn);
	    }
	}

	if (Q > 0) {
	    k = 1;
	    root_start(N_("MA (seasonal)"), prn);
	    for (i=p+P+q; i<p+P+q+Q; i++) {
		if (roots[i].i != 0) {
		    mod = roots[i].r * roots[i].r + roots[i].i * roots[i].i;
		    mod = sqrt(mod);
		} else {
		    mod = fabs(roots[i].r);
		}
		fr = atan2(roots[i].i, roots[i].r) / M_2PI;
		if (i == p + P + q + Q - 1) {
		    hline = 1;
		} else {
		    hline = 0;
		}
		print_root(roots[i].r, roots[i].i, mod, fr, k++, hline, prn);
	    }
	}

	if (plain_format(prn)) {
	    pprintf(prn, "%s\n\n", roots_sep);
	} else if (tex_format(prn)) {
	    pputs(prn, "\\end{tabular}\n");
	} else if (rtf_format(prn)) {
	    pputs(prn, "}\n");
	}
    }
}

static void print_tobit_stats (const MODEL *pmod, PRN *prn)
{
    int cenc = gretl_model_get_int(pmod, "censobs");
    double cenpc = 100.0 * cenc / pmod->nobs;

    if (plain_format(prn)) {
	pprintf(prn, "  %s: %d (%.1f%%)\n", _("Censored observations"), cenc, cenpc);
	pprintf(prn, "  %s = %.*g\n", _("sigma"), GRETL_DIGITS, pmod->sigma);
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s: %d (%.1f%%)\n", I_("Censored observations"), cenc, cenpc);
	pprintf(prn, RTFTAB "%s = %g\n", I_("sigma"), pmod->sigma);
    } else if (tex_format(prn)) {
	char xstr[32];

	pprintf(prn, "%s & \\multicolumn{1}{r}{%.1f\\%%} \\\\\n", 
		I_("Censored observations"), cenpc);
	tex_dcolumn_double(pmod->sigma, xstr);
	pprintf(prn, "$\\hat{\\sigma}$ & %s \\\\\n", xstr);
    }
}

static void print_heckit_stats (const MODEL *pmod, PRN *prn)
{
    int totobs = gretl_model_get_int(pmod, "totobs");
    int cenobs = totobs - pmod->nobs;
    double cenpc = 100.0 * cenobs / totobs;

    if (plain_format(prn)) {
	pprintf(prn, "  %s: %d\n", _("Total observations"), totobs);
	pprintf(prn, "  %s: %d (%.1f%%)\n", _("Censored observations"), cenobs, cenpc);
	pprintf(prn, "  %s = %.*g\n", _("sigma"), GRETL_DIGITS, pmod->sigma);
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s: %d\n", I_("Total observations"), totobs);
	pprintf(prn, RTFTAB "%s: %d (%.1f%%)\n", I_("Censored observations"), cenobs, cenpc);
	pprintf(prn, RTFTAB "%s = %g\n", I_("sigma"), pmod->sigma);
    } else if (tex_format(prn)) {
	char xstr[32];
	
	pprintf(prn, "%s & \\multicolumn{1}{r}{%d} \\\\\n", 
		I_("Total observations"), totobs);
	pprintf(prn, "%s & \\multicolumn{1}{r}{%.1f\\%%} \\\\\n", 
		I_("Censored observations"), cenpc);
	tex_dcolumn_double(pmod->sigma, xstr);
	pprintf(prn, "$\\hat{\\sigma}$ & %s \\\\\n", xstr);
    }
}

static void 
print_poisson_offset (const MODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    int offvar = gretl_model_get_int(pmod, "offset_var");

    if (offvar > 0) {
	char namestr[24];

	sprintf(namestr, "log(%s)", pdinfo->varname[offvar]);
	if (plain_format(prn)) {
	    pprintf(prn, "\n %13s         1.0\n", namestr);
	} else if (rtf_format(prn)) {
	    pputs(prn, RTF_COEFF_ROW);
	    pprintf(prn, "\\ql %s\\cell\\qc 1.0\\cell", namestr);
	    pputs(prn, "\\qc \\cell\\qc \\cell \\qc \\cell \\intbl \\row\n");
	} else if (tex_format(prn)) {
	    char tmp[32];

	    tex_escape(tmp, namestr);
	    pprintf(prn, "{\\rm %s} & \\multicolumn{1}{c}{1.0} \\\\\n", tmp);
	}
    }
}

static void print_arma_stats (const MODEL *pmod, PRN *prn)
{
    double mu = gretl_model_get_double(pmod, "mean_error");
    double var = pmod->sigma * pmod->sigma;

    if (plain_format(prn)) {
	pprintf(prn, "  %s = %.*g\n", _("Mean of innovations"), 
		GRETL_DIGITS, mu);
	pprintf(prn, "  %s = %.*g\n", _("Variance of innovations"), 
		GRETL_DIGITS, var);
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_("Mean of innovations"), mu);
	pprintf(prn, RTFTAB "%s = %g\n", I_("Variance of innovations"), var);
    } else if (tex_format(prn)) {
	char xstr[32];

	tex_dcolumn_double(mu, xstr);
	pprintf(prn, "%s & %s \\\\\n", I_("Mean of innovations"), xstr);
	tex_dcolumn_double(var, xstr);
	pprintf(prn, "%s & %s \\\\\n", I_("Variance of innovations"), xstr);
    }
}

static void plain_print_act_pred (const int *ap, PRN *prn)
{
    int leftlen;
    int numwidth = 1;
    int i, bign = 0;

    for (i=0; i<4; i++) {
	if (ap[i] > bign) {
	    bign = ap[i];
	}
    }

    while (bign /= 10) {
	numwidth++;
    }

    leftlen = strlen(_("Actual")) + 3; /* utflen */

    bufspace(leftlen + 2, prn);
    pputs(prn, _("Predicted"));
    pputc(prn, '\n');
    bufspace(leftlen + 3, prn);
    pprintf(prn, "%*d   %*d\n", numwidth, 0, numwidth, 1);
    bufspace(2, prn);
    pputs(prn, _("Actual"));
    pprintf(prn, " 0  %*d   %*d\n", numwidth, ap[0], numwidth, ap[1]);
    bufspace(leftlen, prn);
    pprintf(prn, "1  %*d   %*d\n", numwidth, ap[2], numwidth, ap[3]);
    pputc(prn, '\n');
}

static void print_binary_statistics (const MODEL *pmod, 
				     const DATAINFO *pdinfo,
				     PRN *prn)
{
    int slopes = !gretl_model_get_int(pmod, "show-pvals");
    double model_chisq = gretl_model_get_double(pmod, "chisq");
    const double *crit = pmod->criterion;
    const int *act_pred;
    int correct = -1;
    double pc_correct = NADBL;
    int i;

    act_pred = gretl_model_get_data(pmod, "discrete_act_pred");
    if (act_pred != NULL) {
	correct = act_pred[0] + act_pred[3];
	pc_correct = 100 * (double) correct / pmod->nobs;
    } 

    if (plain_format(prn)) {
	pprintf(prn, "  %s %s = %.3f\n", _("Mean of"), 
		pdinfo->varname[pmod->list[1]], pmod->ybar);
	if (correct >= 0) {
	    pprintf(prn, "  %s = %d (%.1f%%)\n", 
		    _("Number of cases 'correctly predicted'"), 
		    correct, pc_correct);
	}
	pprintf(prn, "  f(beta'x) %s = %.3f\n", _("at mean of independent vars"), 
		pmod->sdy);
	pseudorsqline(pmod, prn);
	pprintf(prn, "  %s = %g\n", _("Log-likelihood"), pmod->lnL);
	if (pmod->aux != AUX_OMIT && pmod->aux != AUX_ADD) {
	    i = pmod->ncoeff - 1;
	    pprintf(prn, "  %s: %s(%d) = %g (%s %f)\n",
		    _("Likelihood ratio test"), _("Chi-square"), 
		    i, model_chisq, _("p-value"), 
		    chisq_cdf_comp(model_chisq, i));
	}
	pprintf(prn, "  %s (%s) = %g\n", _(aic_str), _(aic_abbrev),
		crit[C_AIC]);
	pprintf(prn, "  %s (%s) = %g\n", _(bic_str), _(bic_abbrev),
		crit[C_BIC]);
	pprintf(prn, "  %s (%s) = %g\n", _(hqc_str), _(hqc_abbrev),
		crit[C_HQC]);

	pputc(prn, '\n');
	if (act_pred != NULL) {
	    plain_print_act_pred(act_pred, prn);
	}
    } else if (rtf_format(prn)) {
	pputc(prn, '\n');
	if (slopes) {
	    pprintf(prn, "\\par {\\super *}%s\n", I_("Evaluated at the mean"));
	}
	pprintf(prn, "\\par %s %s = %.3f\n", I_("Mean of"), 
		pdinfo->varname[pmod->list[1]], pmod->ybar);
	if (correct >= 0) {
	    pprintf(prn, "\\par %s = %d (%.1f%%)\n", 
		    I_("Number of cases 'correctly predicted'"), 
		    correct, pc_correct);
	}
	pprintf(prn, "\\par f(beta'x) %s = %.3f\n", I_("at mean of independent vars"), 
		pmod->sdy);
	pseudorsqline(pmod, prn);
	pprintf(prn, "\\par %s = %g\n", I_("Log-likelihood"), pmod->lnL);
	if (pmod->aux != AUX_OMIT && pmod->aux != AUX_ADD) {
	    i = pmod->ncoeff - 1;
	    pprintf(prn, "\\par %s: %s(%d) = %g (%s %f)\n",
		    I_("Likelihood ratio test"), I_("Chi-square"), 
		    i, model_chisq, I_("p-value"), 
		    chisq_cdf_comp(model_chisq, i));
	}
	pprintf(prn, "\\par %s (%s) = %g\n", I_(aic_str), I_(aic_abbrev),
		crit[C_AIC]);
	pprintf(prn, "\\par %s (%s) = %g\\par\n", I_(bic_str), I_(bic_abbrev),
		crit[C_BIC]);
	pprintf(prn, "\\par %s (%s) = %g\\par\n", I_(hqc_str), I_(hqc_abbrev),
		crit[C_HQC]);
	pputc(prn, '\n');
    } else if (tex_format(prn)) {
	char lnlstr[16];

	if (slopes) {
	    pprintf(prn, "\\begin{center}\n$^*$%s\n\\end{center}\n", 
		    I_("Evaluated at the mean"));
	}

	pputs(prn, "\\vspace{1em}\n\\begin{raggedright}\n");
	pprintf(prn, "%s %s = %.3f\\\\\n", I_("Mean of"), 
		pdinfo->varname[pmod->list[1]], pmod->ybar);
	if (correct >= 0) {
	    pprintf(prn, "%s = %d (%.1f percent)\\\\\n", 
		    I_("Number of cases `correctly predicted'"), 
		    correct, pc_correct);
	}
	pseudorsqline(pmod, prn);
	pprintf(prn, "$f(\\beta'x)$ %s = %.3f\\\\\n", I_("at mean of independent vars"), 
		pmod->sdy);
	tex_float_str(pmod->lnL, lnlstr);
	pprintf(prn, "%s = %s\\\\\n", I_("Log-likelihood"), lnlstr);
	if (pmod->aux != AUX_OMIT && pmod->aux != AUX_ADD) {
	    i = pmod->ncoeff - 1;
	    pprintf(prn, "%s: $\\chi^2_{%d}$ = %.3f (%s %f)\\\\\n",
		    I_("Likelihood ratio test"), 
		    i, model_chisq, I_("p-value"), 
		    chisq_cdf_comp(model_chisq, i));
	}
	pprintf(prn, "%s (%s) = %g\\\\\n", I_(aic_str), I_(aic_abbrev),
		crit[C_AIC]);
	pprintf(prn, "%s (%s) = %g\\\\\n", I_(bic_str), I_(bic_abbrev),
		crit[C_BIC]);
	pprintf(prn, "%s (%s) = %g\\\\\n", I_(tex_hqc_str), I_(hqc_abbrev),
		crit[C_HQC]);
	pputs(prn, "\\end{raggedright}\n");
    }
}


int ols_print_anova (const MODEL *pmod, PRN *prn)
{
    double mst, msr, mse, rss;

    if (pmod->ci != OLS || !pmod->ifc ||
	na(pmod->ess) || na(pmod->tss)) {
	return 1;
    }

    pprintf(prn, "%s:\n\n", _("Analysis of Variance"));

    rss = pmod->tss - pmod->ess;

    pprintf(prn, "%35s %8s %16s\n\n", _("Sum of squares"), _("df"), _("Mean square"));

    msr = rss / pmod->dfn;
    pprintf(prn, "  %-16s %16g %8d %16g\n", _("Regression"), rss, pmod->dfn, msr);

    mse = pmod->ess / pmod->dfd;
    pprintf(prn, "  %-16s %16g %8d %16g\n", _("Residual"), pmod->ess, pmod->dfd, mse);

    mst = pmod->tss / pmod->dfd;
    pprintf(prn, "  %-16s %16g %8d %16g\n", _("Total"), pmod->tss, pmod->nobs - 1, mst);

    pprintf(prn, "\n  R^2 = %g / %g = %.6f\n", rss, pmod->tss, rss / pmod->tss);

    if (pmod->ess == 0) {
	pprintf(prn, "  F(%d, %d) = %g / %g (%s)\n\n", pmod->dfn, pmod->dfd, 
		msr, mse, _("undefined"));
    } else {
	pprintf(prn, "  F(%d, %d) = %g / %g = %g\n\n", 
		pmod->dfn, pmod->dfd, msr, mse, msr / mse);
    }

    return 0;
}
