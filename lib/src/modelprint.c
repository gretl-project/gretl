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

#include "libgretl.h"
#include "libset.h"
#include "system.h"
#include "texprint.h"
#include "usermat.h"
#include "gretl_string_table.h"

#include <glib.h>

static int 
plain_print_coefficients (const MODEL *pmod, const DATAINFO *pdinfo, PRN *prn);
static int 
alt_print_coefficients (const MODEL *pmod, const DATAINFO *pdinfo, PRN *prn);
static void alt_print_rho_terms (const MODEL *pmod, PRN *prn);
static void print_binary_statistics (const MODEL *pmod, 
				     const DATAINFO *pdinfo,
				     PRN *prn);
static void print_arma_roots (const MODEL *pmod, PRN *prn);
static void print_heckit_stats (const MODEL *pmod, PRN *prn);

#define RTFTAB "\\par \\ql \\tab "

#define XDIGITS(m) (((m)->ci == MPOLS)? GRETL_MP_DIGITS : GRETL_DIGITS)

#define FDIGITS(m) (((m)->ci == MPOLS)? GRETL_MP_DIGITS : 5)

#define ordered_model(m) ((m->ci == LOGIT || m->ci == PROBIT) && \
                           gretl_model_get_int(m, "ordered"))

#define binary_model(m) ((m->ci == LOGIT || m->ci == PROBIT) && \
                         !gretl_model_get_int(m, "ordered"))

#define nonlin_model(m) (m->ci == NLS || m->ci == MLE || m->ci == GMM)

#define liml_equation(m) (gretl_model_get_int(m, "method") == SYS_METHOD_LIML)

void model_coeff_init (model_coeff *mc)
{
    mc->b = NADBL;
    mc->se = NADBL;
    mc->tval = NADBL;
    mc->pval = NADBL;
    mc->slope = NADBL; 
    mc->lo = mc->hi = NADBL;
    mc->show_pval = 1;
    mc->df_pval = 0;
    mc->multi = 0;
    mc->name[0] = '\0';
}

static int char_len (const char *s)
{
    if (g_utf8_validate(s, -1, NULL)) {
	return g_utf8_strlen(s, -1);
    } else {
	return strlen(s);
    }
}

static void plain_print_double (char *s, int d, double x, PRN *prn)
{
    if (x < 0 && gretl_print_supports_utf(prn)) {
	char tmp[32];

	*s = '\0';
	strcat(s, "−"); /* U+2212: minus */
	sprintf(tmp, "%.*g", d, -x);
	strcat(s, tmp);
    } else {
	sprintf(s, "%.*g", d, x);
    }
}

/* for use when printing a user-defined model */

static void print_model_stats_table (const double *stats, 
				     const char **names, 
				     int ns, PRN *prn)
{
    char tmp1[32], tmp2[32];
    int i;

    if (plain_format(prn)) {
	pputc(prn, '\n');
    } else if (tex_format(prn)) {
	pputs(prn, "\\medskip\n\n");
	pputs(prn, "\\begin{tabular}{lr@{.}l\n");
    }

    for (i=0; i<ns; i++) {
	if (plain_format(prn)) {
	    plain_print_double(tmp1, GRETL_DIGITS, stats[i], prn);
	    pprintf(prn, "  %s = %s\n", names[i], tmp1);
	} else if (tex_format(prn)) {
	    tex_escape_special(tmp1, names[i]);
	    tex_rl_double(stats[i], tmp2);
	    pprintf(prn, "%s & %s \\\\\n", tmp1, tmp2);
	} else if (rtf_format(prn)) {
	    pprintf(prn, RTFTAB "%s = %g\n", names[i], stats[i]);
	} else if (csv_format(prn)) {
	    pprintf(prn, "\"%s\"%c%.15g\n", names[i], prn_delim(prn), stats[i]);
	}
    }	

    if (tex_format(prn)) {
	pputs(prn, "\\end{tabular}");
    }
}

static void ensure_vsep (PRN *prn)
{
    if (tex_format(prn)) {
	pputs(prn, "\n\\vspace{1ex}\n");
    } else if (rtf_format(prn)) {
	pputs(prn, "\\par\n");
    } 
}

static void garch_variance_line (const MODEL *pmod, PRN *prn)
{
    const char *varstr = N_("Unconditional error variance");
    double v = pmod->sigma * pmod->sigma;

    ensure_vsep(prn);

    if (plain_format(prn)) {  
	pprintf(prn, "%s = %.*g\n\n", _(varstr), GRETL_DIGITS, v);
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_(varstr), v);
    } else if (tex_format(prn)) {
	pprintf(prn, "%s = %g \\\\\n", I_(varstr), v);
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", I_(varstr), prn_delim(prn), v);
    }
}

static void print_intreg_info (const MODEL *pmod, PRN *prn)
{
    const char *nstrs[] = {
	N_("Left-unbounded observations"),
	N_("Right-unbounded observations"),
	N_("Bounded observations"),
	N_("Point observations")
    };
    int nl = gretl_model_get_int(pmod, "n_left");
    int nr = gretl_model_get_int(pmod, "n_right");
    int nb = gretl_model_get_int(pmod, "n_both");
    int np = gretl_model_get_int(pmod, "n_point");

    ensure_vsep(prn);

    if (plain_format(prn)) {  
	pprintf(prn, "%s = %.*g\n", _("sigma"), GRETL_DIGITS, pmod->sigma);
	pprintf(prn, "%s: %d\n", _(nstrs[0]), nl);
	pprintf(prn, "%s: %d\n", _(nstrs[1]), nr);
	pprintf(prn, "%s: %d\n", _(nstrs[2]), nb);
	pprintf(prn, "%s: %d\n", _(nstrs[3]), np);
	pputc(prn, '\n');
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g\n", I_("sigma"), pmod->sigma);
	pprintf(prn, RTFTAB "%s: %d\n", I_(nstrs[0]), nl);
	pprintf(prn, RTFTAB "%s: %d\n", I_(nstrs[1]), nr);
	pprintf(prn, RTFTAB "%s: %d\n", I_(nstrs[2]), nb);
	pprintf(prn, RTFTAB "%s: %d\n", I_(nstrs[3]), np);
    } else if (tex_format(prn)) {
	pprintf(prn, "$\\hat{\\sigma}$ = %g \\\\\n", pmod->sigma);
	pprintf(prn, "%s: %d \\\\\n", I_(nstrs[0]), nl);
	pprintf(prn, "%s: %d \\\\\n", I_(nstrs[1]), nr);
	pprintf(prn, "%s: %d \\\\\n", I_(nstrs[2]), nb);
	pprintf(prn, "%s: %d \\\\\n", I_(nstrs[3]), np);
    } else if (csv_format(prn)) {
	int d = prn_delim(prn);

	pprintf(prn, "%s%c%.15g\n", I_("sigma"), d, pmod->sigma);
	pprintf(prn, "\"%s\"%c%d\n", I_(nstrs[0]), d, nl);
	pprintf(prn, "\"%s\"%c%d\n", I_(nstrs[1]), d, nr);
	pprintf(prn, "\"%s\"%c%d\n", I_(nstrs[2]), d, nb);
	pprintf(prn, "\"%s\"%c%d\n", I_(nstrs[3]), d, np);
    }
}

static void rsqline (const MODEL *pmod, PRN *prn)
{
    if (!na(pmod->rsq) && plain_format(prn)) {
	pprintf(prn, "  %s = %f\n", _("Unadjusted R-squared"), pmod->rsq);
    }
}

static void rssline (const MODEL *pmod, PRN *prn)
{
    if (!na(pmod->ess) && !na(pmod->tss)) {
	if (plain_format(prn)) { 
	    pprintf(prn, "  %s = %.*g\n", _("Explained sum of squares"), 
		XDIGITS(pmod), pmod->tss - pmod->ess);
	} 
    }
}

static void print_liml_equation_data (const MODEL *pmod, PRN *prn)
{
    double lmin = gretl_model_get_double(pmod, "lmin");
    int idf = gretl_model_get_int(pmod, "idf");

    if (!na(lmin)) {
	ensure_vsep(prn);
	if (idf > 0) {
	    double X2 = pmod->nobs * log(lmin);
	    double pv = chisq_cdf_comp(idf, X2);

	    if (tex_format(prn)) {
		pprintf(prn, "%s: $\\chi^2(%d)$ = %g [%.4f] \\\\\n", 
			I_("LR over-identification test"), idf, X2, pv);
	    } else if (rtf_format(prn)) {
		pprintf(prn, "%s: ", I_("LR over-identification test"));
		pprintf(prn, "%s(%d) = %g [%.4f]\n\n", I_("Chi-square"),
			idf, X2, pv);
	    } else {
		pprintf(prn, "%s: ", _("LR over-identification test"));
		pprintf(prn, "%s(%d) = %g [%.4f]\n\n", _("Chi-square"),
			idf, X2, pv);
	    }
	} else if (idf == 0) {
	    pprintf(prn, "%s\n\n", 
		    (plain_format(prn))? _("Equation is just identified") :
		    I_("Equation is just identified"));
	}
    }
}

static void print_panel_AR_test (double z, int order, PRN *prn)
{
    double pv = normal_pvalue_2(z);

    if (na(pv)) {
	return;
    }

    if (plain_format(prn)) {
	pprintf(prn, _("Test for AR(%d) errors:"), order);
    } else {
	pprintf(prn, I_("Test for AR(%d) errors:"), order);
    }

    if (tex_format(prn)) {
	char numstr[32];

	tex_sprint_double_digits(z, numstr, 4);
	pprintf(prn, " & $z$ = %s [%.4f]", numstr, pv);
    } else {
	pprintf(prn, " z = %g [%.4f]", z, pv);
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

    pv = chisq_cdf_comp(df, x);

    if (na(pv)) {
	return;
    }

    if (tex_format(prn)) {
	pprintf(prn, "%s: ", I_(texstrs[j]));
	pprintf(prn, " & $\\chi^2(%d)$ = %g [%.4f]", df, x, pv);
    } else if (plain_format(prn)) {
	if (pmod->ci == GMM) {
	    pputs(prn, "  ");
	}
	pprintf(prn, "%s: ", _(strs[j]));
	pprintf(prn, "%s(%d) = %g [%.4f]", _("Chi-square"), df, x, pv);
    } else {
	pprintf(prn, "%s: ", I_(strs[j]));
	pprintf(prn, "%s(%d) = %g [%.4f]", I_("Chi-square"), df, x, pv);
    }

    gretl_prn_newline(prn);
}

static int GMM_crit_line (const MODEL *pmod, PRN *prn)
{
    double Q = pmod->ess;
    double TQ = pmod->ess * pmod->nobs;

    if (plain_format(prn)) {    
	pprintf(prn, "  %s: Q = %.*g (TQ = %.*g)\n", _("GMM criterion"), 
		XDIGITS(pmod), Q, XDIGITS(pmod), TQ);
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s: Q = %g (TQ = %g)\n", I_("GMM criterion"), 
		Q, TQ);
    } else if (tex_format(prn)) {
	char x1[32], x2[32];

	tex_sprint_double(Q, x1);
	tex_sprint_double(TQ, x2);
	pprintf(prn, "%s, $Q$ = %s ($TQ$ = %s)\\\\\n",
		I_("GMM criterion"), x1, x2);
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", I_("GMM criterion"), 
		prn_delim(prn), Q);
    }	

    return 0;
}

static void print_GMM_stats (const MODEL *pmod, PRN *prn)
{
    double x;

    ensure_vsep(prn);

    GMM_crit_line(pmod, prn);
    x = gretl_model_get_double(pmod, "J_test");
    if (!na(x)) {
	print_GMM_chi2_test(pmod, x, J_TEST, prn);
    }

    gretl_prn_newline(prn);
}

static void print_DPD_stats (const MODEL *pmod, PRN *prn)
{
    double x;

    ensure_vsep(prn);

    if (tex_format(prn)) {
	pputs(prn, "\\begin{tabular}{ll}\n");
    }

    x = gretl_model_get_double(pmod, "AR1");
    if (!na(x)) {
	print_panel_AR_test(x, 1, prn);
    }

    x = gretl_model_get_double(pmod, "AR2");
    if (!na(x)) {
	print_panel_AR_test(x, 2, prn);
    }

    x = gretl_model_get_double(pmod, "sargan");
    if (!na(x)) {
	print_GMM_chi2_test(pmod, x, AB_SARGAN, prn);
    }

    x = gretl_model_get_double(pmod, "wald");
    if (!na(x)) {
	print_GMM_chi2_test(pmod, x, AB_WALD, prn);
    }

    if (tex_format(prn)) {
	pputs(prn, "\\end{tabular}\n");
    } else {
	gretl_prn_newline(prn);
    }
}

static void maybe_print_lad_warning (const MODEL *pmod, PRN *prn)
{
    if (gretl_model_get_int(pmod, "nonunique")) {
	pputs(prn, _("Warning: solution is probably not unique"));
	pputc(prn, '\n');
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

    ensure_vsep(prn);

    if (plain_format(prn)) {
	pprintf(prn, "%s = %g\n", _("'Within' variance"), ws2);
	pprintf(prn, "%s = %g\n", _("'Between' variance"), bs2);
	if (!na(theta)) {
	    pprintf(prn, "%s = %g\n", _("theta used for quasi-demeaning"), theta);
	}
	pputc(prn, '\n');
    } else if (tex_format(prn)) {
	char xstr[32];

	tex_sprint_double(ws2, xstr);
	pprintf(prn, "$\\hat{\\sigma}^2_{\\varepsilon}$ = %s \\\\\n", xstr);
	tex_sprint_double(bs2, xstr);
	pprintf(prn, "$\\hat{\\sigma}^2_u$ = %s \\\\\n", xstr);
	if (!na(theta)) {
	    tex_sprint_double(theta, xstr);
	    pprintf(prn, "$\\theta$ = %s \\\\\n", xstr);
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

static double durbins_h (const MODEL *pmod)
{
    int ldv = gretl_model_get_int(pmod, "ldepvar");
    double se = pmod->sderr[ldv - 2];
    int T = pmod->nobs - 1;
    double h = NADBL;

    if (pmod->ess <= 0.0 || na(se) || (T * se * se) >= 1.0 ||
	na(pmod->rho)) {
	; /* can't calculate h */
    } else {
	h = pmod->rho * sqrt(T / (1 - T * se * se));
	if (xna(h)) {
	    h = NADBL;
	}
    }

    return h;
}

static int least_significant_coeff (const MODEL *pmod)
{
    double x, tmin = 3.2;
    int i, k = 0;

    if (gretl_list_separator_position(pmod->list) > 0) {
	return 0;
    }
    
    for (i=pmod->ifc; i<pmod->ncoeff; i++) {
	if (pmod->sderr[i] > 0) {
	    x = fabs(pmod->coeff[i] / pmod->sderr[i]);
	    if (x < tmin) {
		tmin = x;
		k = i;
	    }
	}
    }

    if (tmin < 3.2) {
	x = coeff_pval(pmod->ci, tmin, pmod->dfd);
	if (!na(x) && x > .10) {
	    return pmod->list[k+2];
	}
    }

    return 0;
}

static void pval_max_line (const MODEL *pmod, const DATAINFO *pdinfo, 
			   PRN *prn)
{
    int k = pmod->ncoeff - pmod->ifc;

    if (k < 3) return;

    if ((k = least_significant_coeff(pmod))) {
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
    } else if (aux == AUX_BP) {
	return N_("Breusch-Pagan test for heteroskedasticity");
    } else if (aux == AUX_HET_1) {
	return N_("Pesaran-Taylor test for heteroskedasticity");
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
#if 1
    else if (ci == HSK)  return N_("WLS"); 
#else
    else if (ci == HSK)  return N_("Heteroskedasticity-corrected");
#endif
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
    else if (ci == INTREG) return N_("Interval");
    else if (ci == ARBOND) {
	if (tex_format(prn)) return N_("Arellano--Bond");
	else return N_("Arellano-Bond");
    } else {
	return "";
    }
}

const char *estimator_string (const MODEL *pmod, PRN *prn)
{
    if (pmod->ci == AR1) {
	if (gretl_model_get_int(pmod, "hilu")) {
	    if (tex_format(prn)) return N_("Hildreth--Lu");
	    else return N_("Hildreth-Lu");
	} else if (gretl_model_get_int(pmod, "pwe")) {
	    if (tex_format(prn)) return N_("Prais--Winsten");
	    else return N_("Prais-Winsten");
	} else {
	    if (tex_format(prn)) return N_("Cochrane--Orcutt");
	    else return N_("Cochrane-Orcutt");
	}
    } else if (pmod->ci == ARMA) {
	if (gretl_model_get_int(pmod, "armax")) {
	    return N_("ARMAX");
	} else if (gretl_model_get_int(pmod, "arima_d") ||
		   gretl_model_get_int(pmod, "arima_D")) {
	    return N_("ARIMA");
	} else {
	    return N_("ARMA");
	}
    } else if (POOLED_MODEL(pmod)) {
	return N_("Pooled OLS");
    } else if (pmod->ci == PANEL) {
	if (gretl_model_get_int(pmod, "fixed-effects")) {
	    return N_("Fixed-effects");
	} else if (gretl_model_get_int(pmod, "random-effects")) {
	    return N_("Random-effects (GLS)");
	} else if (gretl_model_get_int(pmod, "unit-weights")) {
	    if (gretl_model_get_int(pmod, "iters")) {
		return N_("Maximum Likelihood");
	    } else {
		return N_("WLS");
	    }
	} else {
	    return N_("Between-groups");
	}
    } else if (pmod->ci == ARBOND) {
	if (gretl_model_get_int(pmod, "step") == 2) {
	    return N_("2-step Arellano-Bond");
	} else {
	    return N_("1-step Arellano-Bond");
	}
    } else if (pmod->ci == GMM) {
	if (gretl_model_get_int(pmod, "two-step")) {
	    return N_("2-step GMM");
	} else if (gretl_model_get_int(pmod, "iterated")) {
	    return N_("Iterated GMM");
	} else if (gretl_model_get_int(pmod, "step") == 2) {
	    return N_("2-step GMM");
	} else if (gretl_model_get_int(pmod, "step") > 2) {
	    return N_("Iterated GMM");
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
    } else if (pmod->ci == HECKIT) {
	if (gretl_model_get_int(pmod, "two-step")) {
	    return N_("Two-step Heckit");
	} else {
	    return N_("ML Heckit");
	}
    } else if (pmod->ci == LAD) {
	if (gretl_model_get_int(pmod, "rq")) {
	    return N_("Quantile");
	} else {
	    return N_("LAD");
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

	tex_sprint_double(F, x1str);
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
    int i, gotsep = 0;
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

    ccount += char_len(_("Instruments") + 2);

    for (i=2; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    gotsep = 1;
	    continue;
	}
	if (gotsep) {
	    if (tex) {
		tex_escape(vname, pdinfo->varname[list[i]]);
	    } else {
		strcpy(vname, pdinfo->varname[list[i]]);
	    }
	    pprintf(prn, "%s ", vname);
	    ccount += char_len(vname) + 1;
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

    if (ninst == 0) {
	pputs(prn, _("none"));
    }

    if (ccount > 0) {
	gretl_prn_newline(prn);
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
    int k = gretl_model_get_int(pmod, "hac_kernel");
    int h = gretl_model_get_int(pmod, "hac_lag");
    int white = gretl_model_get_int(pmod, "hac_prewhiten");
    int utf = plain_format(prn);
    double bt;

    if (k == KERNEL_QS) {
	bt = gretl_model_get_double(pmod, "qs_bandwidth");
	if (utf) {
	    pprintf(prn, _("HAC standard errors, "
			   "bandwidth %.2f"), bt);
	} else {
	    pprintf(prn, I_("HAC standard errors, "
			    "bandwidth %.2f"), bt);
	}
    } else {
	if (utf) {
	    pprintf(prn, _("HAC standard errors, "
			   "bandwidth %d"), h);
	} else {
	    pprintf(prn, I_("HAC standard errors, "
			    "bandwidth %d"), h);
	}
    }

    pputc(prn, ' ');
    if (utf) {
	pprintf(prn, "(%s", _(kstrs[k]));
    } else {
	pprintf(prn, "(%s", I_(kstrs[k]));
		    
    }
    if (white) {
	pputs(prn, ", ");
	pputs(prn, (utf)? _("prewhitened") : 
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
    const char *s = NULL;

    switch (v) {
    case VCV_HESSIAN:
	s = N_("Standard errors based on Hessian");
	break;
    case VCV_IM:
	s = N_("Standard errors based on Information Matrix");
	break;
    case VCV_OP:
	s = N_("Standard errors based on Outer Products matrix");
	break;
    case VCV_QML:
	s = N_("QML standard errors");
	break;
    case VCV_BW:
	if (tex) {
	    s = N_("Bollerslev--Wooldridge standard errors");
	} else {
	    s = N_("Bollerslev-Wooldridge standard errors");
	}
	break;
    default:
	break;
    }

    if (s != NULL) {
	if (csv_format(prn)) {
	    pprintf(prn, "\"%s\"\n", I_(s));
	} else {
	    pprintf(prn, "%s\n", (utf)? _(s) : I_(s));
	}
    }
}

static void rq_vcv_line (const MODEL *pmod, PRN *prn)
{
    int robust = gretl_model_get_int(pmod, "rq_nid");
    double a = gretl_model_get_double(pmod, "rq_alpha");
    int utf = plain_format(prn);
    int free_s = 0;
    char *s;

    if (!na(a)) {
	if (robust) {
	    s = g_strdup_printf(N_("With robust %g percent confidence intervals"), 
				100 * (1 - a));
	} else {
	    s = g_strdup_printf(N_("With %g percent confidence intervals"), 
				100 * (1 - a));
	}
	free_s = 1;
    } else if (robust) {
	s = N_("Robust (sandwich) standard errors");
    } else {
	s = N_("Asymptotic standard errors assuming IID errors");
    }

    if (csv_format(prn)) {
	pprintf(prn, "\"%s\"", I_(s));
    } else {
	pprintf(prn, "%s", (utf)? _(s) : I_(s));
    }

    gretl_prn_newline(prn);

    if (free_s) {
	g_free(s);
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
    if (pmod->ci == LAD && gretl_model_get_int(pmod, "rq")) {
	rq_vcv_line(pmod, prn);
    } else if (gretl_model_get_int(pmod, "using_hac") ||
	gretl_model_get_int(pmod, "hac_kernel") ||
	gretl_model_get_int(pmod, "hac_lag")) {
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
    } else if (acode & ARMA_LS) {
	pputs(prn, _("Estimated using least squares"));
	pputs(prn, " (");
	pputs(prn, _("conditional ML"));
	pputs(prn, ")\n");
    } else {
	pputs(prn, _("Estimated using BHHH method"));
	pputs(prn, " (");
	pputs(prn, _("conditional ML"));
	pputs(prn, ")\n");
    }	
}

static void godfrey_test_string (int ci, int order, int utf, PRN *prn)
{
    pputc(prn, '\n');

    if (ci == TSLS) {
	if (utf) { 
	    if (order > 1) {
		pprintf(prn, _("Godfrey (1994) test for autocorrelation up to order %d"), 
			order);
	    } else {
		pputs(prn, _("Godfrey (1994) test for first-order autocorrelation"));
	    }
	} else {
	    if (order > 1) {
		pprintf(prn, I_("Godfrey (1994) test for autocorrelation up to order %d"), 
			order);
	    } else {
		pputs(prn, I_("Godfrey (1994) test for first-order autocorrelation"));
	    }
	} 
    } else {
	if (utf) { 
	    if (order > 1) {
		pprintf(prn, _("Breusch-Godfrey test for autocorrelation up to order %d"), 
			order);
	    } else {
		pputs(prn, _("Breusch-Godfrey test for first-order autocorrelation"));
	    }
	} else {
	    if (order > 1) {
		pprintf(prn, I_("Breusch-Godfrey test for autocorrelation up to order %d"), 
			order);
	    } else {
		pputs(prn, I_("Breusch-Godfrey test for first-order autocorrelation"));
	    }
	} 
    }

    pputc(prn, '\n');
}

static void print_intreg_depvar (const MODEL *pmod,
				 const DATAINFO *pdinfo,
				 PRN *prn)
{
    char *lov = (char *) gretl_model_get_data(pmod, "lovar");
    char *hiv = (char *) gretl_model_get_data(pmod, "hivar");

    pprintf(prn, "%s: %s", I_("Lower limit"), lov);
    pprintf(prn, ", %s: %s", I_("Upper limit"), hiv);
}

static void print_model_heading (const MODEL *pmod, 
				 const DATAINFO *pdinfo, 
				 gretlopt opt, 
				 PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN], vname[32];
    int t1 = pmod->t1, t2 = pmod->t2;
    int tex = tex_format(prn);
    int plain = plain_format(prn);
    int utf = plain;
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
    case AUX_BP:
    case AUX_HET_1:	
    case AUX_CHOW:
    case AUX_COINT:
    case AUX_ADF:
    case AUX_DF:
    case AUX_KPSS:
    case AUX_RESET:
    case AUX_GROUPWISE:
	if (plain) {
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
	godfrey_test_string(pmod->ci, order, utf, prn);
	break;	
    case AUX_ARCH:
	order = gretl_model_get_int(pmod, "arch_order");
	pputc(prn, '\n');
	pprintf(prn, (utf)? _("Test for ARCH of order %d") :
		I_("Test for ARCH of order %d"), order);
	pputc(prn, '\n');
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
    case AUX_AUX:
	pputc(prn, '\n');
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
		_(system_short_string(pmod)),
		pmod->nobs, startdate, (tex)? "--" : "-", enddate);
    } else if (!dataset_is_panel(pdinfo)) {
	const char *estr = estimator_string(pmod, prn);
	const char *fmt;

	if (pmod->missmask != NULL) {
	    int mc = model_missval_count(pmod);

	    if (pmod->ci == HECKIT) {
		int Tmax = pmod->t2 - pmod->t1 + 1;
		
		mc = Tmax - gretl_model_get_int(pmod, "totobs");
	    }

	    if (char_len(estr) > 24) {
		fmt = N_("%s estimates %s%s%s (T = %d)");
		pprintf(prn, (utf)? _(fmt) : I_(fmt), (utf)? _(estr) : I_(estr),
			startdate, (tex)? "--" : "-", enddate, pmod->nobs);
	    } else {
		fmt = N_("%s estimates using %d observations from %s%s%s");
		pprintf(prn, (utf)? _(fmt) : I_(fmt), (utf)? _(estr) : I_(estr),
			pmod->nobs, startdate, (tex)? "--" : "-", enddate);
	    }
	    if (mc > 0) {
		gretl_prn_newline(prn);
		pprintf(prn, "%s: %d",
			(utf)? _("Missing or incomplete observations dropped") :
			I_("Missing or incomplete observations dropped"), mc);
	    }
	} else {
	    if (char_len(estr) > 24) {
		fmt = N_("%s estimates %s%s%s (T = %d)");
		pprintf(prn, (utf)? _(fmt) : I_(fmt), (utf)? _(estr) : I_(estr), 
			startdate, (tex)? "--" : "-", enddate, pmod->nobs);

	    } else {
		fmt = N_("%s estimates using the %d observations %s%s%s");
		pprintf(prn, (utf)? _(fmt) : I_(fmt), (utf)? _(estr) : I_(estr), 
			pmod->nobs, startdate, (tex)? "--" : "-", enddate);
	    }
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
    } else if (pmod->aux == AUX_WHITE || pmod->aux == AUX_HET_1) {
	pprintf(prn, "%s: %s", 
		(utf)? _("Dependent variable") : I_("Dependent variable"),
		(tex)? "$\\hat{u}^2$" : "uhat^2");
    } else if (pmod->aux == AUX_BP) {
	const char *fmt;

	if (gretl_model_get_int(pmod, "robust")) {
	    fmt = N_("scaled %s (Koenker robust variant)");
	} else {
	    fmt = N_("scaled %s");
	}
	pprintf(prn, "%s: ", (utf)? _("Dependent variable") : I_("Dependent variable"));
	pprintf(prn, (utf)? _(fmt) : I_(fmt), (tex)? "$\\hat{u}^2$" : "uhat^2");
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
    } else if (pmod->ci == INTREG) {
	print_intreg_depvar(pmod, pdinfo, prn);
    } else { 
	const char *dvname = 
	    gretl_model_get_depvar_name(pmod, pdinfo);

	if (tex) {
	    if (pmod->aux == AUX_VECM) {
		tex_vecm_depvar_name(vname, dvname);
	    } else if (pmod->ci == ARBOND) {
		tex_arbond_depvar_name(vname, dvname);
	    } else {
		tex_escape(vname, dvname);
	    }
	}

	if (pmod->aux == AUX_VAR || pmod->aux == AUX_VECM) {
	    pputs(prn, (tex)? vname : dvname);
	} else {
	    pprintf(prn, "%s: %s", 
		    (utf)? _("Dependent variable") : I_("Dependent variable"),
		    (tex)? vname : dvname);
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

	if (method != SYS_METHOD_FIML && method != SYS_METHOD_LIML) {
	    print_tsls_instruments(pmod->list, pdinfo, prn);
	}
    }

    /* tau for quantile regression */
    else if (pmod->ci == LAD) {
	double tau = gretl_model_get_double(pmod, "tau");

	if (!na(tau)) {
	    if (tex) {
		pprintf(prn, "$\\tau$ = %g", tau);
	    } else {
		pprintf(prn, "tau = %g", tau);
	    }
	    gretl_prn_newline(prn);
	}
    }

    /* VCV variants */
    print_model_vcv_info(pmod, prn);

    /* WLS on panel data */
    if (gretl_model_get_int(pmod, "unit-weights") && !pmod->aux) {
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

    /* rhohat for AR1 (TeX) */
    else if (pmod->ci == AR1) {
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

    if (plain_format(prn) && pmod->ci == LAD) {
	maybe_print_lad_warning(pmod, prn);
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
	} 
	pputs(prn, "\\begin{center}\n");
    } else if (rtf_format(prn)) {
	if (rtf_doc_format(prn)) {
	    pputs(prn, "{\\rtf1\\par\n\\qc ");
	} else {
	    pputs(prn, "\\par\n\\qc ");
	}
    } else if (plain_format(prn)) {
#ifdef ENABLE_NLS
	set_gui_native_printing();
#endif
    }
}

#define RTF_COEFF_ROW  "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                       "\\cellx1900\\cellx3300\\cellx4700\\cellx6100" \
                       "\\cellx7500\\cellx8000\n\\intbl"

#define RTF_BINARY_ROW "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                       "\\cellx1900\\cellx3300\\cellx4700\\cellx6100" \
                       "\\cellx8000\n\\intbl"

#define RTF_INTVL_ROW "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                       "\\cellx1900\\cellx3300\\cellx4700\\cellx6100" \
                       "\n\\intbl"

#define RTF_ROOT_ROW   "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                       "\\cellx500\\cellx1500\\cellx2900\\cellx4300" \
                       "\\cellx5700\\cellx7100\n\\intbl"

/* below: this is used when we're doing something other than a plain
   text print of a model */

static void alt_print_coeff_table_start (const MODEL *pmod, int ci, PRN *prn)
{
    const char *tlabel;
    int use_param = 0;
    int slopes = 0;
    int intervals = 0;
    int seqcols = 0;
    int mp = 0;

    if (ci == MODPRINT || ASYMPTOTIC_MODEL(ci)) {
	tlabel = (tex_format(prn))? N_("$z$-stat") : N_("z-stat");
    } else {
	tlabel = (tex_format(prn))? N_("$t$-ratio") : N_("t-ratio");
    }

    if (pmod != NULL) {
	gretl_matrix *m;

	use_param = nonlin_model(pmod);
	slopes = binary_model(pmod) && !gretl_model_get_int(pmod, "show-pvals");
	intervals = gretl_model_get_data(pmod, "coeff_intervals") != NULL;
	m = gretl_model_get_data(pmod, "rq_sequence");
	seqcols = gretl_matrix_cols(m);
	mp = (pmod->ci == MPOLS);
    }

    if (csv_format(prn)) {
	char d = prn_delim(prn);

	if (mp) {
	    pprintf(prn, "%c\"%s\"%c\"%s\"\n",
		    d, I_("coefficient"), d, I_("std. error"));
	} else if (slopes) {
	    pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    d, I_("coefficient"), d, I_("std. error"),
		    d, I_(tlabel), d, I_("slope at mean"));
	} else if (use_param) {
	    pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    d, I_("estimate"), d, I_("std. error"),
		    d, I_(tlabel), d, I_("p-value"));
	} else if (intervals) {
	    pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    d, I_("coefficient"), d, I_("lower"),
		    d, I_("upper"));
	} else if (seqcols == 3) {
	    pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    d, "tau", d, I_("coefficient"), d, I_("lower"),
		    d, I_("upper"));
	} else if (seqcols == 2) {
	    pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    d, I_("coefficient"), d, I_("std. error"),
		    d, I_(tlabel));
	} else {
	    pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    d, I_("coefficient"), d, I_("std. error"),
		    d, I_(tlabel), d, I_("p-value"));
	}	
    } else {
	const char *cols[6] = { NULL };
	int i;

	cols[0] = " ";
	cols[1] = (use_param)? N_("Estimate") : N_("Coefficient");

	if (intervals) {
	    cols[2] = N_("Lower");
	    cols[3] = N_("Upper");
	} else if (seqcols == 3) {
	    cols[2] = N_("tau");
	    cols[3] = N_("Lower");
	    cols[4] = N_("Upper");
	} else {
	    cols[2] = (tex_format(prn))? N_("Std.\\ Error") : N_("Std. Error");
	    if (!mp) {
		cols[3] = tlabel;
		cols[4] = (slopes)? N_("Slope") : N_("p-value");
	    }
	}

	if (tex_format(prn)) {
	    gretlopt tabopt = (ci == MODPRINT)? OPT_U : OPT_NONE;

	    if (slopes) {
		tabopt |= OPT_B; /* "binary" */
	    }
	    if (mp) {
		tabopt |= OPT_M; /* multiple precision */
	    }	    
	    tex_coeff_table_start(cols, tabopt, prn);
	    return;
	}   

	if (rtf_format(prn)) {
	    pputc(prn, '{');
	    if (intervals) {
		pputs(prn, RTF_INTVL_ROW);
	    } else if (slopes) {
		pputs(prn, RTF_BINARY_ROW);
	    } else {
		pputs(prn, RTF_COEFF_ROW);
	    }
	    for (i=0; cols[i] != NULL; i++) {
		if (slopes && i == 4) {
		    pprintf(prn, " \\qc {\\i %s{\\super *}}\\cell", I_(cols[i]));
		} else {
		    pprintf(prn, " \\qc {\\i %s}\\cell", I_(cols[i]));
		}
	    }
	    if (!slopes && !intervals) {
		pputs(prn, " \\ql \\cell");
	    }
	    pputs(prn, " \\intbl \\row\n");
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
    } else if (plain_format(prn)) {
#ifdef ENABLE_NLS
	unset_gui_native_printing();
#endif
    }
} 

static void addconst_message (const MODEL *pmod, PRN *prn)
{
    if (gretl_model_get_int(pmod, "addconst")) {
	pprintf(prn, "\n\n%s\n", 
		_("WARNING: The constant was present among the regressors but not among the\n"
		  "instruments, so it has been automatically added to the instrument list.\n"
		  "This behavior may change in future versions, so you may want to adjust your\n"
		  "scripts accordingly.\n"));
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

enum {
    ORIG_STATS,
    WTD_STATS,
    RHODIFF_STATS,
    TRANSFORM_STATS
};

static void alternate_stats_message (int i, PRN *prn)
{
    const char *msg[] = {
	N_("Statistics based on the original data"),
	N_("Statistics based on the weighted data"),
	N_("Statistics based on the rho-differenced data"),
	N_("Statistics based on the transformed data")
    };

    if (plain_format(prn)) {
	pprintf(prn, "%s:\n\n", _(msg[i]));
    } else if (tex_format(prn)) {
	pprintf(prn, "\\vspace{1em}%s:\n\n", _(msg[i]));
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"\n", _(msg[i]));
    } else { 
	/* RTF */
	pprintf(prn, "\\par \\qc\n%s:\n\n", I_(msg[i]));	
    }
}

static void print_whites_results (const MODEL *pmod, PRN *prn)
{
    double X = pmod->rsq * pmod->nobs;
    int df = pmod->ncoeff - 1;
    double pv = chisq_cdf_comp(df, X);

    if (plain_format(prn)) {
	pprintf(prn, "\n%s: TR^2 = %f,\n", _("Test statistic"), X);
	pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n", 
		_("with p-value"), _("Chi-square"), df, X, pv);
    } else if (rtf_format(prn)) { /* FIXME */
	pprintf(prn, "\\par \\ql\n%s: TR{\\super 2} = %f,\n", I_("Test statistic"), 
		X);
	pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n", 
		I_("with p-value"), I_("Chi-square"), df, X, pv);
    } else if (tex_format(prn)) {
	pprintf(prn, "\n%s: $TR^2$ = %f,\n", I_("Test statistic"), X);
	pprintf(prn, "%s = $P$($\\chi^2(%d)$ > %f) = %f\n\n",
		I_("with p-value"), df, X, pv);
    }
}

static void print_bp_results (const MODEL *pmod, PRN *prn)
{
    double pv, X = gretl_model_get_double(pmod, "BPLM");
    int df = pmod->ncoeff - 1;

    if (na(X)) {
	return;
    }

    pv = chisq_cdf_comp(df, X);

    if (plain_format(prn)) {
	pprintf(prn, "\n%s: LM = %f,\n", _("Test statistic"), X);
	pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n", 
		_("with p-value"), _("Chi-square"), df, X, pv);
    } else if (rtf_format(prn)) { /* FIXME */
	pprintf(prn, "\\par \\ql\n%s: LM = %f,\n", I_("Test statistic"), 
		X);
	pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n", 
		I_("with p-value"), I_("Chi-square"), df, X, pv);
    } else if (tex_format(prn)) {
	pprintf(prn, "\n%s: LM = %f,\n", I_("Test statistic"), X);
	pprintf(prn, "%s = $P$($\\chi^2(%d)$ > %f) = %f\n\n",
		I_("with p-value"), df, X, pv);
    }
}

static void print_HET_1_results (const MODEL *pmod, PRN *prn)
{
    double z = fabs(pmod->coeff[1]) / pmod->sderr[1];
    double pv = 2.0 * (1 - normal_cdf(z));

    if (plain_format(prn)) {
	pprintf(prn, "\n%s: HET_1 = |%f| / %f = %f,\n", _("Test statistic"), 
		pmod->coeff[1], pmod->sderr[1], z);
	pprintf(prn, "%s = 2 * P(z > %f) = %.3g\n\n", 
		_("with p-value"), z, pv);
    } else if (rtf_format(prn)) { /* FIXME */
	pprintf(prn, "\\par \\ql\n%s: HET_1 = %f,\n", I_("Test statistic"), z);
	pprintf(prn, "%s = 2 * P(z > %f) = %.3g\n\n", 
		I_("with p-value"), z, pv);
    } else if (tex_format(prn)) {
	pprintf(prn, "\n%s: \verb|HET_1| = %f,\n", I_("Test statistic"), z);
	pprintf(prn, "%s = $2 \times P$($z$ > %f) = %f\n\n",
		I_("with p-value"), z, pv);
    }
}

static void maybe_print_jll (const MODEL *pmod, int lldig, PRN *prn)
{
    double jll = gretl_model_get_double(pmod, "jll");

    if (!na(jll)) {
	char jllstr[64];
	char xstr[32];

	sprintf(jllstr, _("Log-likelihood for %s"), 
		(const char *) gretl_model_get_data(pmod, "log-parent"));
	if (lldig > 0) {
	    plain_print_double(xstr, lldig, jll, prn);
	    pprintf(prn, "  (%s = %s)\n", jllstr, xstr);
	} else {
	    plain_print_double(xstr, XDIGITS(pmod), jll, prn);
	    pprintf(prn, "%s = %s\n\n", jllstr, xstr);
	}	    
    }
}

#define fixed_effects_model(m) (m->ci == PANEL && \
                                gretl_model_get_int(m, "fixed-effects"))

#define random_effects_model(m) (m->ci == PANEL && \
                                 gretl_model_get_int(m, "random-effects"))

#define between_model(m) (m->ci == PANEL && \
                          gretl_model_get_int(m, "between"))

#define weighted_model(m) (m->ci == HSK || m->ci == ARCH || \
			   (m->ci == WLS && !gretl_model_get_int(m, "wt_dummy")) || \
                           (m->ci == PANEL && gretl_model_get_int(m, "unit-weights")))

#define panel_ML_model(m) (m->ci == PANEL && \
                           gretl_model_get_int(m, "unit-weights") && \
			   gretl_model_get_int(m, "iters"))

#define non_weighted_panel(m) (m->ci == PANEL && \
			       !gretl_model_get_int(m, "unit-weights"))
#define MID_STATS 14

struct middletab {
    const char **key;     /* stats strings */
    double *val;          /* stats values */
    int minus;            /* variant of minus sign in use */
    int mlen;             /* max length of translated string */
    int ipos[MID_STATS];  /* is i-th stat integer? */
    int nls;              /* translation on? (0/1) */
    int d;                /* CSV field delimiter */
    int multi;            /* Using multiple precision? (0/1) */
    char txt_fmt[36];     /* format for plain text output */
};

#ifdef ENABLE_NLS

static void set_mtab_string_width (struct middletab *mt)
{
    int len, maxlen = 0;
    int i, j;

    for (i=0, j=0; i<7; i++, j+=2) {
	if (!na(mt->val[j])) {
	    len = g_utf8_strlen(_(mt->key[j]), -1);
	    if (len > maxlen) maxlen = len;
	    len = g_utf8_strlen(_(mt->key[j+1]), -1);
	    if (len > maxlen) maxlen = len;
	}
    }

    if (maxlen > 22) {
	fprintf(stderr, "Can't make compact model stats table -- the max\n"
		"length translated string is %d chars, should be < 23\n",
		maxlen);
    }

    mt->mlen = maxlen;
}

#endif

static void middletab_prepare_format (struct middletab *mt, int j)
{
    const char *s1 = mt->key[j];
    const char *s2 = mt->key[j+1];
    int d1 = 0, d2 = 0;
    int len = mt->mlen;

    if (mt->nls) {
	/* calculate bytes minus glyphs */
	d1 = strlen(_(s1)) - g_utf8_strlen(_(s1), -1);
	d2 = strlen(_(s2)) - g_utf8_strlen(_(s2), -1);
    } 

    if (mt->multi) {
	len = (len < 20)? 20 : len + 1;
	sprintf(mt->txt_fmt, "  %%-%ds  %%s\n  %%-%ds  %%s\n", len + d1, len + d2);
    } else if (len < 20) {
	sprintf(mt->txt_fmt, "%%-%ds%%s   %%-%ds%%s\n", 20 + d1, 20 + d2);
    } else if (len < 23) {
	sprintf(mt->txt_fmt, "%%-%ds%%s   %%-%ds%%s\n", len + d1 + 1, len + d2 + 1);
    } else {
	sprintf(mt->txt_fmt, "  %%-%ds  %%s\n  %%-%ds  %%s\n", len + d1, len + d2);
    }
}

static char *print_csv (char *s, double x)
{
    if (na(x)) {
	strcpy(s, "NA");
    } else {
	sprintf(s, "%.15g", x);
    }

    return s;
} 

enum {
    MINUS_HYPHEN,
    MINUS_UTF,
    MINUS_TEX
};

#define RSQ_POS 4

static void mtab_numstart (char *s, double x, int minus)
{
    if (x < 0) {
	if (minus == MINUS_UTF) {
	    strcpy(s, "−"); /* U+2212: minus */
	} else if (minus == MINUS_TEX) {
	    strcpy(s, "$-$");
	} else {
	    strcpy(s, "-"); /* ASCII: use hyphen */
	}
    } else {
	strcpy(s, " "); 
    }
}

/* Try to pack as much as we can of a given number into a fixed width
   of 8 characters (a leading minus, if needed, is a ninth).
*/

static char *print_eight (char *s, struct middletab *mt, int i)
{
    double ax, x = mt->val[i];
    char tmp[16];

    if (mt->ipos[i]) {
	sprintf(s, "%9d", (int) x);
	return s;
    }    

    if (i == RSQ_POS) {
	/* R-squared: don't use scientific notation */
	sprintf(s, "%9.6f", x);
	return s;
    }

    if (na(x)) {
	sprintf(s, "%9s", "NA");
	return s;
    }

    mtab_numstart(s, x, mt->minus);

    ax = fabs(x);

    if (ax < 0.00001 || ax > 99999999) {
	char *p;

	sprintf(tmp, "%#.3g", ax);
	p = strrchr(tmp, (ax < 1)? '-' : '+');
	if (p == NULL) {
	    sprintf(tmp, "%8.6f", ax);
	} else if (strlen(p) == 4) {
	    if (*(p+1) == '0') {
		memcpy(p+1, p+2, 3);
	    } else {
		sprintf(tmp, "%#.2g", ax);
	    }
	}
    } else if (ax < 10.0) {
	sprintf(tmp, "%8.6f", ax);
    } else {
	double lx = log10(ax);
	int ldig = ceil(lx) + (lx == floor(lx));

	ldig = (ldig > 7)? 7 : ldig;
	sprintf(tmp, "%8.*f", 8 - ldig - 1, ax);
    }

    strcat(s, tmp);

    if (mt->minus == MINUS_TEX && strchr(s, 'e') != NULL) {
	tex_modify_exponent(s);
    }

    return s;
} 

static char *print_fifteen (char *s, double x, int minus)
{
    if (na(x)) {
	strcpy(s, " NA");
    } else if (minus == MINUS_TEX) {
	char *p;

	if (x < 0) {
	    sprintf(s, "$-$%.15E", -x);
	} else {
	    sprintf(s, "%.15E", x);
	}

	if ((p = strstr(s, "E-")) != NULL) {
	    char tmp[8];

	    sprintf(tmp, "E--%s", p + 2);
	    strcpy(p, tmp);
	}
    } else if (minus == MINUS_UTF) {
	if (x < 0) {
	    sprintf(s, "−%.15E", -x);
	} else {
	    sprintf(s, "% .15E", x);
	}	
    } else {
	sprintf(s, "% .15E", x);
    }

    return s;
} 

#define RTF_MT_ROW "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                   "\\cellx2500\\cellx3800\\cellx4200\\cellx6700" \
                   "\\cellx8000\n"

#define RTF_MULTI_ROW "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                      "\\cellx2500\\cellx6000\n"

#define RTF_MT_FMT "\\intbl\\ql %s\\cell\\qr %s\\cell \\qc \\cell" \
                   "\\ql %s\\cell\\qr %s\\cell\\intbl\\row\n"

#define RTF_MULTI_FMT "\\intbl\\ql %s\\cell\\qr %s\\cell\\intbl\\row\n"

#define TXT_MT_FMT "%-20s%s   %-20s%s\n"

enum {
    MIDDLE_REGULAR,
    MIDDLE_TRANSFORM,
    MIDDLE_ORIG
};

static void middle_table_row (struct middletab *mt, int j, PRN *prn)
{
    const char *s1 = mt->key[j];
    const char *s2 = mt->key[j+1];
    char x1[48], x2[48];
    int k = j + 1;

    if (tex_format(prn)) {
	if (mt->multi) {
	    pprintf(prn, "%s & %s \\\\\n%s & %s \\\\\n",
		    _(s1), print_fifteen(x1, mt->val[j], mt->minus),
		    _(s2), print_fifteen(x2, mt->val[k], mt->minus));
	} else {
	    pprintf(prn, "%s & %s & %s & %s \\\\\n",
		    _(s1), print_eight(x1, mt, j),
		    _(s2), print_eight(x2, mt, k));
	}
    } else if (rtf_format(prn)) {
	if (mt->multi) {
	    pputs(prn, RTF_MULTI_ROW);
	    pprintf(prn, RTF_MULTI_FMT, I_(s1), 
		    print_fifteen(x1, mt->val[j], mt->minus));
	    pputs(prn, RTF_MULTI_ROW);
	    pprintf(prn, RTF_MULTI_FMT, I_(s2), 
		    print_fifteen(x1, mt->val[k], mt->minus));
	} else {
	    pputs(prn, RTF_MT_ROW);
	    pprintf(prn, RTF_MT_FMT, 
		    I_(s1), print_eight(x1, mt, j),
		    I_(s2), print_eight(x2, mt, k));
	}
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%s%c\"%s\"%c%s\n",
		_(s1), mt->d, print_csv(x1, mt->val[j]), mt->d,
		_(s2), mt->d, print_csv(x2, mt->val[k]));
    } else {
	if (mt->nls || mt->multi) {
	    middletab_prepare_format(mt, j);
	}
	if (mt->multi) {
	    pprintf(prn, mt->txt_fmt, 
		    _(s1), print_fifteen(x1, mt->val[j], mt->minus),
		    _(s2), print_fifteen(x2, mt->val[k], mt->minus));
	} else {
	    pprintf(prn, mt->txt_fmt, 
		    _(s1), print_eight(x1, mt, j),
		    _(s2), print_eight(x2, mt, k));
	}
    }
}

#ifdef ENABLE_NLS

/* If we haven't found a translation for a new-style short string, 
   use the corresponding old-style long string instead.
*/

static void maybe_remedy_translations (const char **S, int n)
{
    const char *old_key[] = {
	N_("Mean of dependent variable"),
	N_("Standard deviation of dep. var."),
	N_("Sum of squared residuals"),
	N_("Standard error of the regression"),
	N_("Unadjusted R-squared"),
	N_("Adjusted R-squared"),
	NULL,
	NULL,
	N_("Log-likelihood"),
	N_("Akaike information criterion"),
	N_("Schwarz Bayesian criterion"),
	NULL,
	NULL,
	NULL
    };
    int i;

    for (i=0; i<n; i++) {
	if (old_key[i] != NULL && !strcmp(S[i], _(S[i]))) {
	    fprintf(stderr, "untranslated: %s -> %s\n", S[i], old_key[i]);
	    S[i] = old_key[i];
	}
    }
}

#endif

/* print the block of statistics that appears beneath of the
   table of coefficients, standard errors, etc.
*/

static void print_middle_table (const MODEL *pmod, PRN *prn, int code)
{
    const char *note = N_("note on model statistics abbreviations here");
    int rtf = rtf_format(prn);
    int tex = tex_format(prn);
    int csv = csv_format(prn);
    char fstr[32], chistr[32];
    const char *key[] = {
	N_("Mean dependent var"),  /* 22: Mean of dependent variable */
	N_("S.D. dependent var"),  /* 22: Standard deviation of dependent var */
	N_("Sum squared resid"),   /* 22: Sum of squared residuals */
	N_("S.E. of regression"),  /* 22: Standard error of the regression */
	N_("R-squared"),           /* 22: */
	N_("Adjusted R-squared"),  /* 22: */
	"F-statistic",             /* will be replaced below */
	N_("P-value(F)"),          /* 22: P-value of F-statistic */	
	N_("Log-likelihood"),      /* 22: */
	N_("Akaike criterion"),    /* 22: Akaike Information Criterion */
	N_("Schwarz criterion"),   /* 22: Schwarz Bayesian Criterion */
	N_("Hannan-Quinn"),        /* 22: Hannan-Quinn Criterion */
	N_("rho"),                 /* 22: 1st-order autocorrelation coeff. */
	N_("Durbin-Watson")        /* 22: Durbin-Watson statistic */
    };
    double val[MID_STATS] = {
	pmod->ybar,
	pmod->sdy,
	pmod->ess,
	pmod->sigma,
	pmod->rsq,
	pmod->adjrsq,
	pmod->fstt,
	NADBL,
	pmod->lnL,
	pmod->criterion[C_AIC],
	pmod->criterion[C_BIC],
	pmod->criterion[C_HQC],
	pmod->rho,
	pmod->dw
    };
    struct middletab mtab;
    int i, j;

    mtab.mlen = 0;
    mtab.d = 0;
    mtab.minus = MINUS_HYPHEN;
    for (i=0; i<MID_STATS; i++) {
	mtab.ipos[i] = 0;
    }
    mtab.nls = 0;
    mtab.multi = (pmod->ci == MPOLS);

#ifdef ENABLE_NLS
    mtab.nls = doing_nls();
    if (mtab.nls) {
	maybe_remedy_translations(key, MID_STATS);
    }
#endif

    if (tex) {
	/* some special strings for TeX output */
	mtab.minus = MINUS_TEX;
	key[4] = "$R^2$";
	key[5] = N_("Adjusted $R^2$");
	key[7] = N_("P-value($F$)");
	key[11] = "Hannan--Quinn";
	key[12] = "$\\hat{\\rho}$";
	key[13] = "Durbin--Watson";
    } else if (!rtf && gretl_print_supports_utf(prn)) {
	/* print a 'real' minus sign? */
	mtab.minus = MINUS_UTF;
    }

    if (pmod->aux == AUX_VECM || pmod->aux == AUX_COINT) {
	/* VECM equation or Engle-Granger test: suppress F-test */
	val[6] = val[7] = NADBL;
    } else if (!na(pmod->fstt)) {
	/* format F-stat and get its p-value */
	if (tex) {
	    sprintf(fstr, "$F(%d, %d)$", pmod->dfn, pmod->dfd);
	} else {
	    sprintf(fstr, "F(%d, %d)", pmod->dfn, pmod->dfd);
	}
	key[6] = fstr;
	val[7] = snedecor_cdf_comp(pmod->dfn, pmod->dfd, pmod->fstt);
    }

    /* special variants of R-squared */
    if (gretl_model_get_int(pmod, "uncentered")) {
	key[4] = (tex)? N_("Uncentered $R^2$") : 
	    N_("Uncentered R-squared"); /* 22: */
	key[5] = (tex)? N_("Centered $R^2$") : 
	    N_("Centered R-squared");  /* 22: */
	val[5] = gretl_model_get_double(pmod, "centered-R2");
    } else if (pmod->ci == POISSON || binary_model(pmod)) {
	key[4] = (tex)? N_("McFadden $R^2$") : 
	    N_("McFadden R-squared");  /* 22: McFadden's pseudo-R-squared */
    }

    if (pmod->ci == ARBOND) {
	for (i=0; i<MID_STATS; i++) {
	    if (i < 2 || i > 3) {
		val[i] = NADBL;
	    }
	}	
    } else if (pmod->aux == AUX_SYS) {
	/* only dep. var. stats, SSR and SE, unless LIML */
	if (liml_equation(pmod) && !gretl_model_get_int(pmod, "restricted")) {
	    for (i=4; i<MID_STATS; i++) {
		if (i < 8 || i > 9) {
		    val[i] = NADBL;
		}
	    }
	    key[9] = N_("Smallest eigenvalue"); /* 22: */
	    val[9] = gretl_model_get_double(pmod, "lmin");
	} else {
	    for (i=4; i<MID_STATS; i++) {
		val[i] = NADBL;
	    }
	}
    } else if (pmod->ci == ARMA) {
	key[2] = N_("Mean of innovations"); /* 22: Mean of ARMA innovations */
	val[2] = gretl_model_get_double(pmod, "mean_error");
	key[3] = N_("S.D. of innovations"); /* 22: Std. dev. of ARMA innovations */
	for (i=4; i<MID_STATS; i++) {
	    if (i < 8 || i > 11) {
		val[i] = NADBL;
	    }
	}	
    } else if (pmod->ci == LAD) {
	key[0] = N_("Median depend. var");  /* 22: Median of dependent variable */
	val[0] = gretl_model_get_double(pmod, "ymedian");
	key[2] = N_("Sum absolute resid");  /* 22: Sum of absolute residuals */
	val[2] = gretl_model_get_double(pmod, "ladsum");
	key[3] = N_("Sum squared resid");
	val[3] = pmod->ess;
	for (i=4; i<MID_STATS; i++) {
	    if (i < 8 || i > 11) {
		val[i] = NADBL;
	    }
	}
    } else if (binary_model(pmod)) {
	val[2] = val[3] = NADBL;
	val[6] = val[7] = NADBL;
	val[12] = val[13] = NADBL;
    } else if (pmod->ci == TOBIT) {
	key[2] = N_("Censored obs"); /* 22: Number of censored observations */
	val[2] = gretl_model_get_int(pmod, "censobs");
	mtab.ipos[2] = 1;
	key[3] = "sigma";
	for (i=4; i<8; i++) {
	    val[i] = NADBL;
	}
    } else if (panel_ML_model(pmod) || 
	       pmod->ci == GARCH ||
	       pmod->ci == HECKIT) {
	for (i=2; i<MID_STATS; i++) {
	    if (i < 8 || i > 11) {
		val[i] = NADBL;
	    }
	}	
    } else if (pmod->ci == MLE || ordered_model(pmod)) {
	for (i=0; i<MID_STATS; i++) {
	    if (i < 8 || i > 11) {
		val[i] = NADBL;
	    }
	}	
    } else if (pmod->ci != VAR && pmod->aux != AUX_VECM && 
	       !na(pmod->rho) && gretl_model_get_int(pmod, "ldepvar")) {
	double h = durbins_h(pmod);

	if (!na(h)) {
	    key[13] = (tex)? N_("Durbin's $h$") : N_("Durbin's h");
	    val[13] = h;
	}
    } else if (pmod->ci == INTREG) {
	double x = gretl_model_get_double(pmod, "overall_test");

	if (na(x)) {
	    val[6] = val[7] = NADBL;
	} else {
	    sprintf(chistr, "%s(%d)", _("Chi-square"), pmod->dfn);
	    key[6] = chistr;  
	    val[6] = gretl_model_get_double(pmod, "overall_test");
	    key[7] = N_("p-value");  
	    val[7] = chisq_cdf_comp(pmod->dfn, val[6]);
	}
	for (i=0; i<MID_STATS; i++) {
	    if (i < 6 || i > 11) {
		val[i] = NADBL;
	    }
	}
    } 

    if (code == MIDDLE_TRANSFORM) {
	/* transformed data: don't print mean, s.d. of dep. var. */
	val[0] = val[1] = NADBL;
    } else if (code == MIDDLE_ORIG) {
	/* print a limited range of stats */
	val[2] = gretl_model_get_double(pmod, "ess_orig");
	val[3] = gretl_model_get_double(pmod, "sigma_orig");
	for (i=4; i<MID_STATS; i++) {
	    val[i] = NADBL;
	}
    }

    /* start the table */
    if (tex) {
	pputs(prn, "\\vspace{1ex}\n");
	if (mtab.multi) {
	    pputs(prn, "\\begin{tabular}{lr}\n");
	} else {
	    pputs(prn, "\\begin{tabular}{lrlr}\n");
	}
    } else if (rtf) {
	if (mtab.multi) {
	    pprintf(prn, "\\par\n{%s", RTF_MULTI_ROW);
	} else {
	    pprintf(prn, "\\par\n{%s", RTF_MT_ROW);
	}
    } else if (csv) {
	mtab.d = prn_delim(prn);
    }

    mtab.key = key;
    mtab.val = val;
    strcpy(mtab.txt_fmt, TXT_MT_FMT);

#ifdef ENABLE_NLS
    if (plain_format(prn) && mtab.nls) {
	set_mtab_string_width(&mtab);
    }
#endif

    /* print the various statistics */
    for (i=0, j=0; i<7; i++, j+=2) {
	if (!na(val[j])) {
	    middle_table_row(&mtab, j, prn);
	}
    }

    if (tex) {
	pputs(prn, "\\end{tabular}\n\n");
    } else if (rtf) {
	pputs(prn, "}\n\n"); /* close RTF table */
    } 

    if (plain_format(prn) && strcmp(note, _(note)) &&
	!string_is_blank(_(note))) {
	pputs(prn, _(note));
	pputc(prn, '\n');
    }

    if (!tex && !rtf) {
	pputc(prn, '\n');
    }
} 

static void print_model_iter_info (const MODEL *pmod, PRN *prn)
{
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

static void set_csv_delim (PRN *prn)
{
    char test[4];

    sprintf(test, "%.1f", 1.0);

    if (test[1] == ',') {
	gretl_print_set_delim(prn, '\t');
    } else {
	gretl_print_set_delim(prn, ',');
    }
}

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
	set_csv_delim(prn);
    }

    if (plain_format(prn)) {
	print_model_iter_info(pmod, prn);
    } 

    model_format_start(prn);

    print_model_heading(pmod, pdinfo, opt, prn);

    if (plain_format(prn)) {
	gotnan = plain_print_coefficients(pmod, pdinfo, prn);
    } else {
	gotnan = alt_print_coefficients(pmod, pdinfo, prn);
    }

    print_coeff_table_end(pmod, prn);

    if (pmod->aux == AUX_ARCH || pmod->aux == AUX_ADF || 
	pmod->aux == AUX_RESET || pmod->aux == AUX_DF || 
	pmod->aux == AUX_KPSS) {
	goto close_format;
    }

    /* other auxiliary regressions */

    if (pmod->aux == AUX_HET_1) {
	rsqline(pmod, prn);
	print_HET_1_results(pmod, prn);
	goto close_format;
    } else if (pmod->aux == AUX_WHITE) { 
	rsqline(pmod, prn);
	print_whites_results(pmod, prn);
	goto close_format;
    } else if (pmod->aux == AUX_BP) { 
	rssline(pmod, prn);
	print_bp_results(pmod, prn);
	goto close_format;
    } else if (pmod->aux == AUX_SQ || 
	       pmod->aux == AUX_LOG || 
	       pmod->aux == AUX_AR) {
	rsqline(pmod, prn);
	goto close_format;
    } 

    if (opt & OPT_S) {
	/* --simple-print */
	goto close_format;
    }

    if (weighted_model(pmod)) {
	alternate_stats_message(WTD_STATS, prn);
	print_middle_table(pmod, prn, MIDDLE_TRANSFORM);
	alternate_stats_message(ORIG_STATS, prn);
	print_middle_table(pmod, prn, MIDDLE_ORIG);
	goto pval_max;
    }

    if (pmod->ci == LOGISTIC) {
	alternate_stats_message(TRANSFORM_STATS, prn);
	print_middle_table(pmod, prn, MIDDLE_TRANSFORM);
	alternate_stats_message(ORIG_STATS, prn);
	print_middle_table(pmod, prn, MIDDLE_ORIG);
	goto pval_max;
    }

    if (pmod->ci == AR || pmod->ci == AR1) {
	alternate_stats_message(RHODIFF_STATS, prn);
    }

    /* print table of model stats */
    if (pmod->ci != GMM) {
	print_middle_table(pmod, prn, MIDDLE_REGULAR);
    }

    /* additional stats/info for some cases */
    if (pmod->aux == AUX_SYS) {
	if (liml_equation(pmod)) {
	    print_liml_equation_data(pmod, prn);
	}
    } else if (pmod->ci == ARMA) {
	print_arma_roots(pmod, prn);
    } else if (pmod->ci == GARCH) {
	garch_variance_line(pmod, prn);
    } else if (pmod->ci == HECKIT) {
	print_heckit_stats(pmod, prn);
    } else if (random_effects_model(pmod)) {
	panel_variance_lines(pmod, prn);
    } else if (pmod->ci == GMM) {
	print_GMM_stats(pmod, prn);
    } else if (pmod->ci == ARBOND) {
	print_DPD_stats(pmod, prn);
    } else if (binary) {
	print_binary_statistics(pmod, pdinfo, prn);
    } else if (pmod->ci == TSLS && plain_format(prn)) {
	addconst_message(pmod, prn);
	r_squared_message(pmod, prn);
    } else if (pmod->ci == INTREG) {
	print_intreg_info(pmod, prn);
    }

    /* FIXME alternate R^2 measures (within, centered) */

    if (plain_format(prn)) {
	maybe_print_jll(pmod, 0, prn);
    }

 pval_max:

    if (plain_format(prn) && pmod->ci != MLE && pmod->ci != PANEL &&
	pmod->ci != ARMA && pmod->ci != NLS && pmod->ci != GMM &&
	pmod->ci != POISSON && pmod->ci != TOBIT && pmod->ci != LAD &&
	pmod->ci != HECKIT && pmod->ci != ARBOND && pmod->ci != GARCH &&
	!ordered_model(pmod) && !pmod->aux) {
	pval_max_line(pmod, pdinfo, prn);
    }

 close_format:

    if (opt & OPT_V) {
	ols_print_anova(pmod, prn);
    }

    if (opt & OPT_O) {
	outcovmx(pmod, pdinfo, prn);
    }

    if (any_tests(pmod) && !(opt & OPT_S)) {
	print_model_tests(pmod, prn);
    }

    model_format_end(prn);
    
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

/* used only for "alternate" model-output formats */

static int 
prepare_model_coeff (const MODEL *pmod, const DATAINFO *pdinfo,
		     int i, int adfnum, model_coeff *mc, PRN *prn)
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
	if (i == adfnum) {
	    mc->pval = gretl_model_get_double(pmod, "dfpval");
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

    if (pmod->ci == MPOLS) {
	mc->multi = 1;
    }

    return gotnan;
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
    if (!na(mc->lo)) {
	pputs(prn, RTF_INTVL_ROW);
    } else {
	pputs(prn, RTF_COEFF_ROW);
    }

    pprintf(prn, "\\ql %s\\cell", mc->name);

    if (na(mc->b)) {
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
    } else {
	rtf_print_double(mc->b, prn);
    }

    if (!na(mc->lo) && !na(mc->hi)) {
	rtf_print_double(mc->lo, prn);
	rtf_print_double(mc->hi, prn);
	goto rtf_finish;
    }

    if (na(mc->se)) {
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	pprintf(prn, " \\qc %s\\cell", I_("undefined"));
	goto rtf_finish;
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

 rtf_finish:

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

    if (!na(mc->lo) && !na(mc->hi)) {
	pprintf(prn, "%c%.15g", d, mc->lo);
	pprintf(prn, "%c%.15g\n", d, mc->hi);
	return;
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

static void alt_print_coeff (const model_coeff *mc, PRN *prn)
{
    if (rtf_format(prn)) {
	rtf_print_coeff(mc, prn);
    } else if (tex_format(prn)) {
	tex_print_coeff(mc, prn);
    } else if (csv_format(prn)) {
	csv_print_coeff(mc, prn);
    }
}

static void print_coeff_separator (const char *s, PRN *prn)
{
    if (tex_format(prn)) {
	/* FIXME */
	pputs(prn, "\\\\ \n");
    } else {
	if (s != NULL && *s != '\0') {
	    pputc(prn, '\n');
	    /* FIXME RTF */
	    print_centered((rtf_format(prn))? I_(s) : _(s), 78, prn);
	    pputc(prn, '\n');
	}
	pputc(prn, '\n');
    }
}

static int 
print_rq_sequence (const MODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    gretl_vector *tauvec = gretl_model_get_data(pmod, "rq_tauvec");
    gretl_matrix *B = gretl_model_get_data(pmod, "rq_sequence");
    double tau, bi, se = NADBL;
    double blo = 0, bhi = 0;
    char test[16];
    int n, bcols, taulen = 5;
    int namelen = 8;
    int offset;
    int ntau, i, j, k = 0;
    const char *headings[] = {
	N_("tau"),
	N_("coefficient"),
	N_("lower"),
	N_("upper"),
	N_("std. error"),
	N_("t-ratio")
    };
    char *head;

    if (tauvec == NULL || B == NULL) {
	return E_DATA;
    }

    ntau = gretl_vector_get_length(tauvec);
    bcols = gretl_matrix_cols(B);

    for (i=2; i<=pmod->list[0]; i++) {
	n = char_len(pdinfo->varname[pmod->list[i]]);
	if (n > namelen) {
	    namelen = n;
	}
    }

    for (j=0; j<ntau; j++) {
	sprintf(test, "%g", gretl_vector_get(tauvec, j));
	n = strlen(test);
	if (n > taulen) {
	    taulen = n;
	}
    }

    /* headings row */

    bufspace(namelen + 4, prn);
    pprintf(prn, " %-*s", taulen, _(headings[0]));
    for (i=1; i<4; i++) {
	j = (i == 1)? 1 : (bcols == 3)? i : i + 2;
	head = _(headings[j]);
	n = char_len(head);
	offset = (i == 1)? 13 - n : 12 - n;
	bufspace(offset, prn);
	pputs(prn, head);
	pputc(prn, ' ');
    }
    pputc(prn, '\n');

    /* separator row */
    pputs(prn, "  ");
    n = namelen + 2 + taulen + 1 + 13 * 3;
    for (i=0; i<n; i++) {
	pputc(prn, '-');
    }
    pputc(prn, '\n');

    for (i=0; i<pmod->ncoeff; i++) {
	pprintf(prn, "  %-*s  ", namelen, pdinfo->varname[pmod->list[i+2]]);
	for (j=0; j<ntau; j++) {
	    tau = gretl_vector_get(tauvec, j);
	    bi  = gretl_matrix_get(B, k, 0);
	    if (bcols == 3) {
		blo = gretl_matrix_get(B, k, 1);
		bhi = gretl_matrix_get(B, k, 2);
	    } else {
		se = gretl_matrix_get(B, k, 1);
	    }
	    if (j > 0) {
		bufspace(namelen + 4, prn);
	    }
	    if (bcols == 3) {
		pprintf(prn, "%.*f  %#12.6g %#12.6g %#12.6g\n", taulen - 2,
			tau, bi, blo, bhi);
	    } else {
		pprintf(prn, "%.*f  %#12.6g %#12.6g %#12.6g\n", taulen - 2,
			tau, bi, se, bi / se);
	    }
	    k++;
	}
	if (i < pmod->ncoeff - 1) {
	    pputc(prn, '\n');
	}
    }

    return 0;
}

#if 0 /* not ready */

static int 
alt_print_rq_sequence (const MODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    gretl_vector *tauvec = gretl_model_get_data(pmod, "rq_tauvec");
    gretl_matrix *B = gretl_model_get_data(pmod, "rq_sequence");
    double tau, bi, se = NADBL;
    double blo = 0, bhi = 0;
    char test[16];
    int n, bcols;
    int offset;
    int ntau, i, j, k = 0;

    if (tauvec == NULL || B == NULL) {
	return E_DATA;
    }

    ntau = gretl_vector_get_length(tauvec);
    bcols = gretl_matrix_cols(B);


    for (i=0; i<pmod->ncoeff; i++) {
	pprintf(prn, "  %-*s  ", namelen, pdinfo->varname[pmod->list[i+2]]);
	for (j=0; j<ntau; j++) {
	    tau = gretl_vector_get(tauvec, j);
	    bi  = gretl_matrix_get(B, k, 0);
	    if (bcols == 3) {
		blo = gretl_matrix_get(B, k, 1);
		bhi = gretl_matrix_get(B, k, 2);
	    } else {
		se = gretl_matrix_get(B, k, 1);
	    }
	    if (bcols == 3) {
		pprintf(prn, "%.*f  %#12.6g %#12.6g %#12.6g\n", taulen - 2,
			tau, bi, blo, bhi);
	    } else {
		pprintf(prn, "%.*f  %#12.6g %#12.6g %#12.6g\n", taulen - 2,
			tau, bi, se, bi / se);
	    }
	    k++;
	}
	if (i < pmod->ncoeff - 1) {
	    pputc(prn, '\n');
	}
    }

    return 0;
}

#endif

struct printval {
    char s[36];
    int lw, rw;
    double x;
};

static void get_number_dims (struct printval *v, int *lmax, int *rmax)
{
    char tmp[36];
    const char *s = v->s;
    char *p;
    int i, n;

    while (*s == ' ') s++;
    n = strlen(s);
    for (i=n-1; i>0 && s[i] == ' '; i--) {
	n--;
    }

    *tmp = '\0';
    strncat(tmp, s, n);

    p = strchr(tmp, '.');
    if (p == NULL) {
	p = strchr(tmp, ',');
    }

    if (p == NULL) {
	v->lw = char_len(tmp);
	v->rw = 0;
    } else {
	n = char_len(tmp);
	*p = '\0';
	v->lw = char_len(tmp);
	v->rw = n - v->lw;
    }

    if (v->lw > *lmax) *lmax = v->lw;
    if (v->rw > *rmax) *rmax = v->rw;
}

static void print_padded_head (const char *s, int w, PRN *prn)
{
    int hlen = char_len(s);
    int offset = (w - hlen) / 2;
    int pad = w - hlen - offset;

    bufspace(offset, prn);
    pputs(prn, s);
    bufspace(pad, prn);    
}

static void print_padded_value (struct printval *val, int w, 
				int lmax, int addoff, PRN *prn)
{
    int offset = lmax - val->lw + addoff;
    int pad = w - (lmax + val->rw + addoff);

    bufspace(offset, prn);
    pputs(prn, val->s);
    bufspace(pad, prn);    
}

static struct printval **allocate_printvals (int n, int m)
{
    struct printval **vals;
    int i, j;

    vals = malloc(n * sizeof *vals);

    for (i=0; i<n && vals != NULL; i++) {
	vals[i] = malloc(m * sizeof **vals);
	if (vals[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(vals[j]);
	    }
	    free(vals);
	    vals = NULL;
	}
    }

    return vals;
}

static int get_ar_data (const MODEL *pmod, double **pb,
			double **pse, int *pk, int *dfd)
{
    double *b, *se;
    int i, k, err = 0;

    if (pmod->arinfo == NULL || 
	pmod->arinfo->arlist == NULL ||
	pmod->arinfo->rho == NULL ||
	pmod->arinfo->sderr == NULL) {
	return E_DATA;
    }

    k = pmod->arinfo->arlist[0];

    if (k > 1) {
	*dfd = pmod->dfd + (pmod->ncoeff - k);
    } 

    b = malloc((pmod->ncoeff + k) * sizeof *b);
    se = malloc((pmod->ncoeff + k) * sizeof *se);

    if (b == NULL || se == NULL) {
	free(b);
	free(se);
	return E_ALLOC;
    }

    for (i=0; i<pmod->ncoeff; i++) {
	b[i] = pmod->coeff[i];
	se[i] = pmod->sderr[i];
    }

    for (i=0; i<k; i++) {
	b[pmod->ncoeff + i] = pmod->arinfo->rho[i];
	se[pmod->ncoeff + i] = pmod->arinfo->sderr[i];
    }

    *pb = b;
    *pse = se;
    *pk = k;

    return err;
}

static int get_arch_data (const MODEL *pmod, double **pb,
			  double **pse, int *pk)
{
    double *ab = gretl_model_get_data(pmod, "arch_coeff");
    double *ase = gretl_model_get_data(pmod, "arch_sderr");
    int order = gretl_model_get_int(pmod, "arch_order");
    double *b, *se;
    int err = 0;

    if (ab != NULL && ase != NULL && order > 0) {
	int i, k = order + 1;

	b = malloc((pmod->ncoeff + k) * sizeof *b);
	se = malloc((pmod->ncoeff + k) * sizeof *se);

	if (b == NULL || se == NULL) {
	    free(b);
	    free(se);
	    return E_ALLOC;
	}

	for (i=0; i<pmod->ncoeff; i++) {
	    b[i] = pmod->coeff[i];
	    se[i] = pmod->sderr[i];
	}

	for (i=0; i<k; i++) {
	    b[pmod->ncoeff + i] = ab[i];
	    se[pmod->ncoeff + i] = ase[i];
	}

	*pb = b;
	*pse = se;
	*pk = k;
    } else {
	err = E_DATA;
    }

    return err;
}

static void figure_colsep (int namelen, int ncols, int *w,
			   int *colsep)
{
    int j, n = namelen + ncols * *colsep;

    for (j=0; j<ncols; j++) {
	n += w[j];
    }   
    n = 68 - (n + 2);
    if (n > ncols) {
	*colsep += 1;
    }
}

static void print_sep_row (int namelen, int ncols, int *w,
			   int colsep, PRN *prn)
{
    int j, n;

    pputc(prn, '\n');
    bufspace(2, prn);
    n = namelen + ncols * colsep;
    for (j=0; j<ncols; j++) {
	n += w[j];
    }
    for (j=0; j<n; j++) {
	pputc(prn, '-');
    }
    pputc(prn, '\n');
}

static void put_asts (double pv, PRN *prn)
{
    if (pv < 0.01) {
	pputs(prn, " ***");
    } else if (pv < 0.05) {
	pputs(prn, " **");
    } else if (pv < 0.10) {
	pputs(prn, " *");
    }
}

static int alt_print_aux_coeffs (const double *b, const double *se,
				 const char **names, int nc, int df, 
				 int ci, PRN *prn)
{
    model_coeff mc;
    int i;

    for (i=0; i<nc; i++) {
	if (xna(b[i])) {
	    return E_NAN;
	}
    }

    alt_print_coeff_table_start(NULL, ci, prn);

    model_coeff_init(&mc);

    /* print row values */
    for (i=0; i<nc; i++) {
	mc.b = b[i];
	mc.se = se[i];
	if (na(se[i]) || se[i] <= 0) {
	    mc.tval = mc.pval = NADBL;
	} else {
	    mc.tval = b[i] / se[i];
	    mc.pval = coeff_pval(ci, mc.tval, df);
	}
	if (tex_format(prn)) {
	    tex_escape_special(mc.name, names[i]);
	} else {
	    *mc.name = '\0';
	    strncat(mc.name, names[i], MC_NAMELEN - 1);
	}
	alt_print_coeff(&mc, prn);
    }

    print_coeff_table_end(NULL, prn);

    return 0;
}

static int plain_print_aux_coeffs (const double *b,
				   const double *se,
				   const char **names,
				   int nc, int df, int ci,
				   PRN *prn)
{
    struct printval **vals, *vij;
    const char *headings[] = { 
	N_("coefficient"), 
	N_("std. error"), 
	N_("t-ratio"),
	N_("p-value")
    };
    const char *head;
    int lmax[4] = {0};
    int rmax[4] = {0};
    int w[4], addoff[4] = {0};
    int hlen;
    double tval, pval;
    int namelen = 0;
    int colsep = 2;
    int n, d, i, j;
    int err = 0;

    vals = allocate_printvals(nc, 4);
    if (vals == NULL) {
	return E_ALLOC;
    }

    if (ci == MODPRINT) {
	headings[2] = N_("z-stat");
    }

    for (i=0; i<nc; i++) {
	if (xna(b[i])) {
	    err = E_NAN;
	    goto bailout;
	}
	n = char_len(names[i]);
	if (n > namelen) {
	    namelen = n;
	}
	if (na(se[i]) || se[i] <= 0) {
	    tval = pval = NADBL;
	} else {
	    tval = b[i] / se[i];
	    pval = coeff_pval(ci, tval, df);
	}
	for (j=0; j<4; j++) {
	    if (j < 2) {
		/* coeff, standard error */
		d = GRETL_DIGITS;
		vals[i][j].x = (j==0)? b[i] : se[i];
	    } else if (j == 2) {
		/* t-ratio */
		d = 4;
		vals[i][j].x = tval;
	    } else {
		/* p-value */
		d = -4;
		vals[i][j].x = pval;
	    }
	    gretl_sprint_fullwidth_double(vals[i][j].x, d, vals[i][j].s, prn);
	    get_number_dims(&vals[i][j], &lmax[j], &rmax[j]);
	}
    }

    if (namelen < 8) {
	namelen = 8;
    }

    /* figure appropriate column separation */

    for (j=0; j<4; j++) {
	w[j] = lmax[j] + rmax[j];
	head = _(headings[j]);
	hlen = char_len(head);
	if (hlen > w[j]) {
	    addoff[j] = (hlen - w[j]) / 2;
	    w[j] = hlen;
	}
    }

    figure_colsep(namelen, 4, w, &colsep);

    /* print headings */

    bufspace(namelen + 2 + colsep, prn);
    for (j=0; j<4; j++) {
	print_padded_head(_(headings[j]), w[j], prn);
	if (j < 3) {
	    bufspace(colsep, prn);
	}
    }

    /* separator row */
    print_sep_row(namelen, 4, w, colsep, prn);

    /* print row values */

    for (i=0; i<nc; i++) {
	pprintf(prn, "  %-*s", namelen, names[i]);
	bufspace(colsep, prn);
	for (j=0; j<4; j++) {
	    vij = &vals[i][j];
	    print_padded_value(vij, w[j], lmax[j], addoff[j], prn);
	    if (j < 3) {
		bufspace(colsep, prn);
	    } else if (!na(vij->x)) {
		put_asts(vij->x, prn);
	    }	
	}
	pputc(prn, '\n');
    }

 bailout:

    for (i=0; i<nc; i++) {
	free(vals[i]);
    }
    free(vals);

    return err;
}

/* Called by external functions that want to print a coefficient
   table that does not reside in a MODEL.
*/

int print_coeffs (const double *b,
		  const double *se,
		  const char **names,
		  int nc, int df, int ci,
		  PRN *prn)
{
    if (plain_format(prn)) {
	return plain_print_aux_coeffs(b, se, names, nc, df, ci, prn);
    } else {
	return alt_print_aux_coeffs(b, se, names, nc, df, ci, prn);
    }
}

static int plain_print_mp_coeffs (const MODEL *pmod, 
				  const DATAINFO *pdinfo, 
				  PRN *prn)
{
    struct printval **vals;
    const char *headings[] = { 
	N_("coefficient"), 
	N_("std. error") 
    };
    const double *b = pmod->coeff;
    const double *se = pmod->sderr;
    char **names = NULL;
    const char *head;
    int lmax[2] = {0};
    int rmax[2] = {0};
    int w[2], addoff[2] = {0};
    int hlen;
    int n, nc = pmod->ncoeff;
    int namelen = 0;
    int colsep = 2;
    int i, j, minus;
    int err = 0;

    vals = allocate_printvals(nc, 2);
    if (vals == NULL) {
	return E_ALLOC;
    }

    names = strings_array_new_with_length(nc, 32);
    if (names == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    minus = (gretl_print_supports_utf(prn))? MINUS_UTF : MINUS_HYPHEN;

    for (i=0; i<nc; i++) {
	if (xna(b[i])) {
	    err = E_NAN;
	    goto bailout;
	}
	gretl_model_get_param_name(pmod, pdinfo, i, names[i]);
	n = char_len(names[i]);
	if (n > namelen) {
	    namelen = n;
	}
	for (j=0; j<2; j++) {
	    vals[i][j].x = (j==0)? b[i] : se[i];
	    print_fifteen(vals[i][j].s, vals[i][j].x, minus);
	    vals[i][j].lw = lmax[j] = 2;
	    vals[i][j].rw = rmax[j] = 19;
	}
    }

    if (namelen < 8) {
	namelen = 8;
    }

    /* figure appropriate column separation */

    for (j=0; j<2; j++) {
	w[j] = lmax[j] + rmax[j];
	head = _(headings[j]);
	hlen = char_len(head);
	if (hlen > w[j]) {
	    addoff[j] = (hlen - w[j]) / 2;
	    w[j] = hlen;
	}
    }

    figure_colsep(namelen, 2, w, &colsep);

    /* print headings */

    bufspace(namelen + 2 + colsep, prn);
    for (j=0; j<2; j++) {
	print_padded_head(_(headings[j]), w[j], prn);
	if (j < 3) {
	    bufspace(colsep, prn);
	}
    }

    /* separator row */
    print_sep_row(namelen, 2, w, colsep, prn);

    /* print row values */

    for (i=0; i<nc; i++) {
	pprintf(prn, "  %-*s", namelen, names[i]);
	bufspace(colsep, prn);
	for (j=0; j<2; j++) {
	    print_padded_value(&vals[i][j], w[j], lmax[j], addoff[j], prn);
	    if (j == 0) {
		bufspace(colsep, prn);
	    } 
	}
	pputc(prn, '\n');
    }

 bailout:

    for (i=0; i<nc; i++) {
	free(vals[i]);
    }
    free(vals);

    free_strings_array(names, nc);

    return err;
}

static char *get_col_heading (const char **S, int j, int slopes,
			      int intervals)
{
    if (j == 3 && slopes) {
	return _(S[j+1]);
    } else if (intervals && (j == 1 || j == 2)) {
	return _(S[j+4]);
    } else {
	return _(S[j]);
    }
}

static void 
print_poisson_offset (const MODEL *pmod, const DATAINFO *pdinfo, 
		      struct printval *val, int namelen, 
		      int colsep, int w, int lmax,
		      int addoff, PRN *prn)
{
    int offvar = gretl_model_get_int(pmod, "offset_var");

    if (offvar > 0) {
	char name[24];
	int n;

	sprintf(name, "log(%s)", pdinfo->varname[offvar]);
	n = strlen(name);
	pprintf(prn, "\n  %-*s", namelen, name);
	if (n > namelen) {
	    colsep -= n - namelen;
	}
	if (colsep > 0) {
	    bufspace(colsep, prn);
	}
	sprintf(val->s, "%.1f", 1.0);
	val->lw = val->rw = 1;
	print_padded_value(val, w, lmax, addoff, prn);
	pputc(prn, '\n');
    }
}

static void print_ar_sum (const MODEL *pmod, PRN *prn)
{
    if (pmod->arinfo != NULL && pmod->arinfo->arlist[0] > 1) {
	double arsum = 0.0;
	int i;

	for (i=0; i<pmod->arinfo->arlist[0]; i++) {
	    arsum += pmod->arinfo->rho[i];
	}
	pprintf(prn, "\n%s = %#g\n", _("Sum of AR coefficients"), arsum);
    }
}

static int plain_print_coeffs (const MODEL *pmod, 
			       const DATAINFO *pdinfo, 
			       PRN *prn)
{
    struct printval **vals, *vij;
    const char *headings[] = { 
	N_("coefficient"), 
	N_("std. error"), 
	N_("t-ratio"),
	N_("p-value"),
	N_("slope"),
	N_("lower"),
	N_("upper")
    };
    const double *b = pmod->coeff;
    const double *se = pmod->sderr;
    double *slopes = NULL;
    char **names = NULL;
    const char *head;
    const char *sepstr = NULL;
    double *xb = NULL;
    double *xse = NULL;
    int seppos = -1;
    int lmax[4] = {0};
    int rmax[4] = {0};
    int w[4], addoff[4] = {0};
    int hlen;
    double tval, pval = 0.0;
    int n, d, nc = pmod->ncoeff;
    int dfd = pmod->dfd;
    int show_slope, adfnum = -1;
    int intervals = 0;
    int namelen = 0;
    int colsep = 2;
    int ncols = 4;
    int i, j;
    int err = 0;

    if (pmod->ci == AR || pmod->ci == ARCH) {
	int k = 0;

	if (pmod->ci == AR) {
	    err = get_ar_data(pmod, &xb, &xse, &k, &dfd);
	} else {
	    err = get_arch_data(pmod, &xb, &xse, &k);
	}
	if (!err) {
	    b = xb;
	    se = xse;
	    nc += k;
	}	
    } else if (pmod->ci == LAD) {
	gretl_matrix *m = gretl_model_get_data(pmod, "coeff_intervals");

	if (m != NULL) {
	    se = m->val;
	    intervals = 1;
	    ncols = 3;
	}
    } else if (nonlin_model(pmod)) {
	headings[0] = N_("estimate");
    }

    if (err) {
	return err;
    }

    vals = allocate_printvals(nc, ncols);
    if (vals == NULL) {
	return E_ALLOC;
    }

    names = strings_array_new_with_length(nc, 32);
    if (names == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    show_slope = binary_model(pmod) && !gretl_model_get_int(pmod, "show-pvals");
    if (show_slope) {
	slopes = gretl_model_get_data(pmod, "slopes");
    }

    if (pmod->aux == AUX_ADF || pmod->aux == AUX_DF) {
	adfnum = gretl_model_get_int(pmod, "dfnum");
    }

    gretl_model_get_coeff_separator(pmod, &sepstr, &seppos);
    if (seppos == -1) {
	if (pmod->ci == GARCH && pmod->list[0] > 4) {
	    seppos = pmod->list[0] - 4;
	} else if (pmod->ci == AR || pmod->ci == ARCH) {
	    seppos = pmod->ncoeff;
	}
    }

    for (i=0; i<nc; i++) {
	if (xna(b[i])) {
	    err = E_NAN;
	    goto bailout;
	}
	gretl_model_get_param_name(pmod, pdinfo, i, names[i]);
	n = char_len(names[i]);
	if (n > namelen) {
	    namelen = n;
	}
	if (na(se[i]) || se[i] <= 0) {
	    tval = pval = NADBL;
	} else {
	    tval = b[i] / se[i];
	    if (slopes != NULL) {
		/* last column: actually slope at mean */
		pval = (pmod->list[i+2] == 0)? 0 : slopes[i];
	    } else if (i == adfnum) {
		/* special Dickey-Fuller p-value */
		pval = gretl_model_get_double(pmod, "dfpval");
	    } else if (!intervals) {
		/* regular p-value */
		pval = coeff_pval(pmod->ci, tval, dfd);
	    } 
	}
	for (j=0; j<ncols; j++) {
	    if (j < 2) {
		/* coeff, standard error or lower c.i. limit */
		d = GRETL_DIGITS;
		vals[i][j].x = (j==0)? b[i] : se[i];
	    } else if (j == 2) {
		/* t-ratio or upper c.i. limit */
		if (intervals) {
		    d = GRETL_DIGITS;
		    vals[i][j].x = se[i + nc];
		} else {
		    d = 4;
		    vals[i][j].x = tval;
		}
	    } else {
		/* p-value or slope */
		d = (show_slope)? 6 : -4;
		vals[i][j].x = pval;
	    }
	    if (show_slope && j == 3 && pmod->list[i+2] == 0) {
		/* don't show 'slope' for constant */
		vals[i][j].s[0] = '\0';
		vals[i][j].lw = vals[i][j].rw = 0;
	    } else {
		gretl_sprint_fullwidth_double(vals[i][j].x, d, vals[i][j].s, prn);
		get_number_dims(&vals[i][j], &lmax[j], &rmax[j]);
	    }
	}
    }

    if (namelen < 8) {
	namelen = 8;
    }

    /* figure appropriate column separation */

    for (j=0; j<ncols; j++) {
	w[j] = lmax[j] + rmax[j];
	head = get_col_heading(headings, j, show_slope, intervals);
	hlen = char_len(head);
	if (hlen > w[j]) {
	    addoff[j] = (hlen - w[j]) / 2;
	    w[j] = hlen;
	}
    }

    figure_colsep(namelen, ncols, w, &colsep);

    /* print headings */

    bufspace(namelen + 2 + colsep, prn);
    for (j=0; j<ncols; j++) {
	head = get_col_heading(headings, j, show_slope, intervals);
	print_padded_head(head, w[j], prn);
	if (j < ncols - 1) {
	    bufspace(colsep, prn);
	}
    }

    /* separator row */
    print_sep_row(namelen, ncols, w, colsep, prn);

    /* print row values */

    for (i=0; i<nc; i++) {
	if (i == seppos) {
	    print_coeff_separator(sepstr, prn);
	}
	pprintf(prn, "  %-*s", namelen, names[i]);
	bufspace(colsep, prn);
	for (j=0; j<ncols; j++) {
	    vij = &vals[i][j];
	    print_padded_value(vij, w[j], lmax[j], addoff[j], prn);
	    if (j < ncols - 1) {
		bufspace(colsep, prn);
	    } else if (!show_slope && !na(vij->x)) {
		put_asts(vij->x, prn);
	    }	
	}
	pputc(prn, '\n');
    }

    if (pmod->ci == AR) {
	print_ar_sum(pmod, prn);
    } else if (pmod->ci == POISSON) {
	print_poisson_offset(pmod, pdinfo, &vals[0][0], namelen, 
			     colsep, w[0], lmax[0], addoff[0],
			     prn);
    }

 bailout:

    for (i=0; i<nc; i++) {
	free(vals[i]);
    }
    free(vals);

    free_strings_array(names, nc);
    free(xb);
    free(xse);

    return err;
}

static int 
plain_print_coefficients (const MODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    if (pmod->ncoeff == 0) {
	return 0;
    } else if (gretl_model_get_data(pmod, "rq_sequence") != NULL) {
	return print_rq_sequence(pmod, pdinfo, prn);
    } else if (pmod->ci == MPOLS) {
	return plain_print_mp_coeffs(pmod, pdinfo, prn);
    } else {
	return plain_print_coeffs(pmod, pdinfo, prn);
    }
}

static void alt_print_arch_terms (const MODEL *pmod, PRN *prn)
{
    double *a = gretl_model_get_data(pmod, "arch_coeff");
    double *se = gretl_model_get_data(pmod, "arch_sderr");
    int order = gretl_model_get_int(pmod, "arch_order");

    if (a != NULL && se != NULL && order > 0) {
	model_coeff mc;
	int i;

	gretl_prn_newline(prn);

	for (i=0; i<=order; i++) {
	    model_coeff_init(&mc);
	    mc.b = a[i];
	    mc.se = se[i];
	    mc.tval = a[i] / se[i];
	    mc.pval = student_pvalue_2(pmod->nobs - (order + 1), mc.tval);

	    if (tex_format(prn)) {
		sprintf(mc.name, "$\\alpha_%d$", i);
	    } else {
		sprintf(mc.name, "alpha(%d)", i);
	    }

	    alt_print_coeff(&mc, prn);
	}
    }
}

static void 
alt_print_poisson_offset (const MODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    int offvar = gretl_model_get_int(pmod, "offset_var");

    if (offvar > 0) {
	char name[24];

	sprintf(name, "log(%s)", pdinfo->varname[offvar]);

	if (plain_format(prn)) {
	    pprintf(prn, "\n  %-13s         1.0\n", name);
	} else if (rtf_format(prn)) {
	    pputs(prn, RTF_COEFF_ROW);
	    pprintf(prn, "\\ql %s\\cell\\qc 1.0\\cell", name);
	    pputs(prn, "\\qc \\cell\\qc \\cell \\qc \\cell \\intbl \\row\n");
	} else if (tex_format(prn)) {
	    char tmp[32];

	    tex_escape(tmp, name);
	    pprintf(prn, "{\\rm %s} & \\multicolumn{1}{c}{1.0} \\\\\n", tmp);
	}
    }
}

static int 
alt_print_coefficients (const MODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *intervals = NULL;
    const char *sepstr = NULL;
    int seppos = -1;
    model_coeff mc;
    int adfnum = -1;
    int nc = pmod->ncoeff;
    int i, err = 0, gotnan = 0;

    if (gretl_model_get_data(pmod, "rq_sequence") != NULL) {
	pputs(prn, "Sorry, not implemented yet!\n");
	return 1;
    }

    alt_print_coeff_table_start(pmod, pmod->ci, prn);

    if (pmod->ci == PANEL) {
	nc = pmod->list[0] - 1;
    }

    gretl_model_get_coeff_separator(pmod, &sepstr, &seppos);
    if (seppos == -1 && pmod->ci == GARCH && pmod->list[0] > 4) {
	seppos = pmod->list[0] - 4;
    }

    if (pmod->ci == LAD) {
	intervals = gretl_model_get_data(pmod, "coeff_intervals");
    }

    if (pmod->aux == AUX_DF || pmod->aux == AUX_ADF) {
	adfnum = gretl_model_get_int(pmod, "dfnum");
    }

    for (i=0; i<nc; i++) {

	err = prepare_model_coeff(pmod, pdinfo, i, adfnum, &mc, prn);
	if (err) gotnan = 1;

	if (i == seppos) {
	    print_coeff_separator(sepstr, prn);
	}

	if (intervals != NULL) {
	    mc.lo = gretl_matrix_get(intervals, i, 0);
	    mc.hi = gretl_matrix_get(intervals, i, 1);
	    mc.show_pval = 0;
	}

	alt_print_coeff(&mc, prn);
    }

    if (pmod->ci == AR) {
	alt_print_rho_terms(pmod, prn);
    } else if (pmod->ci == ARCH) {
	alt_print_arch_terms(pmod, prn);
    } else if (pmod->ci == POISSON) {
	alt_print_poisson_offset(pmod, pdinfo, prn);
    }

    return gotnan;
} 

static void print_rho (const ARINFO *arinfo, int c, int dfd, PRN *prn)
{
    double tval = arinfo->rho[c] / arinfo->sderr[c];
    char ustr[32];

    if (tex_format(prn)) {
	char s1[32], s2[32], s3[32], s4[32];

	tex_rl_double(arinfo->rho[c], s1);
	tex_rl_double(arinfo->sderr[c], s2);
	tex_rl_float(tval, s3, 4);
	tex_rl_float(student_pvalue_2(dfd, tval), s4, 4);

	sprintf(ustr, "$\\hat{u}_{t-%d}$", arinfo->arlist[c+1]);

	pprintf(prn, "%s &\n"
		"  %s &\n"
		"    %s &\n"
		"      %s &\n"
		"        %s \\\\\n",  
		ustr, s1, s2, s3, s4);
    } else if (rtf_format(prn)) {
	char pvalstr[16];
	double pval;

	pputs(prn, RTF_COEFF_ROW);
	pprintf(prn, "\\ql u(-%d)\\cell", arinfo->arlist[c+1]);
	rtf_print_double(arinfo->rho[c], prn);
	rtf_print_double(arinfo->sderr[c], prn);
	pprintf(prn, " \\qc %.4f\\cell", tval);
	pval = student_pvalue_2(dfd, tval);
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

static void alt_print_rho_terms (const MODEL *pmod, PRN *prn)
{
    double xx = 0.0;
    int i, dfd;

    if (pmod->arinfo == NULL || 
	pmod->arinfo->arlist == NULL ||
	pmod->arinfo->rho == NULL ||
	pmod->arinfo->sderr == NULL) {
	return;
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
	pprintf(prn, "%s \\\\ \n", _(tag));
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
	    pprintf(prn, "%s\n%s\n", _(roots_hdr), roots_sep);
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

static void print_heckit_stats (const MODEL *pmod, PRN *prn)
{
    int totobs = gretl_model_get_int(pmod, "totobs");
    int cenobs = totobs - pmod->nobs;
    double cenpc = 100.0 * cenobs / totobs;

    ensure_vsep(prn);

    if (plain_format(prn)) {
	pprintf(prn, "%s: %d\n", _("Total observations"), totobs);
	pprintf(prn, "%s: %d (%.1f%%)\n", _("Censored observations"), 
		cenobs, cenpc);
	pprintf(prn, "%s = %.*g, ", _("sigma"), GRETL_DIGITS, pmod->sigma);
	if (na(pmod->rho)) {
	    pprintf(prn, "%s = NA\n", _("rho"));
	} else {
	    pprintf(prn, "%s = %.*g\n", _("rho"), GRETL_DIGITS, pmod->rho);
	}
	pputc(prn, '\n');
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s: %d\n", I_("Total observations"), totobs);
	pprintf(prn, RTFTAB "%s: %d (%.1f%%)\n", I_("Censored observations"), 
		cenobs, cenpc);
	pprintf(prn, RTFTAB "%s = %g\n", I_("sigma"), pmod->sigma);
	if (!na(pmod->rho)) {
	    pprintf(prn, RTFTAB "%s = %g\n", I_("rho"), pmod->rho);
	}
    } else if (tex_format(prn)) {
	char xstr[32];

	pprintf(prn, "%s: %d \\\\\n", _("Total observations"), totobs);
	pprintf(prn, "%s: %d (%.1f\\%%) \\\\\n", _("Censored observations"), 
		cenobs, cenpc);
	pprintf(prn, "$\\hat{\\sigma}$ = %s \\\\\n", 
		tex_sprint_double(pmod->sigma, xstr));
	if (!na(pmod->rho)) {
	    pprintf(prn, "$\\hat{\\rho}$ = %s \\\\\n", 
		    tex_sprint_double(pmod->rho, xstr));
	}
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

    leftlen = char_len(_("Actual")) + 3;

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
    const int *act_pred;
    int correct = -1;
    double pc_correct = NADBL;
    int i;

    act_pred = gretl_model_get_data(pmod, "discrete_act_pred");
    if (act_pred != NULL) {
	correct = act_pred[0] + act_pred[3];
	pc_correct = 100 * (double) correct / pmod->nobs;
    } 

    ensure_vsep(prn);

    if (plain_format(prn)) {
	if (correct >= 0) {
	    pprintf(prn, "%s = %d (%.1f%%)\n", 
		    _("Number of cases 'correctly predicted'"), 
		    correct, pc_correct);
	}
	pprintf(prn, "f(beta'x) %s = %.3f\n", _("at mean of independent vars"), 
		pmod->sdy);
	if (!na(model_chisq) && pmod->aux != AUX_OMIT && pmod->aux != AUX_ADD) {
	    i = pmod->ncoeff - 1;
	    if (i > 0) {
		pprintf(prn, "%s: %s(%d) = %g [%.4f]\n", 
			_("Likelihood ratio test"), _("Chi-square"), 
			i, model_chisq, chisq_cdf_comp(i, model_chisq));
	    }
	}
	pputc(prn, '\n');
	if (act_pred != NULL) {
	    plain_print_act_pred(act_pred, prn);
	}
    } else if (rtf_format(prn)) {
	pputc(prn, '\n');
	if (slopes) {
	    pprintf(prn, "\\par {\\super *}%s\n", I_("Evaluated at the mean"));
	}
	if (correct >= 0) {
	    pprintf(prn, "\\par %s = %d (%.1f%%)\n", 
		    I_("Number of cases 'correctly predicted'"), 
		    correct, pc_correct);
	}
	pprintf(prn, "\\par f(beta'x) %s = %.3f\n", I_("at mean of independent vars"), 
		pmod->sdy);
	if (pmod->aux != AUX_OMIT && pmod->aux != AUX_ADD) {
	    i = pmod->ncoeff - 1;
	    pprintf(prn, "\\par %s: %s(%d) = %g [%.4f]\n",
		    I_("Likelihood ratio test"), I_("Chi-square"), 
		    i, model_chisq, chisq_cdf_comp(i, model_chisq));
	}
	pputc(prn, '\n');
    } else if (tex_format(prn)) {
	if (slopes) {
	    pprintf(prn, "\\begin{center}\n$^*$%s\n\\end{center}\n", 
		    I_("Evaluated at the mean"));
	}

	pputs(prn, "\\vspace{1em}\n\\begin{raggedright}\n");
	if (correct >= 0) {
	    pprintf(prn, "%s = %d (%.1f percent)\\\\\n", 
		    I_("Number of cases `correctly predicted'"), 
		    correct, pc_correct);
	}
	if (pmod->aux != AUX_OMIT && pmod->aux != AUX_ADD) {
	    i = pmod->ncoeff - 1;
	    pprintf(prn, "%s: $\\chi^2(%d)$ = %.3f [%.4f]\\\\\n",
		    I_("Likelihood ratio test"), 
		    i, model_chisq, chisq_cdf_comp(i, model_chisq));
	}
	pputs(prn, "\\end{raggedright}\n");
    }
}

int ols_print_anova (const MODEL *pmod, PRN *prn)
{
    double mst, msr, mse, rss;

    if (pmod->ci != OLS || !pmod->ifc ||
	na(pmod->ess) || na(pmod->tss)) {
	return E_NOTIMP;
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
	double F = msr / mse;
	double pv = snedecor_cdf_comp(pmod->dfn, pmod->dfd, F);

	pprintf(prn, "  F(%d, %d) = %g / %g = %g ", 
		pmod->dfn, pmod->dfd, msr, mse, F);
	if (pv < .0001) {
	    pprintf(prn, "[%s %.3g]\n\n", _("p-value"), pv);
	} else {
	    pprintf(prn, "[%s %.4f]\n\n", _("p-value"), pv); 
	}
    }

    return 0;
}

static int print_user_model (const gretl_matrix *cs, 
			     const gretl_matrix *adds, 
			     const char *s, gretlopt opt,
			     PRN *prn)
{
    int ncoef = gretl_matrix_rows(cs);
    int nadd = gretl_vector_get_length(adds);
    int ntot = ncoef + nadd;
    char **names = NULL;
    char *tmp;
    const double *b, *se;
    int i, err = 0;

    /* copy the user-defined string 's' before applying strtok */
    tmp = gretl_strdup(s);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    names = malloc(ntot * sizeof *names);
    if (names == NULL) {
	free(tmp);
	return E_ALLOC;
    }

    for (i=0; i<ntot && !err; i++) {
	names[i] = strtok((i == 0)? tmp : NULL, ",");
	if (names[i] == NULL) {
	    free(names);
	    gretl_errmsg_sprintf(_("modprint: expected %d names"), ntot);
	    return E_DATA;
	}
    }

    b = cs->val;
    se = b + ncoef;

    pputc(prn, '\n');

    if (opt & OPT_C) {
	gretl_print_set_format(prn, GRETL_FORMAT_CSV);
	set_csv_delim(prn);
    } else if (opt & OPT_R) {
	gretl_print_set_format(prn, GRETL_FORMAT_RTF);
    } else if (opt & OPT_T) {
	if (opt & OPT_O) {
	    gretl_print_set_format(prn, GRETL_FORMAT_TEX | GRETL_FORMAT_DOC);
	} else {
	    gretl_print_set_format(prn, GRETL_FORMAT_TEX);
	}
    }

    model_format_start(prn);

    print_coeffs(b, se, (const char **) names, ncoef, 0, MODPRINT, prn);

    if (nadd > 0) {
	print_model_stats_table(adds->val, (const char **) names + ncoef, 
				nadd, prn);
    }

    if (plain_format(prn)) {
	pputc(prn, '\n'); 
    } 

    model_format_end(prn);

    free(names);
    free(tmp);

    return err;
}

/**
 * do_modprint:
 * @line: command line.
 * @opt: may contain %OPT_C for CSV, %OPT_R for RTF, or %OPT_T 
 * to use TeX format.  If %OPT_T is given, then %OPT_O calls for
 * a complete LaTeX document, otherwise %OPT_O is ignored.
 * @prn: gretl printer.
 *
 * Prints to @prn the coefficient table and optional additional statistics
 * for a model estimated "by hand". Mainly useful for user-written functions.
 * 
 * The string @line must contain, in order: (1) the name of a k x 2 matrix
 * containing k coefficients and k associated standard errors and (2) the
 * name of a string containing at least k comma-separated names for
 * the coefficients.  Optionally, @line may contain a third element, the 
 * name of a vector containing p additional statistics.  If this argument 
 * is supplied, then argument (2) should contain k + p comma-separated 
 * names, the additional p names to be associated with the additional 
 * statistics. 
 *
 * Returns: 0 on success, non-zero on failure.
 */

int do_modprint (const char *line, gretlopt opt, PRN *prn)
{
    gretl_matrix *coef_se = NULL;
    gretl_matrix *addstats = NULL;
    char *parnames = NULL;
    char **tmp;
    char s[MAXLEN];
    int i, err = 0;

    err = incompatible_options(opt, OPT_C | OPT_R | OPT_T);
    if (err) {
	return err;
    }

    tmp = malloc(4 * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    *s = '\0';
    strncat(s, line, MAXLEN - 1);

    for (i=0; i<4 && !err; i++) {
	tmp[i] = strtok((i == 0)? s : NULL, " ");
	if (tmp[i] == NULL && (i < 3)) { 
	    /* 3rd argument is optional */
	    err = E_PARSE;
	}
    }

    if (!err) {
	coef_se = get_matrix_by_name(tmp[1]);
	parnames = get_string_by_name(tmp[2]);
	if (coef_se == NULL || parnames == NULL) {
	    err = E_DATA;
	}
    }

    if (!err && gretl_matrix_cols(coef_se) != 2) {
	gretl_errmsg_set(_("modprint: the first matrix argument must have 2 columns"));
	err = E_DATA;
    }

    if (!err && tmp[3] != NULL && *tmp[3] != '\0') {
	addstats = get_matrix_by_name(tmp[3]);
	if (addstats == NULL) {
	    err = E_DATA;
	}
    }

    if (!err) {
	err = print_user_model(coef_se, addstats, parnames, opt, prn);
    }

    free(tmp);

    return err;
}

