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

static int 
plain_print_coefficients (const MODEL *pmod, const DATASET *dset, PRN *prn);
static int 
alt_print_coefficients (const MODEL *pmod, const DATASET *dset, PRN *prn);
static void alt_print_rho_terms (const MODEL *pmod, PRN *prn);
static void logit_probit_stats (const MODEL *pmod, PRN *prn);
static void print_arma_roots (const MODEL *pmod, PRN *prn);
static void print_heckit_stats (const MODEL *pmod, PRN *prn);

#define RTFTAB "\\par \\ql \\tab "

#define XDIGITS(m) (((m)->ci == MPOLS)? GRETL_MP_DIGITS : GRETL_DIGITS)

#define FDIGITS(m) (((m)->ci == MPOLS)? GRETL_MP_DIGITS : 5)

#define ordered_model(m) ((m->ci == LOGIT || m->ci == PROBIT) && \
                           gretl_model_get_int(m, "ordered"))

#define binary_model(m) ((m->ci == LOGIT || m->ci == PROBIT) && \
                         !gretl_model_get_int(m, "ordered") && \
                         !gretl_model_get_int(m, "multinom"))

#define multinomial_model(m) (m->ci == LOGIT && gretl_model_get_int(m, "multinom"))

#define logit_probit_model(m) (m->ci == LOGIT || m->ci == PROBIT)

#define re_probit_model(m) (m->ci == PROBIT && (m->opt & OPT_E))

#define liml_equation(m) (gretl_model_get_int(m, "method") == SYS_METHOD_LIML)

#define tsls_model(m) (m->ci == IVREG && \
                       !(m->opt & OPT_L) && \
	               !(m->opt & OPT_G) && \
                       !m->aux)

#define liml_model(m) (m->ci == IVREG && (m->opt & OPT_L))

#define gmm_model(m) (m->ci == GMM || (m->ci == IVREG && (m->opt & OPT_G)))

#define intreg_model(m) (m->ci == INTREG || m->ci == TOBIT)

#define hessian_maybe_fishy(m) (m->ci == ARMA || \
				(m->ci == PROBIT && (m->opt & OPT_E)))

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
    if (x < 0 && gretl_print_has_minus(prn)) {
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
	pputs(prn, "\\begin{tabular}{lr@{.}l}\n");
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
	double LR = gretl_model_get_double(pmod, "garch_LR");
	int LRdf = gretl_model_get_int(pmod, "garch_LR_df");
	
	pprintf(prn, "%s = %.*g\n", _(varstr), GRETL_DIGITS, v);
	if (pmod->opt & OPT_Z) {
	    pprintf(prn, "%s\n", _("The residuals are standardized"));
	}
	if (LR >= 0 && LRdf > 0) {
	    pprintf(prn, "%s:\n", _("Likelihood ratio test for (G)ARCH terms"));
	    pprintf(prn, "  %s(%d) = %g [%g]\n",
		    _("Chi-square"), LRdf, LR, chisq_cdf_comp(LRdf, LR));
	}
	pputc(prn, '\n');
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g\n", A_(varstr), v);
    } else if (tex_format(prn)) {
	pprintf(prn, "%s = %g \\\\\n", A_(varstr), v);
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", A_(varstr), prn_delim(prn), v);
    }
}

static void print_intreg_info (const MODEL *pmod, 
			       const DATASET *dset,
			       PRN *prn)
{
    const char *nstrs[] = {
	N_("Left-unbounded observations"),
	N_("Right-unbounded observations"),
	N_("Bounded observations"),
	N_("Point observations"),
	N_("Pseudo-point observations"),
    };
    gchar *lstr = NULL, *rstr = NULL;
    int nl = gretl_model_get_int(pmod, "n_left");
    int nr = gretl_model_get_int(pmod, "n_right");
    int nb = -1, np = -1, nfp = -1;
    double llim = 0, rlim = NADBL;
    double se_sigma;

    if (pmod->ci == INTREG) {
	nb = gretl_model_get_int(pmod, "n_both");
	np = gretl_model_get_int(pmod, "n_point");
	nfp = gretl_model_get_int(pmod, "n_fpoint");
    } else {
	nstrs[0] = N_("Left-censored observations");
	nstrs[1] = N_("Right-censored observations");
	if (pmod->opt & OPT_L) {
	    llim = gretl_model_get_double(pmod, "llimit");
	    if (tex_format(prn)) {
		lstr = g_strdup_printf(" (%s $\\le$ %g)", 
				       dset->varname[pmod->list[1]],
				       llim);
	    } else {
		lstr = g_strdup_printf(" (%s <= %g)", 
				       dset->varname[pmod->list[1]],
				       llim);
	    }
	}
	if ((pmod->opt & OPT_M) && nr > 0) {
	    rlim = gretl_model_get_double(pmod, "rlimit");
	    if (!na(rlim) && tex_format(prn)) {
		rstr = g_strdup_printf(" (%s $\\ge$ %g)", 
				       dset->varname[pmod->list[1]],
				       rlim);
	    } else if (!na(rlim)) {
		rstr = g_strdup_printf(" (%s >= %g)", 
				       dset->varname[pmod->list[1]],
				       rlim);
	    }
	}	
    }

    ensure_vsep(prn);

    se_sigma = gretl_model_get_double(pmod, "se_sigma");

    if (plain_format(prn)) {  
	pprintf(prn, "%s = %.*g", _("sigma"), GRETL_DIGITS, pmod->sigma);
	if (!na(se_sigma)) {
	    pprintf(prn, " (%g)", se_sigma);
	}
	pputc(prn, '\n');
	pprintf(prn, "%s: %d%s\n", _(nstrs[0]), nl, lstr == NULL ? "" : lstr);
	pprintf(prn, "%s: %d%s\n", _(nstrs[1]), nr, rstr == NULL ? "" : rstr);
	if (nb >= 0 && np >= 0) {
	    pprintf(prn, "%s: %d\n", _(nstrs[2]), nb);
	    pprintf(prn, "%s: %d\n", _(nstrs[3]), np);
	}
	if (nfp > 0) {
	    pprintf(prn, "%s: %d\n", _(nstrs[4]), nfp);
	}
	pputc(prn, '\n');
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s = %g", A_("sigma"), pmod->sigma);
	if (!na(se_sigma)) {
	    pprintf(prn, " (%g)", se_sigma);
	}
	pputc(prn, '\n');
	pprintf(prn, RTFTAB "%s: %d%s\n", A_(nstrs[0]), nl, lstr == NULL ? "" : lstr);
	pprintf(prn, RTFTAB "%s: %d%s\n", A_(nstrs[1]), nr, rstr == NULL ? "" : rstr);
	if (nb >= 0 && np >= 0) {
	    pprintf(prn, RTFTAB "%s: %d\n", A_(nstrs[2]), nb);
	    pprintf(prn, RTFTAB "%s: %d\n", A_(nstrs[3]), np);
	}
	if (nfp > 0) {
	    pprintf(prn, RTFTAB "%s: %d\n", A_(nstrs[4]), nfp);
	}
    } else if (tex_format(prn)) {
	pprintf(prn, "$\\hat{\\sigma}$ = %g", pmod->sigma);
	if (!na(se_sigma)) {
	    pprintf(prn, " (%g)", se_sigma);
	}
	pputs(prn, " \\\\\n");
	pprintf(prn, "%s: %d%s \\\\\n", A_(nstrs[0]), nl, lstr == NULL ? "" : lstr);
	pprintf(prn, "%s: %d%s \\\\\n", A_(nstrs[1]), nr, rstr == NULL ? "" : rstr);
	if (nb >= 0 && np >= 0) {
	    pprintf(prn, "%s: %d \\\\\n", A_(nstrs[2]), nb);
	    pprintf(prn, "%s: %d \\\\\n", A_(nstrs[3]), np);
	}
	if (nfp > 0) {
	    pprintf(prn, "%s: %d \\\\\n", A_(nstrs[4]), nfp);
	}
    } else if (csv_format(prn)) {
	int d = prn_delim(prn);

	pprintf(prn, "%s%c%.15g\n", A_("sigma"), d, pmod->sigma);
	pprintf(prn, "\"%s\"%c%d\n", A_(nstrs[0]), d, nl);
	pprintf(prn, "\"%s\"%c%d\n", A_(nstrs[1]), d, nr);
	if (nb >= 0 && np >= 0) {
	    pprintf(prn, "\"%s\"%c%d\n", A_(nstrs[2]), d, nb);
	    pprintf(prn, "\"%s\"%c%d\n", A_(nstrs[3]), d, np);
	}
	if (nfp > 0) {
	    pprintf(prn, "\"%s\"%c%d\n", A_(nstrs[4]), d, nfp);
	}
    }

    g_free(lstr);
    g_free(rstr);
}

static void rsqline (const MODEL *pmod, PRN *prn)
{
    if (!na(pmod->rsq) && plain_format(prn)) {
	pprintf(prn, "  %s = %f\n", _("Unadjusted R-squared"), pmod->rsq);
    }
}

static void ssrline (const MODEL *pmod, PRN *prn)
{
    if (!na(pmod->ess) && plain_format(prn)) { 
	pprintf(prn, "  %s = %.*g\n", _("Sum of squared residuals"), 
		XDIGITS(pmod), pmod->ess);
    }
}

static void rssline (const MODEL *pmod, PRN *prn)
{
    if (!na(pmod->ess) && !na(pmod->tss) && plain_format(prn)) { 
	pprintf(prn, "  %s = %.*g\n", _("Explained sum of squares"), 
		XDIGITS(pmod), pmod->tss - pmod->ess);
    } 
}

static void print_liml_equation_data (const MODEL *pmod, PRN *prn)
{
    double lmin = gretl_model_get_double(pmod, "lmin");
    int idf = gretl_model_get_int(pmod, "idf");

    if (!na(lmin)) {
	ensure_vsep(prn);
	if (!pmod->aux) {
	    if (tex_format(prn)) {
		pprintf(prn, "%s = %g\\\\\n", A_("Smallest eigenvalue"), lmin);
	    } else {
		pprintf(prn, "%s = %g\n", A_("Smallest eigenvalue"), lmin);
	    }
	}
	if (idf > 0) {
	    double X2 = pmod->nobs * log(lmin);
	    double pv = chisq_cdf_comp(idf, X2);

	    if (tex_format(prn)) {
		pprintf(prn, "%s: $\\chi^2(%d)$ = %g [%.4f] \\\\\n", 
			A_("LR over-identification test"), idf, X2, pv);
	    } else if (rtf_format(prn)) {
		pprintf(prn, "%s: ", A_("LR over-identification test"));
		pprintf(prn, "%s(%d) = %g [%.4f]\n\n", A_("Chi-square"),
			idf, X2, pv);
	    } else {
		pprintf(prn, "%s: ", A_("LR over-identification test"));
		pprintf(prn, "%s(%d) = %g [%.4f]\n\n", A_("Chi-square"),
			idf, X2, pv);
	    }
	} else if (idf == 0) {
	    pprintf(prn, "%s\n\n", A_("Equation is just identified"));
	}
    }
}

static void print_panel_AR_test (double z, int order, PRN *prn)
{
    pprintf(prn, A_("Test for AR(%d) errors:"), order);

    if (na(z)) {
	if (tex_format(prn)) {
	    pputs(prn, " & $z$ = NA");
	} else {
	    pputs(prn, " z = NA");
	}
    } else {
	double pv = normal_pvalue_2(z);

	if (tex_format(prn)) {
	    char numstr[32];

	    tex_sprint_double_digits(z, numstr, 4);
	    pprintf(prn, " & $z$ = %s [%.4f]", numstr, pv);
	} else {
	    pprintf(prn, " z = %g [%.4f]", z, pv);
	}
    }

    gretl_prn_newline(prn);
}

enum {
    AB_SARGAN,
    AB_WALD,
    AB_WALD_TIME,
    J_TEST,
    OVERDISP
};

static void 
print_model_chi2_test (const MODEL *pmod, double x, int j, PRN *prn)
{
    const char *strs[] = {
	N_("Sargan over-identification test"),
	N_("Wald (joint) test"),
	N_("Wald (time dummies)"),
	N_("J test"),
	N_("Overdispersion test")
    };
    const char *texstrs[] = {
	N_("Sargan test"),
	N_("Wald (joint) test"),
	N_("Wald (time dummies)"),
	N_("J test"),
	N_("Overdispersion test")
    };
    double pv;
    int df;

    if (j == AB_SARGAN) {
	df = gretl_model_get_int(pmod, "sargan_df");
    } else if (j == AB_WALD) {
	df = gretl_model_get_int(pmod, "wald_df");
    } else if (j == AB_WALD_TIME) {
	df = gretl_model_get_int(pmod, "wald_time_df");
    } else if (j == J_TEST) {
	df = gretl_model_get_int(pmod, "J_df");
    } else if (j == OVERDISP) {
	df = 1;
    } else {
	return;
    }

    if (na(x)) {
	if (j == AB_SARGAN || j == AB_WALD) {
	    if (df < 0) {
		df = 0;
	    }
	    if (tex_format(prn)) {
		pprintf(prn, "%s: ", A_(texstrs[j]));
		pprintf(prn, " $\\chi^2(%d)$ = NA", df);
	    } else if (plain_format(prn)) {
		if (gmm_model(pmod)) {
		    pputs(prn, "  ");
		}
		pprintf(prn, "%s: ", _(strs[j]));
		pprintf(prn, "%s(%d) = NA", _("Chi-square"), df);
	    } else {
		pprintf(prn, "%s: ", A_(strs[j]));
		pprintf(prn, "%s(%d) = NA", A_("Chi-square"), df);
	    }
	    gretl_prn_newline(prn);
	}
	return;
    }

    pv = chisq_cdf_comp(df, x);

    if (na(pv)) {
	return;
    }

    if (tex_format(prn)) {
	pprintf(prn, "%s: ", A_(texstrs[j]));
	pprintf(prn, " $\\chi^2(%d)$ = %g [%.4f]", df, x, pv);
    } else if (plain_format(prn)) {
	if (gmm_model(pmod)) {
	    pputs(prn, "  ");
	}
	pprintf(prn, "%s: ", _(strs[j]));
	pprintf(prn, "%s(%d) = %g [%.4f]", _("Chi-square"), df, x, pv);
    } else {
	pprintf(prn, "%s: ", A_(strs[j]));
	pprintf(prn, "%s(%d) = %g [%.4f]", A_("Chi-square"), df, x, pv);
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
	pprintf(prn, RTFTAB "%s: Q = %g (TQ = %g)\n", A_("GMM criterion"), 
		Q, TQ);
    } else if (tex_format(prn)) {
	char x1[32], x2[32];

	tex_sprint_double(Q, x1);
	tex_sprint_double(TQ, x2);
	pprintf(prn, "%s, $Q$ = %s ($TQ$ = %s)\\\\\n",
		A_("GMM criterion"), x1, x2);
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%.15g\n", A_("GMM criterion"), 
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
	print_model_chi2_test(pmod, x, J_TEST, prn);
    }
    
    if (!tex_format(prn)) {
	gretl_prn_newline(prn);
    }
}

static void print_DPD_stats (const MODEL *pmod, PRN *prn)
{
    double x;
    int k;

    ensure_vsep(prn);

    if (tex_format(prn)) {
	pputs(prn, "\\begin{tabular}{ll}\n");
    }

    k = gretl_model_get_int(pmod, "ninst");
    if (k > 0) {
	pprintf(prn, A_("Number of instruments = %d"), k);
	gretl_prn_newline(prn);
    }

    x = gretl_model_get_double(pmod, "AR1");
    print_panel_AR_test(x, 1, prn);

    x = gretl_model_get_double(pmod, "AR2");
    print_panel_AR_test(x, 2, prn);

    x = gretl_model_get_double(pmod, "sargan");
    print_model_chi2_test(pmod, x, AB_SARGAN, prn);

    x = gretl_model_get_double(pmod, "wald");
    print_model_chi2_test(pmod, x, AB_WALD, prn);

    x = gretl_model_get_double(pmod, "wald_time");
    if (!na(x)) {
	print_model_chi2_test(pmod, x, AB_WALD_TIME, prn);
    }

    if (tex_format(prn)) {
	pputs(prn, "\\end{tabular}\n");
    } else {
	gretl_prn_newline(prn);
    }
}

static void print_overdisp_test (const MODEL *pmod, PRN *prn)
{
    double x = gretl_model_get_double(pmod, "overdisp");

    if (!na(x)) {
	ensure_vsep(prn);
	print_model_chi2_test(pmod, x, OVERDISP, prn);
	if (!tex_format(prn)) {
	    gretl_prn_newline(prn);
	}
    }
}

static void print_duration_alpha (const MODEL *pmod, PRN *prn)
{
    if ((pmod->opt & OPT_W) && plain_format(prn)) {
	double a = 1.0 / pmod->coeff[pmod->ncoeff - 1];
	double sa = a * a * pmod->sderr[pmod->ncoeff - 1];

	ensure_vsep(prn);
	pprintf(prn, "1/sigma = %g (%g)\n", a, sa);
	if (!tex_format(prn)) {
	    gretl_prn_newline(prn);
	}	
    }
}

static void print_probit_rho (const MODEL *pmod, PRN *prn)
{
    ensure_vsep(prn);

    if (re_probit_model(pmod)) {
	if (tex_format(prn)) {
	    pprintf(prn, "$\\hat{\\sigma}_u$ = %.5f\n", pmod->sigma);
	} else {
	    pprintf(prn, "sigma_u = %g\n", pmod->sigma);
	}
    }

    if (tex_format(prn)) {
	double r = pmod->rho;

	if (r < 0) {
	    pprintf(prn, "$\\hat{\\rho}$ = $-$%.5f\n", fabs(r));
	} else {
	    pprintf(prn, "$\\hat{\\rho}$ = %.5f\n", r);
	}	    
    } else {
	pprintf(prn, "rho = %g\n", pmod->rho);
    }

    gretl_prn_newline(prn);
}

static void maybe_print_lad_warning (const MODEL *pmod, PRN *prn)
{
    if (gretl_model_get_int(pmod, "nonunique")) {
	pputs(prn, A_("Warning: solution is probably not unique"));
	pputc(prn, '\n');
    }
}

static void maybe_print_hessian_warning (const MODEL *pmod, PRN *prn)
{
    if (gretl_model_get_int(pmod, "hess-error")) {
	pputs(prn, A_("Warning: couldn't compute numerical Hessian"));
	pputc(prn, '\n');
    }
    if (gretl_model_get_int(pmod, "non-pd-hess")) {
	pputs(prn, A_("Warning: non-pd Hessian (but still nonsingular)"));
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
	pprintf(prn, RTFTAB "%s = %g", A_("'Within' variance"), ws2);
	pprintf(prn, RTFTAB "%s = %g", A_("'Between' variance"), bs2);
	if (!na(theta)) {
	    pprintf(prn, RTFTAB "%s = %g", A_("theta used for quasi-demeaning"), theta);
	}
    } else if (csv_format(prn)) {
	char d = prn_delim(prn);

	pprintf(prn, "\"%s\"%c%.15g\n", A_("'Within' variance"), d, ws2);
	pprintf(prn, "\"%s\"%c%.15g\n", A_("'Between' variance"), d, bs2);
	if (!na(theta)) {
	    pprintf(prn, "\"%s\"%c%.15g\n", A_("theta used for quasi-demeaning"), 
		    d, theta);
	}
    }	
}

static double durbins_h (const MODEL *pmod)
{
    int ldv = gretl_model_get_int(pmod, "ldepvar");
    double se = pmod->sderr[ldv - 2];
    int T = pmod->nobs;
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

static void pval_max_line (const MODEL *pmod, const DATASET *dset, 
			   PRN *prn)
{
    int k = pmod->ncoeff - pmod->ifc;

    if (k < 3) return;

    if ((k = least_significant_coeff(pmod))) {
	char tmp[128];

	if (pmod->ifc) {
	    sprintf(tmp, _("Excluding the constant, p-value was highest "
			   "for variable %d (%s)"), k, dset->varname[k]);
	} else {
	    sprintf(tmp, _("P-value was highest for variable %d (%s)"), 
		    k, dset->varname[k]);
	}	    
	pprintf(prn, "%s\n\n", tmp);
    }
}

static void print_aux_string (const MODEL *pmod, PRN *prn)
{
    int aux = pmod->aux;
    int plain = plain_format(prn);
    int tex = tex_format(prn);
    int csv = csv_format(prn);
    int close = 1;

    if (plain || tex) {
	pputc(prn, '\n');
    } else if (csv) {
	pputc(prn, '"');
    } 

    if (aux == AUX_SQ) {
	pputs(prn, A_("Auxiliary regression for non-linearity test "
		      "(squared terms)"));
    } else if (aux == AUX_LOG) {
	pputs(prn, A_("Auxiliary regression for non-linearity test "
		      "(log terms)"));
    } else if (aux == AUX_ADD) {
	pputs(prn, A_("Auxiliary regression for added variables"));
    } else if (aux == AUX_WHITE) {
	pputs(prn, A_("White's test for heteroskedasticity"));
	if (pmod->opt & OPT_X) {
	    pprintf(prn, " (%s)", A_("squares only"));
	} 
    } else if (aux == AUX_BP) {
	pputs(prn, A_("Breusch-Pagan test for heteroskedasticity"));
    } else if (aux == AUX_HET_1) {
	pputs(prn, A_("Pesaran-Taylor test for heteroskedasticity"));
    } else if (aux == AUX_CHOW) {
	pputs(prn, A_("Augmented regression for Chow test"));
    } else if (aux == AUX_COINT) {
	if (tex_format(prn)) {
	    pputs(prn, A_("Cointegrating regression -- "));
	} else {
	    pputs(prn, A_("Cointegrating regression - "));
	}
    } else if (aux == AUX_ADF) {
	if (tex_format(prn)) {
	    pputs(prn, A_("Augmented Dickey--Fuller regression"));
	} else {
	    pputs(prn, A_("Augmented Dickey-Fuller regression"));
	}
    } else if (aux == AUX_DF) {
	if (tex_format(prn)) {
	    pputs(prn, A_("Dickey--Fuller regression"));
	} else {
	    pputs(prn, A_("Dickey-Fuller regression"));
	}
    } else if (aux == AUX_KPSS) {
	pputs(prn, A_("KPSS regression"));
    } else if (aux == AUX_RESET) {
	pputs(prn, A_("Auxiliary regression for RESET specification test"));
    } else if (aux == AUX_GROUPWISE) {
	pputs(prn, A_("Groupwise heteroskedasticity"));
    } else if (aux == AUX_COMFAC) {
	pputs(prn, A_("Augmented regression for common factor test"));
    } else {
	close = 0;
    }

    if (close) {
	if (plain || tex) {
	    pputc(prn, '\n');
	} else if (csv) {
	    pputs(prn, "\"\n");
	} else { /* RTF */
	    pputs(prn, "\\par\n");
	}
    }	      
}

static const char *simple_estimator_string (int ci, PRN *prn)
{
    if (ci == OLS || ci == VAR) return N_("OLS");
    else if (ci == WLS)  return N_("WLS"); 
    else if (ci == ARCH) return N_("WLS (ARCH)");
    else if (ci == HSK)  return N_("Heteroskedasticity-corrected");
    else if (ci == AR)   return N_("AR");
    else if (ci == LAD)  return N_("LAD");
    else if (ci == MPOLS) return N_("High-Precision OLS");
    else if (ci == PROBIT) return N_("Probit");
    else if (ci == LOGIT)  return N_("Logit");
    else if (ci == TOBIT)  return N_("Tobit");
    else if (ci == HECKIT) return N_("Heckit");
    else if (ci == POISSON) return N_("Poisson");
    else if (ci == NEGBIN) return N_("Negative Binomial");
    else if (ci == DURATION) return N_("Duration");
    else if (ci == NLS) return N_("NLS");
    else if (ci == MLE) return N_("ML");
    else if (ci == GMM) return N_("GMM");
    else if (ci == LOGISTIC) return N_("Logistic");
    else if (ci == GARCH) return N_("GARCH");
    else if (ci == INTREG) return N_("Interval estimates");
    else if (ci == DPANEL) return N_("Dynamic panel");
    else if (ci == BIPROBIT) return N_("Bivariate probit");
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
	if (pmod->opt & OPT_H) {
	    if (tex_format(prn)) return N_("Hildreth--Lu");
	    else return N_("Hildreth-Lu");
	} else if (pmod->opt & OPT_P) {
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
	if (pmod->opt & OPT_F) {
	    return N_("Fixed-effects");
	} else if (pmod->opt & OPT_U) {
	    return N_("Random-effects (GLS)");
	} else if (pmod->opt & OPT_W) {
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
    } else if (pmod->ci == DPANEL) {
	if (gretl_model_get_int(pmod, "step") == 2) {
	    return N_("2-step dynamic panel");
	} else {
	    return N_("1-step dynamic panel");
	}
    } else if (gmm_model(pmod)) {
	if (pmod->opt & OPT_T) {
	    return N_("2-step GMM");
	} else if (pmod->opt & OPT_I) {
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
	} else if (gretl_model_get_int(pmod, "multinom")) {
	    return N_("Multinomial Logit");
	} else {
	    return N_("Logit");
	}
    } else if (pmod->ci == PROBIT) {
	if (pmod->opt & OPT_E) {
	    return N_("Random-effects probit");
	} else if (gretl_model_get_int(pmod, "ordered")) {
	    return N_("Ordered Probit");
	} else {
	    return N_("Probit");
	}
    } else if (pmod->ci == HECKIT) {
	if (pmod->opt & OPT_T) {
	    return N_("Two-step Heckit");
	} else {
	    return N_("ML Heckit");
	}
    } else if (pmod->ci == LAD) {
	if (gretl_model_get_int(pmod, "rq")) {
	    return N_("Quantile estimates");
	} else {
	    return N_("LAD");
	}
    } else if (pmod->ci == IVREG) {
	if (pmod->opt & OPT_L) {
	    return N_("LIML");
	} else {
	    return N_("TSLS");
	}
    } else if (pmod->ci == NEGBIN) {
	if (pmod->opt & OPT_M) {
	    return N_("Negative Binomial 1");
	} else {
	    return N_("Negative Binomial");
	}
    } else if (pmod->ci == DURATION) {
	if (pmod->opt & OPT_E) {
	    return N_("Duration (exponential)");
	} else if (pmod->opt & OPT_L) {
	    return N_("Duration (log-logistic)");
	} else if (pmod->opt & OPT_Z) {
	    return N_("Duration (log-normal)");
	} else {
	    return N_("Duration (Weibull)");
	}
    } else if (pmod->ci == OLS && 
	       gretl_model_get_int(pmod, "restricted")) {
	return N_("Restricted OLS");
    } else {
	return simple_estimator_string(pmod->ci, prn);
    }
}

static int any_tests (const MODEL *pmod)
{
    if (pmod->ntests > 0) {
	return 1;
    }

    if (pmod->ci == IVREG && gretl_model_get_int(pmod, "stage1-dfn")) {
	return 1;
    }

    return 0;
}

static void get_stock_yogo_critvals (int n, int K2, gretlopt opt,
				     gretl_matrix **p1,
				     gretl_matrix **p2)
{
    gretl_matrix *(*lookup) (int, int, int);
    void *handle;

    lookup = get_plugin_function("stock_yogo_lookup", &handle);

    if (lookup != NULL) {
	if (opt & OPT_L) {
	    /* LIML test size */
	    *p2 = (*lookup) (n, K2, 3);
	} else {
	    /* TSLS relative bias, test size */
	    *p1 = (*lookup) (n, K2, 1);
	    *p2 = (*lookup) (n, K2, 2);
	}
	close_plugin(handle);
    }
}

static void plain_print_sy_vals (gretl_matrix *v, double g, int k, 
				 gretlopt opt, PRN *prn)
{
    int i, gpos = -1;
    double x;

    pputs(prn, "  ");

    if (k == 1) {
	pputs(prn, _("Critical values for TSLS bias relative to OLS:\n"));
    } else if (opt & OPT_L) {
	/* xgettext:no-c-format */
	pputs(prn, _("Critical values for desired LIML maximal size, when running\n"
		     "  tests at a nominal 5% significance level:\n"));
    } else {
	/* xgettext:no-c-format */	
	pputs(prn, _("Critical values for desired TSLS maximal size, when running\n"
		     "  tests at a nominal 5% significance level:\n"));
    }

    pprintf(prn, "\n%9s", (k == 1)? _("bias") : _("size"));
    for (i=0; i<4; i++) {
	pprintf(prn, "%8g%%", 100 * gretl_matrix_get(v, 0, i));
    }
    pprintf(prn, "\n%9s", _("value"));
    for (i=0; i<4; i++) {
	x = gretl_matrix_get(v, 1, i);
	if (gpos < 0 && g > x) {
	    gpos = i;
	}	
	pprintf(prn, "%9.2f", x);
    }    
    pputs(prn, "\n\n  ");

    if (gpos == 0) {
	x = gretl_matrix_get(v, 0, 0);
	if (k == 1) {
	    pprintf(prn, _("Relative bias is probably less than %g%%"), 100 * x);
	} else {
	    pprintf(prn, _("Maximal size is probably less than %g%%"), 100 * x);
	}
    } else {
	if (gpos < 0) {
	    gpos = 3;
	} else {
	    gpos--;
	}
	x = gretl_matrix_get(v, 0, gpos);
	if (k == 1) {
	    pprintf(prn, _("Relative bias may exceed %g%%"), 100 * x);
	} else {	
	    pprintf(prn, _("Maximal size may exceed %g%%"), 100 * x);
	}
    }
    
    pputs(prn, "\n\n");
}

static void maybe_print_weak_insts_test (const MODEL *pmod, PRN *prn)
{
    const char *head = N_("Weak instrument test");
    double F = gretl_model_get_double(pmod, "stage1-F");
    double g = gretl_model_get_double(pmod, "gmin");
    int got_critvals = 0;
    int dfn = 0, dfd = 0;

    if (na(F) && na(g)) {
	return;
    }

    if (plain_format(prn)) {
	pprintf(prn, "%s - \n", _(head));
    } else if (tex_format(prn)) {
	pprintf(prn, "%s -- \\\\\n", A_(head));
    } else if (rtf_format(prn)) {
	pprintf(prn, "%s - \\par\n", A_(head));
    }

    if (!na(F)) {
	/* got first-stage F-test (single endogenous regressor) */
	dfn = gretl_model_get_int(pmod, "stage1-dfn");
	dfd = gretl_model_get_int(pmod, "stage1-dfd");

	if (plain_format(prn)) {
	    pprintf(prn, "  %s (%d, %d) = %.*g\n", _("First-stage F-statistic"),
		    dfn, dfd, GRETL_DIGITS, F);
	} else if (tex_format(prn)) {
	    char x1str[32];

	    tex_sprint_double(F, x1str);
	    pprintf(prn, "\\quad First-stage $F(%d, %d)$ = %s \\\\\n", dfn, dfd, x1str);
	} else if (rtf_format(prn)) {
	    pprintf(prn, "  %s (%d, %d) = %g\n", A_("First-stage F-statistic"), 
		    dfn, dfd, F);
	}
    } else {
	/* got minimum eigenvalue test statistic */
	if (plain_format(prn)) {
	    pprintf(prn, "  %s = %.*g\n", _("Cragg-Donald minimum eigenvalue"),
		    GRETL_DIGITS, g);
	} else if (tex_format(prn)) {
	    char x1str[32];

	    tex_sprint_double(g, x1str);
	    pprintf(prn, "\\quad  %s = %s \\\\\n", A_("Cragg--Donald minimum eigenvalue"), 
		    x1str);
	} else if (rtf_format(prn)) {
	    pprintf(prn, "  %s = %g\n", A_("Cragg-Donald minimum eigenvalue"), g);
	}
    }	

    if (!na(g)) {
	/* print Stock-Yogo critical values, if available */
	gretl_matrix *bvals = NULL;
	gretl_matrix *svals = NULL;
	int n, K2;

	if (na(F)) {
	    n = gretl_model_get_int(pmod, "n");
	    K2 = gretl_model_get_int(pmod, "K2");
	} else {
	    n = 1;
	    K2 = dfn;
	}

	get_stock_yogo_critvals(n, K2, pmod->opt, &bvals, &svals);

	if (bvals != NULL) {
	    got_critvals = 1;
	    if (plain_format(prn)) {
		plain_print_sy_vals(bvals, g, 1, pmod->opt, prn);
	    }
	    gretl_matrix_free(bvals);
	}

	if (svals != NULL) {
	    got_critvals = 1;
	    if (plain_format(prn)) {
		plain_print_sy_vals(svals, g, 2, pmod->opt, prn);
	    }
	    gretl_matrix_free(svals);
	}	
    }

    if (!na(F) && !got_critvals && plain_format(prn)) {
	pprintf(prn, "  %s\n\n", _("A value < 10 may indicate weak instruments"));
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
	maybe_print_weak_insts_test(pmod, prn);
	pputs(prn, "\\end{raggedright}\n");
    } else {
	for (i=0; i<pmod->ntests; i++) {
	    gretl_model_test_print(pmod, i, prn);
	}
	maybe_print_weak_insts_test(pmod, prn);
    }
}

static int 
print_ivreg_instruments (const MODEL *pmod, const DATASET *dset, PRN *prn)
{
    const char *strs[] = {
	N_("Instrumented"),
	N_("Instruments")
    };
    const int *list;
    char vname[16];
    int tex = tex_format(prn);
    int i, j, imin, jmin, vi;

    jmin = (pmod->aux)? 1 : 0;

    for (j=jmin; j<2; j++) {
	int ccount, nv = 0;

	if (j == 0) {
	    list = (const int *) gretl_model_get_data(pmod, "endolist");
	    imin = list == NULL ? -1 : 1;
	} else {
	    list = pmod->list;
	    imin = gretl_list_separator_position(list);
	    if (imin > 0) {
		imin++;
	    } else {
		imin = -1;
	    }
	}

	if (imin < 0) {
	    continue;
	}

	ccount = pprintf(prn, "%s: ", A_(strs[j]));
    
	for (i=imin; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi >= dset->v) {
		fprintf(stderr, "print_ivreg_instruments: bad varnum %d\n", vi);
		continue;
	    }
	    if (tex) {
		tex_escape(vname, dset->varname[vi]);
	    } else {
		strcpy(vname, dset->varname[vi]);
	    }
	    ccount += pprintf(prn, "%s ", vname);
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
	    nv++;
	}
    
	if (nv == 0) {
	    pputs(prn, "none??");
	}

	gretl_prn_newline(prn);
    }

    return 0;
}

static void dpd_asy_vcv_line (PRN *prn)
{
    if (csv_format(prn)) {
	pprintf(prn, "\"%s\"", A_("Asymptotic standard errors"));
    } else {
	pputs(prn, A_("Asymptotic standard errors"));
    }

    pputc(prn, '\n');
}

static void panel_vcv_line (const VCVInfo *vi, PRN *prn)
{
    if (vi->vmin == PANEL_HAC) {
	if (csv_format(prn)) {
	    pprintf(prn, "\"%s\"", A_("Robust (HAC) standard errors"));
	} else {
	    pputs(prn, A_("Robust (HAC) standard errors"));
	}
	pputc(prn, '\n');
    } else if (vi->vmin == PANEL_BK) {
	if (csv_format(prn)) {
	    pprintf(prn, "\"%s\"", A_("Beck-Katz standard errors"));
	} else if (tex_format(prn)) {
	    pputs(prn, A_("Beck--Katz standard errors"));
	} else {
	    pputs(prn, A_("Beck-Katz standard errors"));
	} 	
	pputc(prn, '\n');
    }
}

static void cluster_vcv_line (const MODEL *pmod, const VCVInfo *vi, 
			      const DATASET *dset, PRN *prn)
{
    gchar *cstr;

    if (vi->vmin >= 1 && vi->vmin < dset->v) {
	int n_c = gretl_model_get_int(pmod, "n_clusters");

	cstr = g_strdup_printf(A_("Standard errors clustered by %d values of %s"),
			       n_c, dset->varname[vi->vmin]);
    } else {
	cstr = g_strdup(A_("Clustered standard errors"));
    }

    if (csv_format(prn)) {
	pprintf(prn, "\"%s\"", cstr);
    } else {
	pputs(prn, cstr);
    } 
    pputc(prn, '\n');

    g_free(cstr);
}

static void beck_katz_failed_line (PRN *prn)
{
    if (plain_format(prn)) {
	pputs(prn, A_("Could not compute Beck-Katz standard errors"));
	pputc(prn, '\n');
    }
}

static void hac_vcv_line (const VCVInfo *vi, PRN *prn)
{
    const char *kstrs[] = {
	N_("Bartlett kernel"),
	N_("Parzen kernel"),
	N_("QS kernel")
    };

    if (vi->vmin == KERNEL_QS) {
	pprintf(prn, A_("HAC standard errors, "
			"bandwidth %.2f"), vi->bw);
    } else {
	pprintf(prn, A_("HAC standard errors, "
			"bandwidth %d"), vi->order);
    }

    pprintf(prn, " (%s", A_(kstrs[vi->vmin]));
		    
    if (vi->flags) {
	pprintf(prn, ", %s", A_("prewhitened"));
    }

    pputs(prn, ")\n");
}

static void hc_vcv_line (const VCVInfo *vi, PRN *prn)
{
    int hcv = vi->vmin;
    int jack = 0;

    if (vi->vmin == 4) {
	jack = 1;
	hcv--;
    }

    pprintf(prn, "%s, %s%sHC%d%s", 
	    A_("Heteroskedasticity-robust standard errors"),
	    (jack)? "" : A_("variant"),
	    (jack)? "" : " ",
	    hcv, (jack)? " (jackknife)" : "");

    if (rtf_format(prn)) {
	pputs(prn, "\\par\n");
    } else {
	pputc(prn, '\n');
    }
}

static void ml_vcv_line (const VCVInfo *vi, PRN *prn)
{
    int tex = tex_format(prn);
    const char *s = NULL;

    switch (vi->vmin) {
    case ML_HESSIAN:
	s = N_("Standard errors based on Hessian");
	break;
    case ML_IM:
	s = N_("Standard errors based on Information Matrix");
	break;
    case ML_OP:
	s = N_("Standard errors based on Outer Products matrix");
	break;
    case ML_QML:
	s = N_("QML standard errors");
	break;
    case ML_BW:
	if (tex) {
	    s = N_("Bollerslev--Wooldridge standard errors");
	} else {
	    s = N_("Bollerslev-Wooldridge standard errors");
	}
	break;
    case ML_VCVMAX:
	s = N_("Warning: could not compute standard errors");
	break;
    default:
	break;
    }

    if (s != NULL) {
	if (csv_format(prn)) {
	    pprintf(prn, "\"%s\"\n", A_(s));
	} else {
	    pprintf(prn, "%s\n", A_(s));
	}
    }
}

static void rq_vcv_line (const MODEL *pmod, PRN *prn)
{
    int robust = gretl_model_get_int(pmod, "rq_nid");
    double a = gretl_model_get_double(pmod, "rq_alpha");
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
	pprintf(prn, "\"%s\"", A_(s));
    } else {
	pprintf(prn, "%s", A_(s));
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

static void tex_dpd_depvar_name (char *s, const char *vname)
{
    char vnesc[32];

    tex_escape(vnesc, vname);
    sprintf(s, "$\\Delta$%s", vnesc);
}

void print_model_vcv_info (const MODEL *pmod, const DATASET *dset,
			   PRN *prn)
{
    VCVInfo *vi = NULL;

    if (pmod->ci == LAD && gretl_model_get_int(pmod, "rq")) {
	rq_vcv_line(pmod, prn);
    } else if (gretl_model_get_int(pmod, "panel_bk_failed")) {
	beck_katz_failed_line(prn);
    } else if ((pmod->ci == ARBOND || pmod->ci == DPANEL) && 
	       gretl_model_get_int(pmod, "asy")) {
	dpd_asy_vcv_line(prn);
    } else {
	vi = gretl_model_get_data(pmod, "vcv_info");
    }

    if (vi != NULL) {
	switch (vi->vmaj) {
	case VCV_HC:
	    hc_vcv_line(vi, prn);
	    break;
	case VCV_HAC:
	    hac_vcv_line(vi, prn);
	    break;
	case VCV_ML:
	    ml_vcv_line(vi, prn);
	    break;
	case VCV_PANEL:
	    panel_vcv_line(vi, prn);
	    break;
	case VCV_CLUSTER:
	    cluster_vcv_line(pmod, vi, dset, prn);
	    break;
	default:
	    break;
	}
    }
}

static void print_extra_list (const char *tag, const int *list, 
			      const DATASET *dset, PRN *prn)
{
    int i, v, len;

    len = pputs(prn, A_(tag));

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v < dset->v) {
	    len += pprintf(prn, " %s", dset->varname[v]);
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
				  const DATASET *dset,
				  PRN *prn)
{
    const int *zlist = gretl_model_get_data(pmod, "zerolist");
    const char *tag = N_("Omitted because all values were zero:");

    if (pmod->ci == PANEL && (pmod->opt & OPT_B)) {
	return;
    }

    print_extra_list(tag, zlist, dset, prn);
}

static void print_model_droplist (const MODEL *pmod, 
				  const DATASET *dset,
				  PRN *prn)
{
    const int *dlist = gretl_model_get_data(pmod, "droplist");
    const char *tag = N_("Omitted due to exact collinearity:");

    print_extra_list(tag, dlist, dset, prn);
}

static void print_ivreg_droplist (const MODEL *pmod, 
				  const DATASET *dset,
				  PRN *prn)
{
    const int *dlist = gretl_model_get_data(pmod, "inst_droplist");
    int i, v;

    pputs(prn, A_("Redundant instruments:"));
    for (i=1; i<=dlist[0]; i++) {
	v = dlist[i];
	if (v < dset->v) {
	    pprintf(prn, " %s", dset->varname[v]);
	} else {
	    pprintf(prn, " %d", v);
	}
    }
    pputc(prn, '\n');
}

static void print_arma_depvar (const MODEL *pmod,
			       const DATASET *dset,
			       PRN *prn)
{
    int tex = tex_format(prn);
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
	tex_escape(tmp, dset->varname[yno]);
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
	strcat(vname, dset->varname[yno]);
    }

    pprintf(prn, "%s: %s", A_("Dependent variable"), vname);
}

static void arma_extra_info (const MODEL *pmod, PRN *prn)
{
    int acode = gretl_model_get_int(pmod, "arma_flags");

    if (acode & ARMA_X12A) {
	if (gretl_x12_is_x13()) {
	    pputs(prn, _("Estimated using X-13-ARIMA"));
	} else {
	    pputs(prn, _("Estimated using X-12-ARIMA"));
	}
	pputs(prn, " (");
	pputs(prn, (acode & ARMA_EXACT)? _("exact ML") : _("conditional ML"));
	pputs(prn, ")\n");
    } else if (acode & ARMA_OLS) {
	pputs(prn, _("Estimated using least squares"));
	pputs(prn, " (");
	pputs(prn, _("= MLE"));
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

static void godfrey_test_string (int ci, int order, PRN *prn)
{
    pputc(prn, '\n');

    if (ci == IVREG) {
	if (order > 1) {
	    pprintf(prn, A_("Godfrey (1994) test for autocorrelation up to order %d"), 
		    order);
	} else {
	    pputs(prn, A_("Godfrey (1994) test for first-order autocorrelation"));
	}
    } else {
	if (order > 1) {
	    pprintf(prn, A_("Breusch-Godfrey test for autocorrelation up to order %d"), 
		    order);
	} else {
	    pputs(prn, A_("Breusch-Godfrey test for first-order autocorrelation"));
	}
    }

    pputc(prn, '\n');
}

/* The selection variable should be the first variable following
   the list separator */

static const char *heckit_selvar_name (const MODEL *pmod,
				       const DATASET *dset)
{
    const int *list = pmod->list;
    int pos = gretl_list_separator_position(list);

    if (pos > 0 && pos < list[0] && list[pos+1] < dset->v) {
	return dset->varname[list[pos+1]];
    } else {
	return NULL;
    }
}

static void print_intreg_depvar (const MODEL *pmod,
				 const DATASET *dset,
				 PRN *prn)
{
    int lov = gretl_model_get_int(pmod, "lovar");
    int hiv = gretl_model_get_int(pmod, "hivar");

    if (lov < dset->v && hiv < dset->v) {
	pprintf(prn, "%s: %s", A_("Lower limit"), dset->varname[lov]);
	pprintf(prn, ", %s: %s", A_("Upper limit"), dset->varname[hiv]);
    }
}

static void maybe_print_T (const MODEL *pmod, 
			   const DATASET *dset,
			   const char *start,
			   PRN *prn)
{
    if (pmod->ci == HECKIT) {
	return;
    } else if (!strcmp(start, "1") && !model_has_missing_obs(pmod)) {
	return;
    } else {
	int xsect = dataset_is_cross_section(dset);

	if (model_has_missing_obs(pmod) || !xsect || strcmp(start, "1")) {
	    const char *nstrs[] = {
		/* TRANSLATORS: 'n' denotes sample size */
		N_("n"),
		/* TRANSLATORS: 'T' denotes time-series sample size */
		N_("T")
	    };
	    const char *nstr = xsect ? nstrs[0] : nstrs[1];

	    if (tex_format(prn)) {
		pprintf(prn, " ($%s$ = %d)", A_(nstr), pmod->nobs);
	    } else {
		pprintf(prn, " (%s = %d)", A_(nstr), pmod->nobs);
	    }
	} 
    }
}

static void make_obs_sep (char *targ, const char *obs, int tex)
{
    if (tex) {
	strcpy(targ, "--");
    } else if (strchr(obs, '-') != NULL) {
	strcpy(targ, ":");
    } else {
	strcpy(targ, "-");
    }
}

static void print_model_heading (const MODEL *pmod, 
				 const DATASET *dset, 
				 gretlopt opt, 
				 PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN], vname[32];
    char datesep[4];
    int t1 = pmod->t1, t2 = pmod->t2;
    int tex = tex_format(prn);
    int csv = csv_format(prn);
    int dvnl = 1;
    int order = 0;

    if (pmod->aux != AUX_VAR && pmod->aux != AUX_VECM) {
	ntodate(startdate, t1, dset);
	ntodate(enddate, t2, dset);
	make_obs_sep(datesep, startdate, tex); 
    }

    switch (pmod->aux) {
    case AUX_SQ:
    case AUX_LOG:
    case AUX_ADD:
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
    case AUX_COMFAC:
    case AUX_BIPROB:
	print_aux_string(pmod, prn);
	break;
    case AUX_AR:
	order = gretl_model_get_int(pmod, "BG_order");
	godfrey_test_string(pmod->ci, order, prn);
	break;	
    case AUX_ARCH:
	order = gretl_model_get_int(pmod, "arch_order");
	pputc(prn, '\n');
	pprintf(prn, A_("Test for ARCH of order %d"), order);
	pputc(prn, '\n');
	break;	
    case AUX_SYS:
	pprintf(prn, "%s %d: ", A_("Equation"), pmod->ID + 1);
	break;	
    case AUX_VAR:
	pprintf(prn, "\n%s %d: ", A_("Equation"), pmod->ID);
	break;
    case AUX_VECM:
	pprintf(prn, "%s %d: ", A_("Equation"), pmod->ID);
	break;
    case AUX_AUX:
	pputc(prn, '\n');
	break;
    default:
	if (pmod->ID < 0 || (opt & OPT_S)) {
	    if (!csv) {
		pputc(prn, '\n');
	    }
	} else if (pmod->name) {
	    if (csv) {
		pprintf(prn, "\"%s:\"\n", pmod->name);
	    } else if (strlen(pmod->name) > 8) {
		pprintf(prn, "\n%s:\n", pmod->name);
	    } else {
		pprintf(prn, "\n%s: ", pmod->name);
	    }
	} else {
	    if (csv) {
		pprintf(prn, "\"%s %d: ", A_("Model"), pmod->ID);
	    } else {
		pprintf(prn, "\n%s %d: ", A_("Model"), pmod->ID);
	    }
	}
	break;
    }

    if (pmod->aux == AUX_VAR || pmod->aux == AUX_VECM) {
	;
    } else if (pmod->aux == AUX_SYS) {
	pprintf(prn, A_("%s, using observations %s%s%s"),
		A_(system_short_string(pmod)),
		startdate, datesep, enddate);
	maybe_print_T(pmod, dset, startdate, prn);
    } else if (!dataset_is_panel(dset)) {
	int mc, Tmax = pmod->t2 - pmod->t1 + 1;
	const char *estr = estimator_string(pmod, prn);
	const char *fmt;

	if (char_len(estr) > 32) {
	    fmt = N_("%s, obs. %s%s%s");
	    pprintf(prn, A_(fmt), A_(estr), startdate, 
		    (tex)? "--" : "-", enddate);
	    maybe_print_T(pmod, dset, startdate, prn);
	} else {
	    fmt = N_("%s, using observations %s%s%s");
	    pprintf(prn, A_(fmt), A_(estr), startdate, 
		    datesep, enddate);
	    maybe_print_T(pmod, dset, startdate, prn);
	}

	if (pmod->ci == HECKIT) {
	    mc = Tmax - gretl_model_get_int(pmod, "totobs");
	} else {
	    mc = Tmax - pmod->nobs;
	}

	if (mc > 0) {
	    gretl_prn_newline(prn);
	    pprintf(prn, "%s: %d", A_("Missing or incomplete observations dropped"), 
		    mc);
	}
    } else {
	/* panel data */
	int effn = gretl_model_get_int(pmod, "n_included_units");
	int Tmin = gretl_model_get_int(pmod, "Tmin");
	int Tmax = gretl_model_get_int(pmod, "Tmax");

	pprintf(prn, A_("%s, using %d observations"),
		A_(estimator_string(pmod, prn)), pmod->nobs);
	if ((pmod->opt & OPT_U) && (pmod->opt & OPT_N)) {
	    gretl_prn_newline(prn);
	    pprintf(prn, A_("Using Nerlove's transformation"));
	}
	if (effn > 0) {
	    gretl_prn_newline(prn);
	    pprintf(prn, A_("Included %d cross-sectional units"), effn);
	}
	if (Tmin > 0 && Tmax > 0) {
	    gretl_prn_newline(prn);
	    if (Tmin == Tmax) {
		pprintf(prn, A_("Time-series length = %d"), Tmin);
	    } else {
		pprintf(prn, A_("Time-series length: minimum %d, maximum %d"), 
			Tmin, Tmax);
	    }
	}
	if (pmod->ci == DPANEL) {
	    if (pmod->opt & OPT_L) {
		gretl_prn_newline(prn);
		pputs(prn, "Including equations in levels");
	    }
	    if (pmod->opt & OPT_X) {
		gretl_prn_newline(prn);
		pputs(prn, "H-matrix as per Ox/DPD");
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
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG || pmod->aux == AUX_ADD) {
	pprintf(prn, "%s: %s", A_("Dependent variable"),
		(tex)? "$\\hat{u}$" : "uhat");
    } else if (pmod->aux == AUX_WHITE || pmod->aux == AUX_HET_1) {
	pprintf(prn, "%s: %s", A_("Dependent variable"),
		(tex)? "$\\hat{u}^2$" : "uhat^2");
    } else if (pmod->aux == AUX_BP) {
	const char *fmt;

	if (pmod->opt & OPT_R) {
	    fmt = N_("scaled %s (Koenker robust variant)");
	} else {
	    fmt = N_("scaled %s");
	}
	pprintf(prn, "%s: ", A_("Dependent variable"));
	pprintf(prn, A_(fmt), (tex)? "$\\hat{u}^2$" : "uhat^2");
    } else if (pmod->aux == AUX_ARCH) {
	pprintf(prn, "%s: %s", A_("Dependent variable"),
		(tex)? "$u_t^2$" : "ut^2");
    } else if (pmod->ci == NLS || pmod->ci == MLE || pmod->ci == GMM) {
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
	print_arma_depvar(pmod, dset, prn);
    } else if (pmod->ci == INTREG) {
	print_intreg_depvar(pmod, dset, prn);
    } else if (pmod->ci == BIPROBIT) {
	dvnl = 0;
    } else { 
	const char *dvname = 
	    gretl_model_get_depvar_name(pmod, dset);

	if (tex) {
	    if (pmod->aux == AUX_VECM) {
		tex_vecm_depvar_name(vname, dvname);
	    } else if (pmod->ci == ARBOND || pmod->ci == DPANEL) {
		tex_dpd_depvar_name(vname, dvname);
	    } else {
		tex_escape(vname, dvname);
	    }
	}

	if (pmod->aux == AUX_VAR || pmod->aux == AUX_VECM) {
	    pputs(prn, (tex)? vname : dvname);
	} else {
	    pprintf(prn, "%s: %s", A_("Dependent variable"),
		    (tex)? vname : dvname);
	}
    }

    if (csv) pputc(prn, '"');

    if (dvnl) {
	gretl_prn_newline(prn);
    }

    /* supplementary strings below the estimator and sample info */

    if (pmod->ci == IVREG) {
	/* list of instruments for IV estimation */
	int method = gretl_model_get_int(pmod, "method");

	if (method != SYS_METHOD_FIML && method != SYS_METHOD_LIML) {
	    print_ivreg_instruments(pmod, dset, prn);
	}
    } else if (pmod->ci == LAD) {
	/* tau for quantile regression */
	double tau = gretl_model_get_double(pmod, "tau");

	if (!na(tau)) {
	    if (tex) {
		pprintf(prn, "$\\tau$ = %g", tau);
	    } else {
		pprintf(prn, "tau = %g", tau);
	    }
	    gretl_prn_newline(prn);
	}
    } else if (pmod->ci == HECKIT) {
	/* selection variable for Heckit */
	const char *selvar = heckit_selvar_name(pmod, dset);

	if (selvar != NULL) {
	    if (csv) pputc(prn, '"');
	    if (tex) {
		tex_escape(vname, selvar);
	    }
	    pprintf(prn, "%s: %s", A_("Selection variable"), 
		    (tex)? vname : selvar);
	    if (csv) pputc(prn, '"');
	    pputc(prn, '\n');	    
	}	
    } else if (pmod->ci == AR1) {
	if (pmod->arinfo != NULL) {
	    double r = pmod->arinfo->rho[0];

	    if (tex) {
		pprintf(prn, "$\\rho$ = %g", r);
	    } else {
		pprintf(prn, "rho = %g", r);
	    }
	    gretl_prn_newline(prn);
	}
    } else if (pmod->ci == PROBIT && (pmod->opt & OPT_E)) {
	int qp = gretl_model_get_int(pmod, "quadpoints");

	if (qp > 0) {
	    pprintf(prn, A_("Using %d quadrature points"), qp);
	    gretl_prn_newline(prn);
	}
    }

    /* VCV variants */
    print_model_vcv_info(pmod, dset, prn);

    if (pmod->ci == PANEL && (pmod->opt & OPT_W) && !pmod->aux) {
	/* WLS on panel data */
	if (tex) {
	    pputs(prn, "\\\\\n");
	}
	if (gretl_model_get_int(pmod, "iters")) {
	    pprintf(prn, A_("Allowing for groupwise heteroskedasticity"));
	} else {
	    pprintf(prn, A_("Weights based on per-unit error variances"));
	}
	pputc(prn, '\n');
    } else if (pmod->ci == WLS && !pmod->aux) {
	/* weight variable for WLS */
	if (tex) {
	    tex_escape(vname, dset->varname[pmod->nwt]);
	}
	if (csv) pputc(prn, '"');
	pprintf(prn, "%s: %s", A_("Variable used as weight"), 
		(tex)? vname : dset->varname[pmod->nwt]);
	if (csv) pputc(prn, '"');
	pputc(prn, '\n');
    } else if (pmod->ci == ARCH) {
	/* weight variable for ARCH */
	if (csv) pputc(prn, '"');
	pprintf(prn, "%s: %s", A_("Variable used as weight"), 
		(tex)? "$1/\\hat{\\sigma}_t$" : "1/sigma");
	if (csv) pputc(prn, '"');
	pputc(prn, '\n');
    } else if (pmod->ci == AR1) {
	/* rhohat for AR1 (TeX only) */
	if (tex) {
	    pprintf(prn, "$\\hat{\\rho}$ = %g\n", 
		    gretl_model_get_double(pmod, "rho_in"));
	}
    } else if (pmod->ci == LOGISTIC) {
	/* y-hat formula for logistic regression */
	if (tex) {
	    pprintf(prn, "$\\hat{y} = %g / (1 + e^{-X\\hat{\\beta}})$\n", 
		    gretl_model_get_double(pmod, "lmax"));  
	} else {
	    pprintf(prn, "yhat = %g / (1 + exp(-X*b))\n",  
		    gretl_model_get_double(pmod, "lmax"));
	}
    }

    /* IVREG: message about redundant instruments */
    if (plain_format(prn) && pmod->ci == IVREG &&
	gretl_model_get_data(pmod, "inst_droplist") != NULL) {
	print_ivreg_droplist(pmod, dset, prn);
    }  

    /* messages about collinear and/or zero regressors */
    if (plain_format(prn)) {
	if (gretl_model_get_data(pmod, "zerolist") != NULL) {
	    print_model_zerolist(pmod, dset, prn);
	}
	if (gretl_model_get_data(pmod, "droplist") != NULL) {
	    print_model_droplist(pmod, dset, prn);
	}
    } 

    if (plain_format(prn) && pmod->ci == LAD) {
	maybe_print_lad_warning(pmod, prn);
    } 

    if (plain_format(prn) && hessian_maybe_fishy(pmod)) {
	maybe_print_hessian_warning(pmod, prn);
    }    

    if (pmod->missmask == NULL && gretl_model_get_int(pmod, "wt_dummy")) { 
	/* FIXME alt formats */
	pprintf(prn, "%s %d\n", 
		A_("Weight var is a dummy variable, effective obs ="),
		pmod->nobs);
    } 

    if (rtf_format(prn)) {
	pputs(prn, "\\par\n");
    } else {
	pputc(prn, '\n');
    }
}

static int use_zscore (const MODEL *pmod)
{
    if (pmod == NULL) {
	/* modprint */
	return 1;
    } else if (0 && gretl_model_get_int(pmod, "dfcorr")) {
	/* override ASYMPTOTIC_MODEL in system case? */
	return 0;
    } else if (ASYMPTOTIC_MODEL(pmod->ci)) {
	return 1;
    } else {
	return 0;
    }
}

static void model_format_start (PRN *prn)
{
    set_alt_gettext_mode(prn);

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
    } 
}

#define RTF_MULTICOL  "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                      "\\cellx8000\n\\intbl"

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
   text print of a model. When using TeX, returns the number of
   columns in the table, otherwise just returns 0.
*/

static int alt_print_coeff_table_start (const MODEL *pmod, int ci, PRN *prn)
{
    const char *tlabel;
    int use_param = 0;
    int slopes = 0;
    int intervals = 0;
    int seqcols = 0;
    int mp = 0, ret = 0;

    if (use_zscore(pmod)) {
	tlabel = (tex_format(prn))? N_("$z$") : N_("z");
    } else {
	tlabel = (tex_format(prn))? N_("$t$-ratio") : N_("t-ratio");
    }

    if (pmod != NULL) {
	gretl_matrix *m;

	use_param = NONLIST_MODEL(pmod->ci);
	slopes = binary_model(pmod) && !(pmod->opt & OPT_P);
	intervals = gretl_model_get_data(pmod, "coeff_intervals") != NULL;
	m = gretl_model_get_data(pmod, "rq_sequence");
	seqcols = gretl_matrix_cols(m);
	mp = (pmod->ci == MPOLS);
    }

    if (csv_format(prn)) {
	char d = prn_delim(prn);

	if (mp) {
	    pprintf(prn, "%c\"%s\"%c\"%s\"\n",
		    d, A_("coefficient"), d, A_("std. error"));
	} else if (slopes) {
	    pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    d, A_("coefficient"), d, A_("std. error"),
		    d, A_(tlabel), d, A_("slope at mean"));
	} else if (use_param) {
	    pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    d, A_("estimate"), d, A_("std. error"),
		    d, A_(tlabel), d, A_("p-value"));
	} else if (intervals) {
	    pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    d, A_("coefficient"), d, A_("lower"),
		    d, A_("upper"));
	} else if (seqcols == 3) {
	    pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    d, "tau", d, A_("coefficient"), d, A_("lower"),
		    d, A_("upper"));
	} else if (seqcols == 2) {
	    pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    d, A_("coefficient"), d, A_("std. error"),
		    d, A_(tlabel));
	} else {
	    pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n",
		    d, A_("coefficient"), d, A_("std. error"),
		    d, A_(tlabel), d, A_("p-value"));
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
	    ret = tex_coeff_table_start(cols, tabopt, prn);
	} else if (rtf_format(prn)) {
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
		    pprintf(prn, " \\qc {\\i %s{\\super *}}\\cell", A_(cols[i]));
		} else {
		    pprintf(prn, " \\qc {\\i %s}\\cell", A_(cols[i]));
		}
	    }
	    if (!slopes && !intervals) {
		pputs(prn, " \\ql \\cell");
	    }
	    pputs(prn, " \\intbl \\row\n");
	} 
    }

    return ret;
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
	pprintf(prn, "\\vspace{1em}%s:\n\n", A_(msg[i]));
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"\n", A_(msg[i]));
    } else { 
	/* RTF */
	pprintf(prn, "\\par \\qc\n%s:\n\n", A_(msg[i]));	
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
	pprintf(prn, "\\par \\ql\n%s: TR{\\super 2} = %f,\n", A_("Test statistic"), 
		X);
	pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n", 
		A_("with p-value"), A_("Chi-square"), df, X, pv);
    } else if (tex_format(prn)) {
	pprintf(prn, "\n%s: $TR^2$ = %f,\n", A_("Test statistic"), X);
	pprintf(prn, "%s = $P$($\\chi^2(%d)$ > %f) = %f\n\n",
		A_("with p-value"), df, X, pv);
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
	pprintf(prn, "\\par \\ql\n%s: LM = %f,\n", A_("Test statistic"), 
		X);
	pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n", 
		A_("with p-value"), A_("Chi-square"), df, X, pv);
    } else if (tex_format(prn)) {
	pprintf(prn, "\n%s: LM = %f,\n", A_("Test statistic"), X);
	pprintf(prn, "%s = $P$($\\chi^2(%d)$ > %f) = %f\n\n",
		A_("with p-value"), df, X, pv);
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
	pprintf(prn, "\\par \\ql\n%s: HET_1 = %f,\n", A_("Test statistic"), z);
	pprintf(prn, "%s = 2 * P(z > %f) = %.3g\n\n", 
		A_("with p-value"), z, pv);
    } else if (tex_format(prn)) {
	pprintf(prn, "\n%s: \verb|HET_1| = %f,\n", A_("Test statistic"), z);
	pprintf(prn, "%s = $2 \times P$($z$ > %f) = %f\n\n",
		A_("with p-value"), z, pv);
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

#define fixed_effects_model(m) (m->ci == PANEL && (m->opt & OPT_F))

#define random_effects_model(m) (m->ci == PANEL && (m->opt & OPT_U))

#define between_model(m) (m->ci == PANEL && (m->opt & OPT_B))

#define weighted_model(m) (m->ci == HSK || m->ci == ARCH || \
			   (m->ci == WLS && !gretl_model_get_int(m, "wt_dummy")) || \
                           (m->ci == PANEL && (m->opt & OPT_W)))

#define panel_ML_model(m) (m->ci == PANEL && (m->opt & OPT_W) && \
			   gretl_model_get_int(m, "iters"))

#define non_weighted_panel(m) (m->ci == PANEL && !(m->opt & OPT_W))

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

static void set_mtab_string_width (struct middletab *mt)
{
    int len, badkey = 0, maxlen = 0;
    int i, j;

    for (i=0, j=0; i<7; i++, j+=2) {
	if (!na(mt->val[j])) {
	    len = g_utf8_strlen(_(mt->key[j]), -1);
	    if (len > maxlen) {
		maxlen = len;
		badkey = j;
	    }
	    len = g_utf8_strlen(_(mt->key[j+1]), -1);
	    if (len > maxlen) {
		maxlen = len;
		badkey = j+1;
	    }
	}
    }

    if (maxlen > 22) {
	fprintf(stderr, "Can't make compact model stats table -- the max\n"
		"length translated string is %d chars, should be < 23\n",
		maxlen);
	fprintf(stderr, "offending string: '%s' ->\n '%s'\n", 
		mt->key[badkey], _(mt->key[badkey]));
    }

    mt->mlen = maxlen;
}

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

static void mtab_numstart (char *s, double x, int d7, int minus)
{
    if (x < 0) {
	if (minus == MINUS_UTF) {
	    strcpy(s, (d7)? " −" : "−"); /* U+2212: minus */
	} else if (minus == MINUS_TEX) {
	    strcpy(s, "$-$");
	} else {
	    strcpy(s, (d7)? " -" : "-"); /* ASCII: use hyphen */
	}
    } else {
	strcpy(s, " "); 
    }
}

#define seven_digits(x) (x > 999999 && x < 1.0e8)

/* Try to pack as much as we can of a given number into a fixed width
   of 8 characters (a leading minus, if needed, is a ninth).
*/

static char *print_eight (char *s, struct middletab *mt, int i)
{
    double ax, x = mt->val[i];
    int d7 = (x < 0 && seven_digits(-x));
    char tmp[16];

    if (mt->ipos[i]) {
	sprintf(s, "%9d", (int) x);
	return s;
    }    

    if (na(x)) {
	sprintf(s, "%9s", "NA");
	return s;
    }

    if (i == RSQ_POS || i == RSQ_POS + 1) {
	/* R-squared: don't use scientific notation */
	sprintf(s, "%9.6f", x);
	return s;
    }    

    mtab_numstart(s, x, d7, mt->minus);

    ax = fabs(x);

    if (d7) {
	sprintf(tmp, "%.0f", ax);
    } else if (ax < 0.00001 || ax > 99999999) {
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
		    A_(s1), print_fifteen(x1, mt->val[j], mt->minus),
		    A_(s2), print_fifteen(x2, mt->val[k], mt->minus));
	} else {
	    pprintf(prn, "%s & %s & %s & %s \\\\\n",
		    A_(s1), print_eight(x1, mt, j),
		    A_(s2), print_eight(x2, mt, k));
	}
    } else if (rtf_format(prn)) {
	if (mt->multi) {
	    pputs(prn, RTF_MULTI_ROW);
	    pprintf(prn, RTF_MULTI_FMT, A_(s1), 
		    print_fifteen(x1, mt->val[j], mt->minus));
	    pputs(prn, RTF_MULTI_ROW);
	    pprintf(prn, RTF_MULTI_FMT, A_(s2), 
		    print_fifteen(x1, mt->val[k], mt->minus));
	} else {
	    pputs(prn, RTF_MT_ROW);
	    pprintf(prn, RTF_MT_FMT, 
		    A_(s1), print_eight(x1, mt, j),
		    A_(s2), print_eight(x2, mt, k));
	}
    } else if (csv_format(prn)) {
	pprintf(prn, "\"%s\"%c%s%c\"%s\"%c%s\n",
		A_(s1), mt->d, print_csv(x1, mt->val[j]), mt->d,
		A_(s2), mt->d, print_csv(x2, mt->val[k]));
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

static int string_is_translated (const char *s)
{
    if (strcmp(s, "Hannan-Quinn") && strcmp(s, "rho")) {
	return strcmp(s, _(s));
    } else {
	/* give the benefit of the doubt */
	return 1;
    }
}

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
	NULL,
	NULL,
	NULL,
	NULL,
	N_("Akaike information criterion"),
	N_("Schwarz Bayesian criterion"),
	N_("Hannan-Quinn Information Criterion"),
	NULL,
	NULL
    };
    int i;

    for (i=0; i<n; i++) {
	if (old_key[i] != NULL && !string_is_translated(S[i])) {
	    /* new-style string is not translated */
	    if (strcmp(old_key[i], _(old_key[i]))) {
		/* but the old-style one is, so we'll use it */
		S[i] = old_key[i];
	    }
	}
    }
}

/* print the block of statistics that appears beneath of the
   table of coefficients, standard errors, etc.
*/

static void print_middle_table (const MODEL *pmod, PRN *prn, int code)
{
    const char *note = 
	/* TRANSLATORS: please do not translate literally: this is for 
	   your use in describing locale-specific abbreviations. It is
	   OK to leave it untranslated, in which case it will not be
	   printed.
	*/
	N_("note on model statistics abbreviations here");
    int rtf = rtf_format(prn);
    int tex = tex_format(prn);
    int csv = csv_format(prn);
    char teststr[32];
    const char *key[] = {
	/* TRANSLATORS: maximum length of string is 22 characters */
	N_("Mean dependent var"),  /* Mean of dependent variable */
	/* TRANSLATORS: maximum length of string is 22 characters */
	N_("S.D. dependent var"),  /* Standard deviation of dependent var */
	/* TRANSLATORS: maximum length of string is 22 characters */
	N_("Sum squared resid"),   /* Sum of squared residuals */
	/* TRANSLATORS: maximum length of string is 22 characters */
	N_("S.E. of regression"),  /* Standard error of the regression */
	/* TRANSLATORS: maximum length of string is 22 characters */
	N_("R-squared"),           
	/* TRANSLATORS: maximum length of string is 22 characters */
	N_("Adjusted R-squared"),  
	"F-statistic",             /* will be replaced below */
	/* TRANSLATORS: maximum length of string is 22 characters */
	N_("P-value(F)"),          /* P-value of F-statistic */	
	/* TRANSLATORS: maximum length of string is 22 characters */
	N_("Log-likelihood"),      
	/* TRANSLATORS: maximum length of string is 22 characters */
	N_("Akaike criterion"),    /* Akaike Information Criterion */
	/* TRANSLATORS: maximum length of string is 22 characters */
	N_("Schwarz criterion"),   /* Schwarz Bayesian Criterion */
	/* TRANSLATORS: maximum length of string is 22 characters */
	N_("Hannan-Quinn"),        /* Hannan-Quinn Criterion */
	/* TRANSLATORS: maximum length of string is 22 characters */
	N_("rho"),                 /* 1st-order autocorrelation coeff. */
	/* TRANSLATORS: maximum length of string is 22 characters */
	N_("Durbin-Watson")        /* Durbin-Watson statistic */
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

    mtab.nls = doing_nls();
    if (mtab.nls) {
	maybe_remedy_translations(key, MID_STATS);
    }

    if (tex) {
	/* some special strings for TeX output */
	mtab.minus = MINUS_TEX;
	key[4] = "$R^2$";
	key[5] = N_("Adjusted $R^2$");
	key[7] = N_("P-value($F$)");
	key[11] = "Hannan--Quinn";
	key[12] = "$\\hat{\\rho}$";
	key[13] = "Durbin--Watson";
    } else if (!rtf && gretl_print_has_minus(prn)) {
	/* print a 'real' minus sign? */
	mtab.minus = MINUS_UTF;
    }

    if (pmod->aux == AUX_VECM || pmod->aux == AUX_COINT) {
	/* VECM equation or Engle-Granger test: suppress F-test */
	val[6] = val[7] = NADBL;
    } else if (!na(pmod->fstt)) {
	/* format F-stat and get its p-value */
	if (tex) {
	    sprintf(teststr, "$F(%d, %d)$", pmod->dfn, pmod->dfd);
	} else if (pmod->ci == PANEL && (pmod->opt & OPT_F)) {
	    sprintf(teststr, A_("LSDV F(%d, %d)"), pmod->dfn, pmod->dfd);
	} else {
	    sprintf(teststr, "F(%d, %d)", pmod->dfn, pmod->dfd);
	}
	key[6] = teststr;
	val[7] = snedecor_cdf_comp(pmod->dfn, pmod->dfd, pmod->fstt);
    } else if (!na(pmod->chisq)) {
	/* alternative: chi-square and its p-value */
	sprintf(teststr, "%s(%d)", A_("Chi-square"), pmod->dfn);
	key[6] = teststr;  
	val[6] = pmod->chisq;
	key[7] = N_("p-value");  
	val[7] = chisq_cdf_comp(pmod->dfn, val[6]);
    }

    /* special variants of R-squared */
    if (gretl_model_get_int(pmod, "uncentered")) {
	key[4] = (tex)? N_("Uncentered $R^2$") : 
	    N_("Uncentered R-squared"); /* 22: */
	key[5] = (tex)? N_("Centered $R^2$") : 
	    N_("Centered R-squared");  /* 22: */
	val[5] = gretl_model_get_double(pmod, "centered-R2");
    } else if (COUNT_MODEL(pmod->ci) || binary_model(pmod)) {
	key[4] = (tex)? N_("McFadden $R^2$") : 
	    N_("McFadden R-squared");  /* 22: McFadden's pseudo-R-squared */
    } else if (pmod->ci == PANEL && (pmod->opt & OPT_F)) {
	/* 22: panel, fixed effects */
	key[4] = (tex)? N_("LSDV $R^2$") : N_("LSDV R-squared");
	key[5] = (tex)? N_("Within $R^2$") : N_("Within R-squared");
    }

    if (pmod->ci == ARBOND || pmod->ci == DPANEL) {
	for (i=0; i<MID_STATS; i++) {
	    if (i < 2 || i > 3) {
		val[i] = NADBL;
	    }
	}
    } else if (pmod->ci == IVREG && (pmod->opt & OPT_G)) {
	/* IVREG via GMM */
	for (i=2; i<MID_STATS; i++) {
	    val[i] = NADBL;
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
	    for (i=6; i<MID_STATS; i++) {
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
    } else if (logit_probit_model(pmod)) {
	val[2] = val[3] = NADBL;
	val[6] = val[7] = NADBL;
	val[12] = val[13] = NADBL;
    } else if (pmod->ci == BIPROBIT) {
	val[2] = val[3] = NADBL;
	val[6] = val[7] = NADBL;
	val[12] = val[13] = NADBL;
    } else if (pmod->ci == HECKIT) {
	key[2] = (tex)? "$\\hat{\\sigma}$" : N_("sigma");
	val[2] = pmod->sigma;
	key[3] = (tex)? "$\\hat{\\rho}$" : N_("rho");
	val[3] = pmod->rho;
	for (i=4; i<MID_STATS; i++) {
	    /* may add an R-squared (item 4)? */
	    if (i < 8 || i > 11) {
		val[i] = NADBL;
	    }
	}	
    } else if (panel_ML_model(pmod) || pmod->ci == GARCH) {
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
    } else if (intreg_model(pmod)) {
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
	pputs(prn, "\\par\n{");
    } else if (csv) {
	mtab.d = prn_delim(prn);
    }

    mtab.key = key;
    mtab.val = val;
    strcpy(mtab.txt_fmt, TXT_MT_FMT);

    if (plain_format(prn) && mtab.nls) {
	set_mtab_string_width(&mtab);
    }

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
	pputc(prn, '\n');
	pprintf(prn, A_("Convergence achieved after %d iterations\n"), iters);
    } else {
	int fncount = gretl_model_get_int(pmod, "fncount");
	int grcount = gretl_model_get_int(pmod, "grcount");

	if (fncount > 0) {
	    pputc(prn, '\n');
	    pprintf(prn, A_("Function evaluations: %d\n"), fncount);
	    pprintf(prn, A_("Evaluations of gradient: %d\n"), grcount);
	}
    }
}

static void aux_print_info_criteria (const MODEL *pmod, PRN *prn)
{
    const char *istrs[] = {
	N_("AIC"), N_("BIC"), N_("HQC")
    };
    int i, n = 0;

    for (i=0; i<C_MAX; i++) {
	if (!na(pmod->criterion[i])) {
	    if (n > 0) {
		pputc(prn, ' ');
	    }	
	    pprintf(prn, "  %s: %g", A_(istrs[i]), pmod->criterion[i]);
	    n++;
	}
    }

    if (n > 0) {
	if (gretl_model_get_int(pmod, "eg-resids")) {
	    pputc(prn, '\n');
	} else {
	    pputs(prn, "\n\n");
	}
    }
}

static void set_csv_delim (PRN *prn)
{
    char test[4];

    sprintf(test, "%.1f", 1.0);

    if (test[1] == ',') {
	gretl_print_set_delim(prn, ';');
    } else {
	gretl_print_set_delim(prn, ',');
    }
}

/**
 * printmodel:
 * @pmod: pointer to gretl model.
 * @dset: data information struct.
 * @opt: may contain %OPT_O to print covariance matrix, %OPT_S
 * to get a "simple" print (just coefficients and standard
 * errors).
 * @prn: gretl printing struct.
 *
 * Print to @prn the estimates in @pmod plus associated statistics.
 * 
 * Returns: 0 on success, 1 if some of the values to print were %NaN.
 */

int printmodel (MODEL *pmod, const DATASET *dset, gretlopt opt, 
		PRN *prn)
{
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

    print_model_heading(pmod, dset, opt, prn);

    if (plain_format(prn)) {
	gotnan = plain_print_coefficients(pmod, dset, prn);
    } else {
	gotnan = alt_print_coefficients(pmod, dset, prn);
    }

    print_coeff_table_end(pmod, prn);

    if (pmod->aux == AUX_DF || pmod->aux == AUX_ADF || 
	pmod->aux == AUX_KPSS) {
	aux_print_info_criteria(pmod, prn);
	goto close_format;
    }

    if (pmod->aux == AUX_ARCH || pmod->aux == AUX_RESET) {
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
	       pmod->aux == AUX_AR ||
	       pmod->aux == AUX_ADD) {
	rsqline(pmod, prn);
	goto close_format;
    } else if (pmod->aux == AUX_COMFAC) {
	ssrline(pmod, prn);
	goto close_format;
    }

    if (opt & OPT_S) {
	/* --simple-print */
	if (pmod->ci == OLS && !na(pmod->rsq) && plain_format(prn)) {
	    pprintf(prn, "%s = %g, %s = %f\n\n", A_("SSR"), 
		    pmod->ess, A_("R-squared"), pmod->rsq);
	}
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

    /* random effects probit */
    if (re_probit_model(pmod)) {
	print_probit_rho(pmod, prn);
    }

    /* additional stats/info for some cases */
    if ((pmod->aux == AUX_SYS && liml_equation(pmod)) || 
	liml_model(pmod)) {
	print_liml_equation_data(pmod, prn);
    } else if (pmod->ci == ARMA) {
	print_arma_roots(pmod, prn);
    } else if (pmod->ci == GARCH) {
	garch_variance_line(pmod, prn);
    } else if (pmod->ci == HECKIT) {
	print_heckit_stats(pmod, prn);
    } else if (random_effects_model(pmod)) {
	panel_variance_lines(pmod, prn);
    } else if (gmm_model(pmod)) {
	print_GMM_stats(pmod, prn);
    } else if (pmod->ci == ARBOND || pmod->ci == DPANEL) {
	print_DPD_stats(pmod, prn);
    } else if (logit_probit_model(pmod)) {
	if (!pmod->aux) {
	    logit_probit_stats(pmod, prn);
	}
    } else if (tsls_model(pmod) && plain_format(prn)) {
	addconst_message(pmod, prn);
    } else if (intreg_model(pmod)) {
	print_intreg_info(pmod, dset, prn);
    } else if (pmod->ci == POISSON) {
	print_overdisp_test(pmod, prn);
    } else if (pmod->ci == DURATION) {
	print_duration_alpha(pmod, prn);
    } else if (pmod->ci == BIPROBIT) {
	print_probit_rho(pmod, prn);
    }

    /* FIXME alternate R^2 measures (within, centered) */

    if (plain_format(prn) && !pmod->aux) {
	maybe_print_jll(pmod, 0, prn);
    }

 pval_max:

    if (plain_format(prn) && pmod->ci != MLE && pmod->ci != PANEL &&
	pmod->ci != ARMA && pmod->ci != NLS && pmod->ci != GMM &&
	pmod->ci != LAD && pmod->ci != HECKIT && 
	pmod->ci != ARBOND && pmod->ci != DPANEL && pmod->ci != GARCH && 
	pmod->ci != DURATION && !ordered_model(pmod) && !multinomial_model(pmod) && 
	!COUNT_MODEL(pmod->ci) && !intreg_model(pmod) &&
	pmod->ci != BIPROBIT && !pmod->aux) {
	pval_max_line(pmod, dset, prn);
    }

 close_format:

    if (opt & OPT_V) {
	ols_print_anova(pmod, prn);
    }

    if (opt & OPT_O) {
	outcovmx(pmod, dset, prn);
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
prepare_model_coeff (const MODEL *pmod, const DATASET *dset,
		     int i, int adfnum, model_coeff *mc, PRN *prn)
{
    int gotnan = 0;

    model_coeff_init(mc);

    mc->show_pval = !binary_model(pmod) || (pmod->opt & OPT_P);

    if (tex_format(prn)) {
	make_tex_coeff_name(pmod, dset, i, mc->name);
    } else {
	gretl_model_get_param_name(pmod, dset, i, mc->name);
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
	pprintf(prn, " \\qc %s\\cell", A_("undefined"));
    } else {
	rtf_print_double(mc->b, prn);
    }

    if (!na(mc->lo) && !na(mc->hi)) {
	rtf_print_double(mc->lo, prn);
	rtf_print_double(mc->hi, prn);
	goto rtf_finish;
    }

    if (na(mc->se)) {
	pprintf(prn, " \\qc %s\\cell", A_("undefined"));
	pprintf(prn, " \\qc %s\\cell", A_("undefined"));
	pprintf(prn, " \\qc %s\\cell", A_("undefined"));
	goto rtf_finish;
    } 

    rtf_print_double(mc->se, prn); 

    if (!na(mc->tval)) {
	pprintf(prn, " \\qc %.4f\\cell", mc->tval);
    } else {
	pprintf(prn, " \\qc %s\\cell", A_("undefined"));
    }

    if (!na(mc->slope)) {
	rtf_print_double(mc->slope, prn);
    } else if (mc->show_pval) {
	if (na(mc->pval)) {
	    if (mc->df_pval) {
		pprintf(prn, " \\qc %s\\cell", A_("unknown"));
	    } else {
		pprintf(prn, " \\qc %s\\cell", A_("undefined"));
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
	pprintf(prn, "%c\"%s\"", d, A_("undefined"));
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
	pprintf(prn, "%c\"%s\"\n", d, A_("undefined"));
	return;
    }

    pprintf(prn, "%c%.15g", d, mc->se);

    if (!na(mc->tval)) {
	pprintf(prn, "%c%.15g", d, mc->tval);
    } else {
	pprintf(prn, "%c\"%s\"\n", d, A_("undefined"));
    }

    if (!na(mc->slope)) {
	/* slope for binary models */
	pprintf(prn, "%c%.15g", d, mc->slope);
    } else if (mc->show_pval) {
	if (na(mc->pval)) {
	    if (mc->df_pval) {
		pprintf(prn, "%c\"%s\"\n", d, A_("unknown"));
	    } else {
		pprintf(prn, "%c\"%s\"\n", d, A_("undefined"));
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

static void print_coeff_separator (const char *s, int n, PRN *prn)
{
    int havestr = (s != NULL && *s != '\0');

    if (plain_format(prn)) {
	if (havestr) {
	    if (n > 0) {
		pputs(prn, "\n  "); 
		print_centered(_(s), n, prn);
	    } else {
		pputs(prn, "  "); 
		pputs(prn, _(s));
	    }
	    pputc(prn, '\n'); 
	}
	if (n > 0) {
	    pputc(prn, '\n');
	}
    } else if (tex_format(prn)) {
	if (havestr) {
	    pputs(prn, "\\\\ [-8pt]\n");
	    pprintf(prn, "\\multicolumn{%d}{c}{%s} \\\\[1ex]\n", n, A_(s));
	} else {
	    pputs(prn, "\\\\ \n");
	}
    } else if (rtf_format(prn)) {
	pputs(prn, RTF_MULTICOL);
	if (havestr) {
	    pprintf(prn, "\\qc %s", A_(s));
	}
	pputs(prn, "\\cell\\intbl\\row\n");
    } else if (csv_format(prn)) {
	if (havestr) {
	    pprintf(prn, "\n\"%s\"\n", A_(s));
	} else {
	    pputc(prn, '\n');
	}
    }
}

static void print_coeff_left_string (const char *s, PRN *prn)
{
    if (plain_format(prn)) {
	pprintf(prn, " %s:\n", s);
    } else if (tex_format(prn)) {
	char tmp[32];

	tex_escape(tmp, s);
	pputs(prn, "\\\\ [-8pt]\n");
	
	pprintf(prn, "%s \\\\[1ex]\n", tmp);
    } else if (rtf_format(prn)) {
	pputs(prn, RTF_MULTICOL);
	pprintf(prn, "\\ql %s", s);
	pputs(prn, "\\cell\\intbl\\row\n");
    } else if (csv_format(prn)) {
	pprintf(prn, "\n\"%s\"\n", s);
    }
}

static int 
print_rq_sequence (const MODEL *pmod, const DATASET *dset, PRN *prn)
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
    const char *head;

    if (tauvec == NULL || B == NULL) {
	return E_DATA;
    }

    ntau = gretl_vector_get_length(tauvec);
    bcols = gretl_matrix_cols(B);

    for (i=2; i<=pmod->list[0]; i++) {
	n = char_len(dset->varname[pmod->list[i]]);
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
	pprintf(prn, "  %-*s  ", namelen, dset->varname[pmod->list[i+2]]);
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
alt_print_rq_sequence (const MODEL *pmod, const DATASET *dset, PRN *prn)
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
	pprintf(prn, "  %-*s  ", namelen, dset->varname[pmod->list[i+2]]);
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

static int print_sep_row (int namelen, int ncols, int *w,
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

    return n;
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
	headings[2] = N_("z");
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
				  const DATASET *dset, 
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

    minus = (gretl_print_has_minus(prn))? MINUS_UTF : MINUS_HYPHEN;

    for (i=0; i<nc; i++) {
	if (xna(b[i])) {
	    err = E_NAN;
	    goto bailout;
	}
	gretl_model_get_param_name(pmod, dset, i, names[i]);
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

    strings_array_free(names, nc);

    return err;
}

static const char *
get_col_heading (const char **S, int j, int slopes,
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
print_count_offset (const MODEL *pmod, const DATASET *dset, 
		    struct printval *val, int namelen, 
		    int colsep, int w, int lmax,
		    int addoff, PRN *prn)
{
    int offvar = gretl_model_get_int(pmod, "offset_var");

    if (offvar > 0) {
	char name[24];
	int n;

	sprintf(name, "log(%s)", dset->varname[offvar]);
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

static void mn_logit_coeffsep (char *sep, const MODEL *pmod, 
			       const DATASET *dset, int i)
{
    const char *vname = gretl_model_get_depvar_name(pmod, dset);
    const gretl_matrix *y = gretl_model_get_data(pmod, "yvals");
    int val = (y != NULL)? y->val[i] : i;

    sprintf(sep, "%s = %d", vname, val);
}

static int plain_print_coeffs (const MODEL *pmod, 
			       const DATASET *dset, 
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
    int seppos = -1, cblock = 0;
    int lmax[4] = {0};
    int rmax[4] = {0};
    int w[4], addoff[4] = {0};
    int hlen;
    double tval, pval = 0.0;
    int n, d, nc = pmod->ncoeff;
    int dfd = pmod->dfd;
    int show_slope, adfnum = -1;
    int intervals = 0;
    int dotlen, namelen = 0;
    int colsep = 2;
    int ncols = 4;
    int i, j, k;
    int err = 0;

    if (use_zscore(pmod)) {
	headings[2] = N_("z");
    } 

    if (pmod->ci == AR || pmod->ci == ARCH) {
	k = 0;
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
    } else if (NONLIST_MODEL(pmod->ci)) {
	headings[0] = N_("estimate");
    }

    if (err) {
	return err;
    }

    nc -= gretl_model_get_int(pmod, "skipdums");

    vals = allocate_printvals(nc, ncols);
    if (vals == NULL) {
	return E_ALLOC;
    }

    names = strings_array_new_with_length(nc, 32);
    if (names == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    show_slope = binary_model(pmod) && !(pmod->opt & OPT_P);
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
	} else if (multinomial_model(pmod)) {
	    cblock = gretl_model_get_int(pmod, "cblock");
	}
    }

    for (i=0; i<nc; i++) {
	if (xna(b[i])) {
	    err = E_NAN;
	    goto bailout;
	}
	gretl_model_get_param_name(pmod, dset, i, names[i]);
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
	    } else if (pmod->ci == DURATION && j > 1 && 
		       !strcmp(names[i], "sigma")) {
		/* suppress result for H0: sigma = 0 */
		vals[i][j].x = NADBL;
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
    } else if (namelen > NAMETRUNC) {
	namelen = NAMETRUNC;
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
    dotlen = print_sep_row(namelen, ncols, w, colsep, prn);

    /* biprobit special: name of first dependent variable */
    if (pmod->ci == BIPROBIT) {
	print_coeff_left_string(gretl_model_get_depvar_name(pmod, dset), 
				prn);
    }

    /* print row values */

    k = 0;
    for (i=0; i<nc; i++) {
	char tmp[NAMETRUNC];

	if (i == seppos) {
	    if (pmod->ci == BIPROBIT) {
		pputc(prn, '\n');
		print_coeff_left_string(sepstr, prn);
	    } else {
		print_coeff_separator(sepstr, dotlen, prn);
	    }
	} else if (cblock > 0 && i % cblock == 0) {
	    char mnlsep[32];

	    mn_logit_coeffsep(mnlsep, pmod, dset, ++k);
	    print_coeff_separator(mnlsep, 0, prn);
	}
	maybe_trim_varname(tmp, names[i]);
	pprintf(prn, "  %-*s", namelen, tmp);
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
    } else if (COUNT_MODEL(pmod->ci)) {
	print_count_offset(pmod, dset, &vals[0][0], namelen, 
			   colsep, w[0], lmax[0], addoff[0],
			   prn);
    }

 bailout:

    for (i=0; i<nc; i++) {
	free(vals[i]);
    }
    free(vals);

    strings_array_free(names, nc);
    free(xb);
    free(xse);

    return err;
}

static int 
plain_print_coefficients (const MODEL *pmod, const DATASET *dset, PRN *prn)
{
    if (pmod->ncoeff == 0) {
	return 0;
    } else if (gretl_model_get_data(pmod, "rq_sequence") != NULL) {
	return print_rq_sequence(pmod, dset, prn);
    } else if (pmod->ci == MPOLS) {
	return plain_print_mp_coeffs(pmod, dset, prn);
    } else {
	return plain_print_coeffs(pmod, dset, prn);
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
alt_print_count_offset (const MODEL *pmod, const DATASET *dset, PRN *prn)
{
    int offvar = gretl_model_get_int(pmod, "offset_var");

    if (offvar > 0) {
	char name[24];

	sprintf(name, "log(%s)", dset->varname[offvar]);

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
alt_print_coefficients (const MODEL *pmod, const DATASET *dset, PRN *prn)
{
    gretl_matrix *intervals = NULL;
    const char *sepstr = NULL;
    int seppos = -1;
    model_coeff mc;
    int adfnum = -1;
    int nc = pmod->ncoeff;
    int cols, gotnan = 0;
    int i, err = 0;

    if (gretl_model_get_data(pmod, "rq_sequence") != NULL) {
	pputs(prn, "Sorry, not implemented yet!\n");
	return 1;
    }

    cols = alt_print_coeff_table_start(pmod, pmod->ci, prn);

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

    if (pmod->ci == BIPROBIT) {
	print_coeff_left_string(gretl_model_get_depvar_name(pmod, dset), 
				prn);
    }

    for (i=0; i<nc; i++) {

	err = prepare_model_coeff(pmod, dset, i, adfnum, &mc, prn);
	if (err) gotnan = 1;

	if (i == seppos) {
	    if (pmod->ci == BIPROBIT) {
		print_coeff_left_string(sepstr, prn);
	    } else {		
		print_coeff_separator(sepstr, cols, prn);
	    }
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
    } else if (COUNT_MODEL(pmod->ci)) {
	alt_print_count_offset(pmod, dset, prn);
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
		 A_("Root"), i, rx, ix, mod, fr);
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
		 A_("Root"), i, rx, ix, mod, fr);
     }
}

static void root_start (const char *tag, PRN *prn)
{
    if (plain_format(prn)) {
	pprintf(prn, "  %s\n", _(tag));
    } else if (tex_format(prn)) {
	pprintf(prn, "%s \\\\ \n", A_(tag));
    } else if (rtf_format(prn)) {
	pputs(prn, RTF_ROOT_ROW);
	pprintf(prn, "\\ql %s\\cell\\ql \\cell\\ql \\cell\\ql \\cell\\ql \\cell"
		"\\ql\\cell \\intbl \\row\n", A_(tag));
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
		    A_("Real"), A_("Imaginary"), A_("Modulus"), A_("Frequency"));
	} else if (rtf_format(prn)) {
	    pputs(prn, "\n\\par\n{" RTF_ROOT_ROW);
	    pprintf(prn, "\\qr \\cell \\qc \\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell"
		    " \\qc {\\i %s}\\cell \\intbl \\row\n",
		    A_("Real"), A_("Imaginary"), A_("Modulus"), A_("Frequency"));
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
    double cenpc = (100.0 * cenobs) / totobs;

    ensure_vsep(prn);

    if (plain_format(prn)) {
	pprintf(prn, "%s: %d\n", _("Total observations"), totobs);
	pprintf(prn, "%s: %d (%.1f%%)\n", _("Censored observations"), 
		cenobs, cenpc);
	pputc(prn, '\n');
    } else if (rtf_format(prn)) {
	pprintf(prn, RTFTAB "%s: %d\n", A_("Total observations"), totobs);
	pprintf(prn, RTFTAB "%s: %d (%.1f%%)\n", A_("Censored observations"), 
		cenobs, cenpc);
    } else if (tex_format(prn)) {
	pprintf(prn, "%s: %d \\\\\n", A_("Total observations"), totobs);
	pprintf(prn, "%s: %d (%.1f\\%%) \\\\\n", A_("Censored observations"), 
		cenobs, cenpc);
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

static int limdep_df (const MODEL *pmod)
{
    int df = gretl_model_get_int(pmod, "lr_df");

    if (df == 0) {
	df = pmod->dfn;
    }

    return df;
}

static void logit_probit_stats (const MODEL *pmod, PRN *prn)
{
    const int *act_pred = NULL;
    int binary, slopes, correct = 0;
    double pc_correct;
    int df = 0;

    if ((pmod->opt & OPT_M) || gretl_model_get_int(pmod, "ordered")) {
	/* ordered logit/probit or multinomial logit */
	binary = slopes = 0;
    } else {
	binary = 1;
	slopes = !(pmod->opt & OPT_P);
    }

    /* for overall likelihood ratio test */
    if (!na(pmod->chisq)) {
	df = limdep_df(pmod);
    }

    if (binary) {
	act_pred = gretl_model_get_data(pmod, "discrete_act_pred");
	if (act_pred != NULL) {
	    correct = act_pred[0] + act_pred[3];
	}
    } else {	
	correct = gretl_model_get_int(pmod, "correct");
    }

    pc_correct = 100 * (double) correct / pmod->nobs;

    ensure_vsep(prn);

    if (plain_format(prn)) {
	if (correct > 0) {
	    pprintf(prn, "%s = %d (%.1f%%)\n", 
		    _("Number of cases 'correctly predicted'"), 
		    correct, pc_correct);
	}
	if (binary) {
	    double fXb = gretl_model_get_double(pmod, "fXb");

	    if (!na(fXb)) {
		pprintf(prn, "f(beta'x) %s = %.3f\n", _("at mean of independent vars"), 
			fXb);
	    }
	}
	if (df) {
	    pprintf(prn, "%s: %s(%d) = %g [%.4f]\n", 
		    _("Likelihood ratio test"), _("Chi-square"), 
		    df, pmod->chisq, chisq_cdf_comp(df, pmod->chisq));
	}
	pputc(prn, '\n');
	if (act_pred != NULL) {
	    plain_print_act_pred(act_pred, prn);
	}
    } else if (rtf_format(prn)) {
	pputc(prn, '\n');
	if (slopes) {
	    pprintf(prn, "\\par {\\super *}%s\n", A_("Evaluated at the mean"));
	}
	if (correct > 0) {
	    pprintf(prn, "\\par %s = %d (%.1f%%)\n", 
		    A_("Number of cases 'correctly predicted'"), 
		    correct, pc_correct);
	}
	if (binary) {
	    pprintf(prn, "\\par f(beta'x) %s = %.3f\n", A_("at mean of independent vars"), 
		    pmod->sdy);
	}
	if (df) {
	    pprintf(prn, "\\par %s: %s(%d) = %g [%.4f]\n",
		    A_("Likelihood ratio test"), A_("Chi-square"), 
		    df, pmod->chisq, chisq_cdf_comp(df, pmod->chisq));
	}
	pputc(prn, '\n');
    } else if (tex_format(prn)) {
	if (slopes) {
	    pprintf(prn, "\\begin{center}\n$^*$%s\n\\end{center}\n", 
		    A_("Evaluated at the mean"));
	}
	if (correct > 0 || df) {
	    pputs(prn, "\\vspace{1em}\n\\begin{raggedright}\n");
	    if (correct > 0) {
		pprintf(prn, "%s = %d (%.1f %s)\\\\\n", 
			A_("Number of cases `correctly predicted'"), 
			correct, pc_correct, A_("percent"));
	    }
	    if (df) {
		pprintf(prn, "%s: $\\chi^2(%d)$ = %.3f [%.4f]\\\\\n",
			A_("Likelihood ratio test"), 
			df, pmod->chisq, chisq_cdf_comp(df, pmod->chisq));
	    }
	    pputs(prn, "\\end{raggedright}\n");
	}
    }
}

int ols_print_anova (const MODEL *pmod, PRN *prn)
{
    double mst, msr, mse, rss;
    int n, c1, c2, c3;

    if (pmod->ci != OLS || !pmod->ifc ||
	na(pmod->ess) || na(pmod->tss)) {
	return E_NOTIMP;
    }

    pprintf(prn, "%s:\n\n", _("Analysis of Variance"));

    if (pmod->dfn == 0) {
	/* degenerate model: const only */
	rss = 0.0;
    } else {
	rss = pmod->tss - pmod->ess;
    }

    c1 = g_utf8_strlen(_("Sum of squares"), -1);
    c2 = g_utf8_strlen(_("df"), -1);
    c3 = g_utf8_strlen(_("Mean square"), -1);

    c1 = (c1 < 35)? 35 : c1;
    c2 = (c2 > 8)? c2 + 1 : (c2 < 8)? 8 : c2;
    c3 = (c3 > 16)? c3 + 1 : (c3 < 16)? 16 : c3;

    /* header strings are right-aligned */
    n = g_utf8_strlen(_("Sum of squares"), -1);
    bufspace(c1 - n, prn);
    pputs(prn, _("Sum of squares"));
    n = g_utf8_strlen(_("df"), -1);
    bufspace(c2 + 1 - n, prn);
    pputs(prn, _("df"));
    n = g_utf8_strlen(_("Mean square"), -1);
    bufspace(c3 + 1 - n, prn);
    pputs(prn, _("Mean square"));
    pputs(prn, "\n\n");
    c1 = 16;

    /* Mean Square, regression */
    msr = rss / pmod->dfn;
    /* string left-aligned with initial offset of 2 */
    n = g_utf8_strlen(_("Regression"), -1);
    bufspace(2, prn);
    pputs(prn, _("Regression"));
    bufspace(16 - n, prn);
    if (pmod->dfn == 0) {
	pprintf(prn, " %*g %*d %*s\n", c1, rss, c2, pmod->dfn, c3, _("undefined"));
    } else {
	pprintf(prn, " %*g %*d %*g\n", c1, rss, c2, pmod->dfn, c3, msr);
    }

    /* Mean Square, errors */
    mse = pmod->ess / pmod->dfd;
    /* string left-aligned with initial offset of 2 */
    n = g_utf8_strlen(_("Residual"), -1);
    bufspace(2, prn);
    pputs(prn, _("Residual"));
    bufspace(16 - n, prn);
    pprintf(prn, " %*g %*d %*g\n", c1, pmod->ess, c2, pmod->dfd, c3, mse);

    /* Mean Square, total */
    mst = pmod->tss / (pmod->nobs - 1);
    /* string left-aligned with initial offset of 2 */
    n = g_utf8_strlen(_("Total"), -1);
    bufspace(2, prn);
    pputs(prn, _("Total"));
    bufspace(16 - n, prn);
    pprintf(prn, " %*g %*d %*g\n", c1, pmod->tss, c2, pmod->nobs - 1, c3, mst);

    pprintf(prn, "\n  R^2 = %g / %g = %.6f\n", rss, pmod->tss, rss / pmod->tss);

    if (pmod->dfn == 0) {
	pprintf(prn, "  F(%d, %d) %s\n\n", pmod->dfn, pmod->dfd, _("undefined"));
	return 0;
    }

    if (pmod->ess == 0 || rss == 0.0) {
	pprintf(prn, "  F(%d, %d) = %g / %g (%s)\n\n", pmod->dfn, pmod->dfd, 
		msr, mse, _("undefined"));
    } else {
	double F = msr / mse;
	double pv = snedecor_cdf_comp(pmod->dfn, pmod->dfd, F);

	pprintf(prn, "  F(%d, %d) = %g / %g = %g ", 
		pmod->dfn, pmod->dfd, msr, mse, F);
	if (pv < .0001) {
	    pprintf(prn, "[%s %.3g]\n\n", _("p-value"), pv);
	} else if (!na(pv)) {
	    pprintf(prn, "[%s %.4f]\n\n", _("p-value"), pv); 
	}
    }

    return 0;
}

/**
 * print_model_from_matrices:
 * @cs: k x 2 matrix containing coefficients and standard errors.
 * @adds: matrix containing an additional p statistics, or NULL.
 * @s: string containing all required names, space-separated.
 * @prn: gretl printer.
 *
 * Prints to @prn the coefficient table and optional additional statistics
 * for a model estimated "by hand". Mainly useful for user-written functions.
 *
 * The number of space-separated "words" in the string @s must be k + p, where 
 * k is the number of coefficients and p the number of additional statistics 
 * given in @adds.
 * 
 * Returns: 0 on success, non-zero on failure.
 */

int print_model_from_matrices (const gretl_matrix *cs, 
			       const gretl_matrix *adds, 
			       const char *s, PRN *prn)
{
    int ncoef = gretl_matrix_rows(cs);
    int nadd = gretl_vector_get_length(adds);
    int ntot = ncoef + nadd;
    const char *sep = ",";
    char *tmp, **names = NULL;
    int i, err = 0;

    /* copy the user-defined string @s before applying strtok */
    tmp = gretl_strdup(s);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    names = malloc(ntot * sizeof *names);
    if (names == NULL) {
	free(tmp);
	return E_ALLOC;
    }

    if (strchr(s, ',') == NULL) {
	sep = " ";
    }

    for (i=0; i<ntot && !err; i++) {
	char *name = strtok((i == 0)? tmp : NULL, sep);

	if (name == NULL) {
	    gretl_errmsg_sprintf(_("modprint: expected %d names"), ntot);
	    err = E_DATA;
	} else {
	    while (isspace(*name)) {
		name++;
	    }
	    names[i] = name;
	}
    }

    if (!err) {
	const double *b = cs->val;
	const double *se = b + ncoef;

	set_alt_gettext_mode(prn);

	if (plain_format(prn)) {
	    /* newline here is useless for TeX and makes RTF choke */
	    pputc(prn, '\n'); 
	} else if (csv_format(prn)) {
	    set_csv_delim(prn);
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
    }

    free(names);
    free(tmp);

    return err;
}
