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

/*  guiprint.c - RTF and LaTeX generation for gretl, plus native
    printing */

#include "gretl.h"
#include "selector.h"
#include "winstack.h"
#include "textutil.h"
#include "treeutils.h"
#include "forecast.h"
#include "texprint.h"
#include "guiprint.h"
#include "gui_recode.h"
#include "graph_page.h"

#include "uservar.h"
#include "libset.h"
#include "gretl_xml.h"

#ifdef G_OS_WIN32
# include <windows.h>
# include "gretlwin32.h"
#else
# include "clipboard.h"
#endif

#ifdef __APPLE__
# include "osx_open.h"
#endif

#define PAGE_LINES 47

#ifdef G_OS_WIN32

#undef WGRDEBUG

int win32_print_graph (char *emfname)
{
    HENHMETAFILE hemf;
    HDC dc;
    PRINTDLG pdlg;
    DOCINFO di;
    int printok;
# ifdef WGRDEBUG
    FILE *fp = fopen("debug.txt", "w");
# endif

    hemf = GetEnhMetaFile(emfname);
    if (hemf == NULL) {
	file_read_errbox(emfname);
	return 1;
    }

    memset(&pdlg, 0, sizeof pdlg);
    pdlg.lStructSize = sizeof pdlg;
    pdlg.Flags = PD_RETURNDC | PD_NOPAGENUMS;

    printok = PrintDlg(&pdlg);
    if (!printok) {
	/* canceled */
	DeleteEnhMetaFile(hemf);
	return 0;
    }

    dc = pdlg.hDC;

    memset(&di, 0, sizeof di);
    di.cbSize = sizeof di;
    di.lpszDocName = "gretl";

    printok = StartDoc(dc, &di);

    if (printok) {
	RECT rect;
	float hfrac = 0.8, vfrac;
	float hpx, vpx;
	float hppi, vppi;
	float hsize, vsize;
	float hmarg, vmarg;

	StartPage(dc);

	hpx = (float) GetDeviceCaps(dc, HORZRES);
	vpx = (float) GetDeviceCaps(dc, VERTRES);
	hppi = (float) GetDeviceCaps(dc, LOGPIXELSX);
	vppi = (float) GetDeviceCaps(dc, LOGPIXELSY);

	hsize = hfrac * hpx;
	hmarg = ((1.0 - hfrac) / 2.0) * hpx;

	vsize = hsize * 0.75 * (vppi / hppi);
	vfrac = vsize / vpx;
	vmarg = ((1.0 - vfrac) / 3.0) * vpx;

	rect.left = (long) hmarg;
	rect.top = (long) vmarg;
	rect.right = (long) (hmarg + hsize);
	rect.bottom = (long) (vmarg + vsize);

# ifdef WGRDEBUG
	fprintf(fp, "hpx=%g, vpx=%g\n", hpx, vpx);
	fprintf(fp, "hsize=%g, vsize=%g\n", hsize, vsize);
	fprintf(fp, "hfrac=%g, vfrac=%g\n", hfrac, vfrac);
	fprintf(fp, "rect = %ld, %ld, %ld, %ld\n",
		rect.left, rect.top, rect.right, rect.bottom);
	fclose(fp);
# endif

	PlayEnhMetaFile(dc, hemf, &rect);
	printok = (EndPage(dc) > 0);
    }

    if (printok) {
        EndDoc(dc);
    } else {
        AbortDoc(dc);
    }

    DeleteDC(dc);
    GlobalFree(pdlg.hDevMode);
    GlobalFree(pdlg.hDevNames);

    DeleteEnhMetaFile(hemf);

    return !printok;
}

#endif /* G_OS_WIN32 */

#define GRETL_PNG_TMP "gretltmp.png"

void rtf_print_obs_marker (int t, const DATASET *pdinfo, PRN *prn)
{
    char tmp[OBSLEN] = {0};
    const char *obs;

    if (pdinfo->markers) {
	obs = pdinfo->S[t];
    } else {
	ntolabel(tmp, t, pdinfo);
	obs = tmp;
    }

    pprintf(prn, "\\intbl \\ql %s\\cell", obs);
}

/* row format specifications for RTF "tables" */

#define STATS_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
    "\\cellx2700\\cellx4000\\cellx6700\\cellx8000\n\\intbl"

static void printf_rtf (double x, PRN *prn, int endrow)
{
    if (na(x)) {
	if (endrow) {
	    pprintf(prn, "\\qc %s\\cell\\intbl \\row\n",
		    _("undefined"));
	} else {
	    pprintf(prn, "\\qc %s\\cell", _("undefined"));
	}
    } else if (endrow) {
	pprintf(prn, "\\qc %#.*g\\cell \\intbl \\row\n",
		get_gretl_digits(), x);
    } else {
	pprintf(prn, "\\qc %#.*g\\cell", get_gretl_digits(), x);
    }
}

static void printk_rtf (int k, PRN *prn, int endrow)
{
    if (endrow) {
	pprintf(prn, "\\qc %d\\cell\\intbl \\row\n", k);
    } else {
	pprintf(prn, "\\qc %d\\cell", k);
    }
}

#define SUMM_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262"	\
    "\\cellx1600\\cellx3200\\cellx4800\\cellx6400"			\
    "\\cellx8000\n"

#define VAR_SUMM_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262"	\
    "\\cellx2000\\cellx4000\\cellx6000\\cellx8000\n"

#define SUMM_ROW_S  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262"	\
    "\\cellx1600\\cellx2800\\cellx4000\\cellx5200"			\
    "\\cellx6400\\cellx7200\n"

static void
rtfprint_simple_summary (const Summary *summ, const DATASET *pdinfo, PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN];
    int save_digits = get_gretl_digits();
    int i, vi;

    ntolabel(date1, pdinfo->t1, pdinfo);
    ntolabel(date2, pdinfo->t2, pdinfo);

    pputs(prn, "{\\rtf1\\par\n\\qc ");
    pprintf(prn, _("Summary Statistics, using the observations %s - %s"),
	    date1, date2);
    pputs(prn, "\\par\n");

    if (summary_has_missing_values(summ)) {
	pprintf(prn, "%s\\par\n\n", _("(missing values were skipped)"));
    }
    pprintf(prn, "{" SUMM_ROW_S
	    "\\intbl \\qc %s\\cell", _("Variable"));

    pprintf(prn,
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    _("Mean"), _("Median"), _("S.D."), _("Min"), _("Max"));

    set_gretl_digits(3);

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[vi]);
	printf_rtf(summ->mean[i], prn, 0);
	printf_rtf(summ->median[i], prn, 0);
	printf_rtf(summ->sd[i], prn, 0);
	printf_rtf(summ->low[i], prn, 0);
	printf_rtf(summ->high[i], prn, 1);
    }

    set_gretl_digits(save_digits);

    pputs(prn, "}}\n");
}

static void
rtfprint_summary_full (const Summary *summ, const DATASET *pdinfo, PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN];
    int save_digits = get_gretl_digits();
    int i, vi;

    ntolabel(date1, pdinfo->t1, pdinfo);
    ntolabel(date2, pdinfo->t2, pdinfo);

    pputs(prn, "{\\rtf1\\par\n\\qc ");
    pprintf(prn, _("Summary Statistics, using the observations %s - %s"),
	    date1, date2);
    pputs(prn, "\\par\n");

    if (summ->list[0] == 1) {
	pprintf(prn, _("for the variable %s (%d valid observations)"),
		pdinfo->varname[summ->list[1]], summ->n);
	pputs(prn, "\\par\n\n");
	pputs(prn, "{" VAR_SUMM_ROW "\\intbl ");
    } else {
	if (summary_has_missing_values(summ)) {
	    pputs(prn, _("(missing values were skipped)"));
	    pputs(prn, "\\par\n\n");
	}
	pprintf(prn, "{" SUMM_ROW
		"\\intbl \\qc %s\\cell", _("Variable"));
    }

    if (save_digits > 5) {
	set_gretl_digits(5);
    }

    pprintf(prn,
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    _("Mean"), _("Median"), _("Minimum"), _("Maximum"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[vi]);
	}
	printf_rtf(summ->mean[i], prn, 0);
	printf_rtf(summ->median[i], prn, 0);
	printf_rtf(summ->low[i], prn, 0);
	printf_rtf(summ->high[i], prn, 1);
    }

    if (summ->list[0] > 1) pprintf(prn, "\\intbl \\qc %s\\cell",
				   _("Variable"));

    pprintf(prn,
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    _("Std. Dev."), _("C.V."), _("Skewness"), _("Ex. kurtosis"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[vi]);
	}
	printf_rtf(summ->sd[i], prn, 0);
	printf_rtf(summ->cv[i], prn, 0);
	printf_rtf(summ->skew[i], prn, 0);
	printf_rtf(summ->xkurt[i], prn, 1);
    }

    if (summ->list[0] > 1) pprintf(prn, "\\intbl \\qc %s\\cell",
				   _("Variable"));

    pprintf(prn,
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    _("5% Perc."), _("95% Perc."), _("IQ range"), _("Missing obs."));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[vi]);
	}
	printf_rtf(summ->perc05[i], prn, 0);
	printf_rtf(summ->perc95[i], prn, 0);
	printf_rtf(summ->iqr[i], prn, 0);
	printk_rtf(summ->misscount[i], prn, 1);
    }

    set_gretl_digits(save_digits);

    pputs(prn, "}}\n");
}

int text_equation_ok (const MODEL *pmod)
{
    if (pmod->ci == OLS || pmod->ci == WLS ||
	pmod->ci == HSK || pmod->ci == AR1 ||
	pmod->ci == ARCH || pmod->ci == IVREG ||
	pmod->ci == PANEL || pmod->ci == LOGIT ||
	pmod->ci == PROBIT || pmod->ci == TOBIT ||
	pmod->ci == LOGISTIC) {
	return 1;
    } else if (pmod->ci == LAD) {
	if (gretl_model_get_data(pmod, "rq_tauvec") != NULL) {
	    /* multi-tau estimation: can't show equation */
	    return 0;
	} else {
	    return 1;
	}
    } else {
	return 0;
    }
}

static char *eqn_numstr (double x, char *s)
{
    sprintf(s, "%#.3g", x);
    return gretl_fix_exponent(s);
}

int text_print_equation (const MODEL *pmod, const DATASET *pdinfo,
			 gretlopt opt, PRN *prn)
{
    double x;
    char vname[32], xstr[12];
    int pad, c, cols[6] = {0};
    int k, totk = pmod->ncoeff;
    int ii, ii0, rem, maxper = 5;
    int n_lines, i, j;

    /* dependent variable */
    pputc(prn, '\n');
    c = pprintf(prn, "^%s = ", gretl_model_get_depvar_name(pmod, pdinfo));

    /* how many lines are needed, at @maxper coeffs per line? */
    n_lines = totk / maxper + (totk % maxper ? 1 : 0);

    /* coeffs remaining to print */
    rem = totk;

    /* index of first coeff on line */
    ii0 = 0;

    for (j=0; j<n_lines; j++) {
	k = rem > maxper ? maxper : rem;
	rem -= k; /* reduce for bext iteration */
	ii = ii0;

	/* coefficients times indep vars */
	for (i=0; i<k; i++) {
	    cols[i] = (i == 0)? (c - 1) : (c + 2); /* FIXME */
	    if (ii == 0 && pmod->ifc) {
		eqn_numstr(pmod->coeff[ii], xstr);
		c += pputs(prn, xstr);
	    } else {
		eqn_numstr(fabs(pmod->coeff[ii]), xstr);
		c += pprintf(prn, " %c %s", (pmod->coeff[ii] < 0.0)? '-' : '+',
			     xstr);
		gretl_model_get_param_name(pmod, pdinfo, ii, vname);
		c += pprintf(prn, "*%s", vname);
	    }
	    ii++;
	}
	pputc(prn, '\n');

	/* find our starting point in line @j */
	c = cols[0];
	bufspace(cols[0], prn);
	ii = ii0;

	/* standard errors or t-stats */
	for (i=0; i<k; i++) {
	    if (i == 0 && j > 0) {
		/* starting second or subsequent line */
		bufspace(3, prn);
		c += 3;
	    } else if (i > 0) {
		pad = cols[i] - c;
		if (pad > 0) {
		    bufspace(pad, prn);
		    c += pad;
		}
	    }
	    if (na(pmod->sderr[ii])) {
		c += pprintf(prn, "(NA)");
	    } else if (opt & OPT_T) {
		x = pmod->coeff[ii] / pmod->sderr[ii];
		eqn_numstr(x, xstr);
		c += pprintf(prn, "(%s)", xstr);
	    } else {
		eqn_numstr(pmod->sderr[ii], xstr);
		c += pprintf(prn, "(%s)", xstr);
	    }
	    ii++;
	}

	/* prepare for next line? */
	if (j < n_lines - 1) {
	    int p;

	    pputs(prn, "\n\n");
	    ii0 += k;
	    cols[0] = 3;
	    for (p=0; p<6; p++) {
		cols[p] = 0;
	    }
	    c = pputs(prn, "  ");
	}
    }

    pputs(prn, "\n\n");

    if (dataset_is_time_series(pdinfo)) {
	pprintf(prn, "T = %d", pmod->nobs);
    } else {
	pprintf(prn, "n = %d", pmod->nobs);
    }

    /* additional info (R^2 etc) */
    if (pmod->ci == LAD) {
	x = gretl_model_get_double(pmod, "ladsum");
	if (!na(x)) {
	    eqn_numstr(x, xstr);
	    pprintf(prn, ", %s = %s ", _("sum of abs. residuals"), xstr);
	}
    } else {
	if (!na(pmod->adjrsq)) {
	    pprintf(prn, ", %s = %.3f ", _("R-squared"), pmod->rsq);
	} else if (!na(pmod->lnL)) {
	    eqn_numstr(pmod->lnL, xstr);
	    pprintf(prn, ", %s = %s ", _("loglikelihood"), xstr);
	}
	x = gretl_model_get_double(pmod, "rho_gls");
	if (!na(x)) {
	    pprintf(prn, ", rho = %.3f", x);
	}
    }

    pprintf(prn, "\n(%s)\n",
	    (opt & OPT_T)? _("t-ratios in parentheses") :
	    _("standard errors in parentheses"));

    return 0;
}

int text_print_x_y_fitted (int vx, int vy, const double *f,
			   const DATASET *dset, PRN *prn)
{
    char obs1[OBSLEN], obs2[OBSLEN];
    char label[VNAMELEN];
    const double *x = dset->Z[vx];
    const double *y = dset->Z[vy];
    int obslen = max_obs_marker_length(dset);
    int t1 = dset->t1;
    int t2 = dset->t2;
    int t, pmaxx, pmaxy;
    int err = 0;

    for (t=t1; t<=dset->t2; t++) {
	if (na(x[t])) {
	    t1++;
	} else {
	    break;
	}
    }

    for (t=t2; t>dset->t1; t--) {
	if (na(x[t])) {
	    t2--;
	} else {
	    break;
	}
    }

    ntolabel(obs1, t1, dset);
    ntolabel(obs2, t2, dset);
    pprintf(prn, _("Model estimation range: %s - %s"), obs1, obs2);
    pputs(prn, "\n\n");
    bufspace(obslen, prn);

    for (t=0; t<3; t++) {
	if (t == 0) strcpy(label, dset->varname[vx]);
	if (t == 1) strcpy(label, dset->varname[vy]);
	if (t == 2) strcpy(label, _("fitted"));
	pprintf(prn, "%*s", UTF_WIDTH(label, 13), label);
    }

    pputs(prn, "\n\n");

    pmaxx = get_precision(x, dset->n, 6);
    pmaxy = get_precision(y, dset->n, 6);

    for (t=t1; t<=t2; t++) {
	print_obs_marker(t, dset, obslen, prn);
	if (na(x[t])) {
	    /* nothing to print */
	    pputc(prn, '\n');
	} else if (na(y[t])) {
	    /* y missing but x and fitted should be OK */
	    if (pmaxx == PMAX_NOT_AVAILABLE || pmaxy == PMAX_NOT_AVAILABLE) {
		pprintf(prn, "%13g%13s%13g\n", x[t], "NA", f[t]);
	    } else {
		pprintf(prn, "%13.*f%13s%13.*f\n", pmaxx, x[t], "NA", pmaxy, f[t]);
	    }
	} else {
	    /* got all values */
	    if (pmaxx == PMAX_NOT_AVAILABLE || pmaxy == PMAX_NOT_AVAILABLE) {
		pprintf(prn, "%13g%13g%13g\n", x[t], y[t], f[t]);
	    } else {
		pprintf(prn, "%13.*f%13.*f%13.*f\n",
			pmaxx, x[t], pmaxy, y[t], pmaxy, f[t]);
	    }
	}
    }

    pputc(prn, '\n');

    return err;
}

/* print value in (non-correlation) matrix */

static void tex_matnum (double x, PRN *prn)
{
    char s[32];

    tex_sprint_double_digits(x, s, 5);
    pprintf(prn, "%s & ", s);
}

static void printf_tex (double x, PRN *prn, int endrow)
{
    if (na(x)) {
	if (endrow) {
	    pprintf(prn, "\\multicolumn{2}{c}{%s}\\\\", _("undefined"));
	} else {
	    pprintf(prn, "\\multicolumn{2}{c}{%s} & ", _("undefined"));
	}
    } else {
	char s[32];

	tex_rl_double(x, s);
	if (endrow) {
	    pprintf(prn, "%s\\\\", s);
	} else {
	    pprintf(prn, "%s & ", s);
	}
    }
}

static void printk_tex (int k, PRN *prn, int endrow)
{
    if (endrow) {
	pprintf(prn, "\\multicolumn{2}{c}{%d}\\\\", k);
    } else {
	pprintf(prn, "\\multicolumn{2}{c}{%d} & ", k);
    }
}

static void
texprint_simple_summary (const Summary *summ, const DATASET *pdinfo, PRN *prn)
{
    char pt = get_local_decpoint();
    char date1[OBSLEN], date2[OBSLEN], vname[2*VNAMELEN];
    int save_digits = get_gretl_digits();
    int i, vi;

    ntolabel(date1, pdinfo->t1, pdinfo);
    ntolabel(date2, pdinfo->t2, pdinfo);

    pputs(prn, "\\begin{center}\n");
    pprintf(prn, _("Summary Statistics, using the observations %s--%s"),
	    date1, date2);
    pputs(prn, "\\\\\n");

    if (summary_has_missing_values(summ)) {
	pprintf(prn, "%s\\\\[8pt]\n\n", _("(missing values were skipped)"));
    } else {
	pputs(prn, "\n\\vspace{8pt}\n\n");
    }
    pprintf(prn, "\\begin{tabular}{lr@{%c}lr@{%c}lr@{%c}lr@{%c}lr@{%c}l}\n",
	    pt, pt, pt, pt, pt);
    pprintf(prn, "%s &", _("Variable"));

    pprintf(prn, " \\multicolumn{2}{c}{%s}\n"
	    " & \\multicolumn{2}{c}{%s}\n"
	    "  & \\multicolumn{2}{c}{%s}\n"
	    "   & \\multicolumn{2}{c}{%s}\n"
	    "    & \\multicolumn{2}{c}{%s} \\\\[1ex]\n",
	    _("Mean"), _("Median"), _("S.D."), _("Min"), _("Max"));

    set_gretl_digits(3);

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	tex_escape(vname, pdinfo->varname[vi]);
	pprintf(prn, "%s & ", vname);
	printf_tex(summ->mean[i], prn, 0);
	printf_tex(summ->median[i], prn, 0);
	printf_tex(summ->sd[i], prn, 0);
	printf_tex(summ->low[i], prn, 0);
	printf_tex(summ->high[i], prn, 1);
	pputc(prn, '\n');
    }

    set_gretl_digits(save_digits);

    pputs(prn, "\\end{tabular}\n\\end{center}\n");
}

static void
texprint_summary_full (const Summary *summ, const DATASET *pdinfo, PRN *prn)
{
    char pt = get_local_decpoint();
    char date1[OBSLEN], date2[OBSLEN];
    char vname[VNAMELEN*2];
    int save_digits = get_gretl_digits();
    int i, vi;

    ntolabel(date1, pdinfo->t1, pdinfo);
    ntolabel(date2, pdinfo->t2, pdinfo);

    pputs(prn, "\\begin{center}\n");
    pprintf(prn, _("Summary Statistics, using the observations %s--%s"),
	    date1, date2);
    pputs(prn, "\\\\\n");

    if (summ->list[0] == 1) {
	tex_escape(vname, pdinfo->varname[summ->list[1]]);
	pprintf(prn, _("for the variable %s (%d valid observations)"),
		vname, summ->n);
	pputs(prn, "\\\\[8pt]\n\n");
	pprintf(prn, "\\begin{tabular}{r@{%c}lr@{%c}lr@{%c}lr@{%c}l}\n",
		pt, pt, pt, pt);
    } else {
	if (summary_has_missing_values(summ)) {
	    pprintf(prn, "%s\\\\[8pt]\n\n", _("(missing values were skipped)"));
	} else {
	    pputs(prn, "\n\\vspace{8pt}\n\n");
	}
	pprintf(prn, "\\begin{tabular}{lr@{%c}lr@{%c}lr@{%c}lr@{%c}l}\n",
		pt, pt, pt, pt);
	pprintf(prn, "%s &", _("Variable"));
    }

    if (save_digits > 5) {
	set_gretl_digits(5);
    }

    pprintf(prn, " \\multicolumn{2}{c}{%s}\n"
	    " & \\multicolumn{2}{c}{%s}\n"
	    "  & \\multicolumn{2}{c}{%s}\n"
	    "   & \\multicolumn{2}{c}{%s} \\\\[1ex]\n",
	    _("Mean"), _("Median"), _("Minimum"), _("Maximum"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    tex_escape(vname, pdinfo->varname[vi]);
	    pprintf(prn, "%s & ", vname);
	}
	printf_tex(summ->mean[i], prn, 0);
	printf_tex(summ->median[i], prn, 0);
	printf_tex(summ->low[i], prn, 0);
	printf_tex(summ->high[i], prn, 1);
	if (i == summ->list[0] - 1) {
	    pputs(prn, "[10pt]\n\n");
	} else {
	    pputc(prn, '\n');
	}
    }

    if (summ->list[0] > 1) {
	pprintf(prn, "%s & ", _("Variable"));
    }

    pprintf(prn, " \\multicolumn{2}{c}{%s}\n"
	    " & \\multicolumn{2}{c}{%s}\n"
	    "  & \\multicolumn{2}{c}{%s}\n"
	    "   & \\multicolumn{2}{c}{%s} \\\\[1ex]\n",
	    _("Std.\\ Dev."), _("C.V."), _("Skewness"), _("Ex.\\ kurtosis"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    tex_escape(vname, pdinfo->varname[vi]);
	    pprintf(prn, "%s & ", vname);
	}
	printf_tex(summ->sd[i], prn, 0);
	printf_tex(summ->cv[i], prn, 0);
	printf_tex(summ->skew[i], prn, 0);
	printf_tex(summ->xkurt[i], prn, 1);
	if (i == summ->list[0] - 1) {
	    pputs(prn, "[10pt]\n\n");
	} else {
	    pputc(prn, '\n');
	}
    }

    if (summ->list[0] > 1) {
	pprintf(prn, "%s & ", _("Variable"));
    }

    pprintf(prn, " \\multicolumn{2}{c}{%s}\n"
	    " & \\multicolumn{2}{c}{%s}\n"
	    "  & \\multicolumn{2}{c}{%s}\n"
	    "   & \\multicolumn{2}{c}{%s} \\\\[1ex]\n",
	    /* xgettext:no-c-format */
	    _("5\\% perc."),
	    /* xgettext:no-c-format */
	    _("95\\% perc."),
	    _("IQ Range"),
	    _("Missing obs."));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    tex_escape(vname, pdinfo->varname[vi]);
	    pprintf(prn, "%s & ", vname);
	}
	printf_tex(summ->perc05[i], prn, 0);
	printf_tex(summ->perc95[i], prn, 0);
	printf_tex(summ->iqr[i], prn, 0);
	printk_tex(summ->misscount[i], prn, 1);
	pputc(prn, '\n');
    }

    set_gretl_digits(save_digits);

    pputs(prn, "\\end{tabular}\n\\end{center}\n");
}

void special_print_summary (const Summary *summ, const DATASET *pdinfo,
			    PRN *prn)
{
    if (tex_format(prn)) {
	if (summ->opt & OPT_S) {
	    texprint_simple_summary(summ, pdinfo, prn);
	} else {
	    texprint_summary_full(summ, pdinfo, prn);
	}
    } else if (rtf_format(prn)) {
	if (summ->opt & OPT_S) {
	    rtfprint_simple_summary(summ, pdinfo, prn);
	} else {
	    rtfprint_summary_full(summ, pdinfo, prn);
	}
    }
}

static void tex_outxx (double xx, PRN *prn)
{
    if (na(xx)) {
	pprintf(prn, "%s & ", _("undefined"));
    } else {
	pprintf(prn, "$%.4f$ & ", xx);
    }
}

static void rtf_outxx (double xx, PRN *prn)
{
    if (na(xx)) {
	pprintf(prn, "\\qc %s\\cell ", _("undefined"));
    } else {
	pprintf(prn, "\\qc %.4f\\cell ", xx);
    }
}

static void rtf_vmat_row (int lo, PRN *prn)
{
    int i, w = 1400;
    int cmax = (lo + 1 > 6)? 6 : lo + 1;

    pputs(prn, "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262");

    for (i=1; i<=cmax; i++) {
	pprintf(prn, "\\cellx%d", w * i);
    }

    pputs(prn, "\n\\intbl ");
}

static void rtf_table_pad (int pad, PRN *prn)
{
    while (pad--) pputs(prn, "\\cell ");
}

static void rtf_vmat_blank_row (int lo, int n, PRN *prn)
{
    rtf_vmat_row(lo, prn);
    while (n--) pputs(prn, "\\cell ");
    pputs(prn, "\\intbl \\row\n");
}

#define FIELDS 5

static void
rtfprint_vmatrix (const VMatrix *vmat, const DATASET *pdinfo, PRN *prn)
{
    register int i, j;
    int n = vmat->t2 - vmat->t1 + 1;
    int blockmax = vmat->dim / FIELDS;
    int nf, li2, p, k, idx, ij2;

    if (vmat->ci == CORR) {
	char date1[OBSLEN], date2[OBSLEN];

	ntolabel(date1, vmat->t1, pdinfo);
	ntolabel(date2, vmat->t2, pdinfo);

	pputs(prn, "{\\rtf1\\par\n\\qc ");
	pprintf(prn, _("Correlation coefficients, using the observations "
		       "%s - %s"), date1, date2);
	pputs(prn, "\\par\n");
	if (vmat->missing) {
	    pprintf(prn, "%s\\par\n", _("(missing values were skipped)"));
	}
	pprintf(prn, _("5%% critical value (two-tailed) = %.4f for n = %d"),
		rhocrit(n, 0.05), n);
	pputs(prn, "\\par\n\\par\n{");
    } else {
	pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n\\par\n{",
		_("Coefficient covariance matrix"));
    }

    for (i=0; i<=blockmax; i++) {
	int pad;

	nf = i * FIELDS;
	li2 = vmat->dim - nf;
	p = (li2 > FIELDS) ? FIELDS : li2;
	if (p == 0) break;

	pad = (vmat->dim > FIELDS)? FIELDS - p : vmat->dim - p;
	rtf_vmat_row(vmat->dim, prn);
	if (pad) rtf_table_pad(pad, prn);

	/* print the varname headings */
	for (j=0; j<p; j++)  {
	    pprintf(prn, "%s\\cell %s", vmat->names[j + nf],
		    (j == p - 1)? "\\cell \\intbl \\row\n" : "");
	}

	/* print rectangular part, if any, of matrix */
	for (j=0; j<nf; j++) {
	    pputs(prn, "\\intbl ");
	    if (pad) {
		rtf_table_pad(pad, prn);
	    }
	    for (k=0; k<p; k++) {
		idx = ijton(j, nf+k, vmat->dim);
		if (vmat->ci == CORR) {
		    rtf_outxx(vmat->vec[idx], prn);
		} else {
		    printf_rtf(vmat->vec[idx], prn, 0);
		}
	    }
	    pprintf(prn, "\\ql %s\\cell \\intbl \\row\n", vmat->names[j]);
	}

	/* print upper triangular part of matrix */
	for (j=0; j<p; ++j) {
	    pputs(prn, "\\intbl ");
	    rtf_table_pad(pad + j, prn);
	    ij2 = nf + j;
	    for (k=j; k<p; k++) {
		idx = ijton(ij2, nf+k, vmat->dim);
		if (vmat->ci == CORR) {
		    rtf_outxx(vmat->vec[idx], prn);
		} else {
		    printf_rtf(vmat->vec[idx], prn, 0);
		}
	    }
	    pprintf(prn, "\\ql %s\\cell \\intbl \\row\n", vmat->names[ij2]);
	}

	if (i < blockmax) {
	    rtf_vmat_blank_row(vmat->dim, pad + p + 1, prn);
	}
    }

    pputs(prn, "}}\n");
}

static void
texprint_vmatrix (const VMatrix *vmat, const DATASET *pdinfo, PRN *prn)
{
    register int i, j;
    int n = vmat->t2 - vmat->t1 + 1;
    int lo, nf, li2, p, k, idx, ij2;
    char vname[2*VNAMELEN];
    int fields = 5;

    lo = vmat->dim;

    if (vmat->ci == CORR) {
	char date1[OBSLEN], date2[OBSLEN];

	ntolabel(date1, vmat->t1, pdinfo);
	ntolabel(date2, vmat->t2, pdinfo);

	pputs(prn, "\\begin{center}\n");
	pprintf(prn, _("Correlation coefficients, using the observations "
		       "%s--%s"), date1, date2);
	pputs(prn, "\\\\\n");
	if (vmat->missing) {
	    pputs(prn, _("(missing values were skipped)"));
	    pputs(prn, "\\\\\n");
	}
	pprintf(prn, _("5\\%% critical value (two-tailed) = %.4f for n = %d"),
		rhocrit(n, 0.05), n);
	pputs(prn, "\\\\\n");
    } else {
	pprintf(prn, "\\begin{center}\n%s\\\\\n",
		_("Coefficient covariance matrix"));
    }

    pputs(prn, "\\vspace{8pt}\n");

    for (i=0; i<=lo/fields; i++) {
	nf = i * fields;
	li2 = lo - nf;
	/* p = number of cols we'll print */
	p = (li2 > fields) ? fields : li2;
	if (p == 0) break;

	pputs(prn, "\\begin{tabular}{");
	for (j=0; j<p; j++) {
	    pputc(prn, 'r');
	}
	pputs(prn, "l}\n");

	/* print the varname headings */
	for (j=0; j<p; j++)  {
	    tex_escape(vname, vmat->names[j + nf]);
	    if (vmat->ci == CORR) {
		pprintf(prn, "%s%s", vname,
			(j == p - 1)? " &\\\\\n" : " & ");
	    } else {
		pprintf(prn, "\\multicolumn{1}{c}{%s}%s", vname,
			(j == p - 1)? " &\\\\\n" : " &\n");
	    }
	}

	/* print rectangular part, if any, of matrix */
	for (j=0; j<nf; j++) {
	    for (k=0; k<p; k++) {
		idx = ijton(j, nf+k, lo);
		if (vmat->ci == CORR) {
		    tex_outxx(vmat->vec[idx], prn);
		} else {
		    tex_matnum(vmat->vec[idx], prn);
		}
	    }
	    tex_escape(vname, vmat->names[j]);
	    pprintf(prn, "%s\\\\\n", vname);
	}

	/* print upper triangular part of matrix */
	for (j=0; j<p; ++j) {
	    ij2 = nf + j;
	    for (k=0; k<j; k++) {
		pputs(prn, " & ");
	    }
	    for (k=j; k<p; k++) {
		idx = ijton(ij2, nf+k, lo);
		if (vmat->ci == CORR) {
		    tex_outxx(vmat->vec[idx], prn);
		} else {
		    tex_matnum(vmat->vec[idx], prn);
		}
	    }
	    tex_escape(vname, vmat->names[ij2]);
	    pprintf(prn, "%s\\\\\n", vname);
	}

	pputs(prn, "\\end{tabular}\n\n");
    }

    pputs(prn, "\\end{center}\n");
}

void special_print_vmatrix (const VMatrix *vmat, const DATASET *pdinfo,
			    PRN *prn)
{
    if (tex_format(prn)) {
	texprint_vmatrix(vmat, pdinfo, prn);
    } else if (rtf_format(prn)) {
	rtfprint_vmatrix(vmat, pdinfo, prn);
    }
}

static int texprint_fcast_stats (const FITRESID *fr,
				 gretlopt opt,
				 PRN *prn)
{
    const char *strs[] = {
	N_("Mean Error"),
	N_("Mean Squared Error"),
	N_("Root Mean Squared Error"),
	N_("Mean Absolute Error"),
	N_("Mean Percentage Error"),
	N_("Mean Absolute Percentage Error"),
	N_("Theil's $U_1$"),
	N_("Bias proportion, $U^M$"),
	N_("Regression proportion, $U^R$"),
	N_("Disturbance proportion, $U^D$")
    };
    const char *U2_str = N_("Theil's $U_2$");
    gretl_matrix *m;
    double x;
    int i, j, t1, t2;
    int n_used = 0;
    int len, err = 0;

    fcast_get_continuous_range(fr, &t1, &t2);

    if (t2 - t1 + 1 <= 0) {
	return E_MISSDATA;
    }

    m = forecast_stats(fr->actual, fr->fitted, t1, t2, &n_used,
		       opt, &err);
    if (err) {
	return err;
    }

    len = gretl_vector_get_length(m);

    pputs(prn, _("Forecast evaluation statistics"));
    pprintf(prn, " (T = %d)", n_used);
    pputs(prn, "\\\\[1ex]\n\n");
    pputs(prn, "\\begin{tabular}{ll}\n");

    j = 0;
    for (i=0; i<len; i++) {
	const char *sj;

	x = gretl_vector_get(m, i);
	if (!isnan(x)) {
	    sj = strs[j];
	    if ((opt & OPT_T) && !strncmp(sj, "Theil", 5)) {
		sj = U2_str;
	    }
	    pprintf(prn, "%s & %s%.5g \\\\\n", _(sj), (x < 0)? "$-$" : "",
		    fabs(x));
	    if (i == 1) {
		pprintf(prn, "%s & %.5g \\\\\n", _(strs[j+1]), sqrt(x));
	    }
	}
	j += (i == 1)? 2 : 1;
    }

    pputs(prn, "\\end{tabular}\n");

    gretl_matrix_free(m);

    return err;
}

static
void tex_fit_resid_head (const FITRESID *fr, const DATASET *pdinfo,
			 PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN];

    ntolabel(date1, fr->t1, pdinfo);
    ntolabel(date2, fr->t2, pdinfo);

    pputs(prn, "\\begin{raggedright}\n");
    pputs(prn, _("Model estimation range:"));
    pprintf(prn, " %s--%s \\\\ \n", date1, date2);

    pprintf(prn, _("Standard error of residuals = %g"), fr->sigma);
    pputs(prn, "\n\\end{raggedright}\n");
}

static
void rtf_fit_resid_head (const FITRESID *fr, const DATASET *pdinfo,
			 PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN];

    ntolabel(date1, fr->t1, pdinfo);
    ntolabel(date2, fr->t2, pdinfo);

    pputs(prn, "{\\rtf1\\par\n\\qc ");
    pputs(prn, _("Model estimation range:"));
    pprintf(prn, " %s - %s\\par\n", date1, date2);

    pputs(prn, "\\qc ");
    pprintf(prn, _("Standard error of residuals = %g"), fr->sigma);
    pputs(prn, "\\par\n\\par\n");
}

static void tex_print_x (double x, int pmax, PRN *prn)
{
    if (x < 0) {
	pputs(prn, "$-$");
    }

    x = fabs(x);

    if (pmax != PMAX_NOT_AVAILABLE) {
	pprintf(prn, "%.*f", pmax, x);
    } else {
	pprintf(prn, "%g", x);
    }

    pputs(prn, " & ");
}

static void texprint_fit_resid (const FITRESID *fr,
				const DATASET *dset,
				PRN *prn)
{
    gretlopt fc_opt = OPT_NONE;
    int t, anyast = 0;
    double xx;
    char vname[2*VNAMELEN];

    tex_fit_resid_head(fr, dset, prn);

    tex_escape(vname, fr->depvar);

    pprintf(prn, "\n\\begin{center}\n"
	    "\\begin{longtable}{rrrrl}\n"
	    " & \n"
	    " \\multicolumn{1}{c}{%s} & \n"
	    "  \\multicolumn{1}{c}{%s} & \n"
	    "   \\multicolumn{1}{c}{%s}\\\\\n",
	    vname, _("fitted"), _("residual"));

    for (t=fr->t1; t<=fr->t2; t++) {
	tex_print_obs_marker(t, dset, prn);
	pputs(prn, " & ");

	if (na(fr->actual[t])) {
	    ;
	} else if (na(fr->fitted[t])) {
	    tex_print_x(fr->actual[t], fr->pmax, prn);
	} else {
	    int ast;

	    xx = fr->actual[t] - fr->fitted[t];
	    ast = (fabs(xx) > 2.5 * fr->sigma);
	    if (ast) anyast = 1;
	    tex_print_x(fr->actual[t], fr->pmax, prn);
	    tex_print_x(fr->fitted[t], fr->pmax, prn);
	    tex_print_x(xx, fr->pmax, prn);
	    if (ast) {
		pputs(prn, " *");
	    }
	}
	pputs(prn, " \\\\\n");
    }

    pputs(prn, "\\end{longtable}\n\n");

    if (anyast) {
	pputs(prn, _("\\textit{Note}: * denotes a residual "
		     "in excess of 2.5 standard errors\n\n"));
    }

    if (dataset_is_time_series(dset)) {
	fc_opt |= OPT_T;
    }
    texprint_fcast_stats(fr, fc_opt, prn);

    pputs(prn, "\\end{center}\n\n");
}

#define FR_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
    "\\cellx800\\cellx2400\\cellx4000\\cellx5600"		\
    "\\cellx6100\n"

static void rtfprint_fit_resid (const FITRESID *fr,
				const DATASET *pdinfo,
				PRN *prn)
{
    double xx;
    int anyast = 0;
    int t;

    rtf_fit_resid_head(fr, pdinfo, prn);

    pputs(prn, "{" FR_ROW "\\intbl ");
    pprintf(prn,
	    " \\qc \\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\ql \\cell"
	    " \\intbl \\row\n",
	    fr->depvar, _("fitted"), _("residual"));

    for (t=fr->t1; t<=fr->t2; t++) {
	rtf_print_obs_marker(t, pdinfo, prn);
	if (na(fr->actual[t])) {
	    pputs(prn, "\\qc \\cell \\qc \\cell \\qc \\cell \\ql \\cell"
		  " \\intbl \\row\n");
	} else if (na(fr->fitted[t])) {
	    printf_rtf(fr->actual[t], prn, 0);
	    pputs(prn, "\\qc \\cell \\qc \\cell \\ql \\cell"
		  " \\intbl \\row\n");
	} else {
	    int ast;

	    xx = fr->actual[t] - fr->fitted[t];
	    ast = (fabs(xx) > 2.5 * fr->sigma);
	    if (ast) {
		anyast = 1;
	    }
	    printf_rtf(fr->actual[t], prn, 0);
	    printf_rtf(fr->fitted[t], prn, 0);
	    printf_rtf(xx, prn, 0);
	    pprintf(prn, "\\ql %s\\cell \\intbl \\row\n",
		    (ast)? "*" : "");
	}
    }

    pputs(prn, "}\n");
    if (anyast) {
	pprintf(prn, "\\par\n\\qc %s \\par\n",
		_("Note: * denotes a residual in excess of 2.5 standard errors"));
    }
    pputs(prn, "}\n");
}

void special_print_fit_resid (const FITRESID *fr,
			      const DATASET *pdinfo,
			      PRN *prn)
{
    if (tex_format(prn)) {
	texprint_fit_resid(fr, pdinfo, prn);
    } else if (rtf_format(prn)) {
	rtfprint_fit_resid(fr, pdinfo, prn);
    }
}

static void texprint_fcast_x (double x, int places, char *str)
{
    if (places != PMAX_NOT_AVAILABLE && !na(x)) {
	tex_rl_float(x, str, places);
    } else {
	tex_rl_double(x, str);
    }
}

static void texprint_fcast_without_errs (const FITRESID *fr,
					 const DATASET *pdinfo,
					 PRN *prn)
{
    char actual[32], fitted[32];
    char vname[2*VNAMELEN];
    char pt = get_local_decpoint();
    int t;

    pputs(prn, "%% The table below needs the \"longtable\" package\n\n");

    pprintf(prn, "\\begin{center}\n"
	    "\\begin{longtable}{%%\n"
	    "r%% col 1: obs\n"
	    "  l%% col 2: varname\n"
	    "    r@{%c}l}%% col 3: fitted\n",
	    pt);

    tex_escape(vname, fr->depvar);

    pprintf(prn, "%s & %s & \\multicolumn{1}{c}{%s} \\\\ [4pt] \n",
	    _("Obs"), vname, _("prediction"));

    for (t=fr->t1; t<=fr->t2; t++) {
	texprint_fcast_x(fr->actual[t], fr->pmax, actual);
	texprint_fcast_x(fr->fitted[t], fr->pmax, fitted);
	tex_print_obs_marker(t, pdinfo, prn);
	pprintf(prn, " & %s & %s \\\\\n",
		actual, fitted);
    }

    pputs(prn, "\\end{longtable}\n\n");
    texprint_fcast_stats(fr, OPT_D, prn);
    pputs(prn, "\\end{center}\n\n");
}

static void texprint_fcast_with_errs (const FITRESID *fr,
				      const DATASET *pdinfo,
				      PRN *prn)
{
    double maxerr, tval = 0;
    double conf = 100 * (1 - fr->alpha);
    int pmax = fr->pmax;
    int errpmax = fr->pmax;
    char actual[32], fitted[32], sderr[32], lo[32], hi[32];
    char vname[2*VNAMELEN];
    gchar *tmp = NULL;
    char pt = get_local_decpoint();
    int t;

    pputs(prn, "\\begin{center}\n");

    if (fr->asymp) {
	tval = normal_critval(fr->alpha / 2);
	pprintf(prn, _("For %g\\%% confidence intervals, $z(%g) = %.2f$\n\n"),
		conf, fr->alpha / 2, tval);
    } else {
	tval = student_critval(fr->df, fr->alpha / 2);
	pprintf(prn, _("For %g\\%% confidence intervals, $t(%d, %g) = %.3f$\n\n"),
		conf, fr->df, fr->alpha / 2, tval);
    }

    pputs(prn, "\\end{center}\n");

    pputs(prn, "%% The table below needs the "
	  "\"longtable\" package\n\n");

    pprintf(prn, "\\begin{center}\n"
	    "\\begin{longtable}{%%\n"
	    "r%% col 1: obs\n"
	    "  r@{%c}l%% col 2: actual\n"
	    "    r@{%c}l%% col 3: fitted\n"
	    "      r@{%c}l%% col 4: std error\n"
	    "        r@{%c}l%% col 5: conf int lo\n"
	    "         r@{%c}l}%% col 5: conf int hi\n",
	    pt, pt, pt, pt, pt);

    tex_escape(vname, fr->depvar);
    tmp = g_strdup_printf(_("%g\\%% interval"), conf);

    pprintf(prn, "%s & \\multicolumn{2}{c}{%s} "
	    " & \\multicolumn{2}{c}{%s}\n"
	    "  & \\multicolumn{2}{c}{%s}\n"
	    "   & \\multicolumn{4}{c}{%s} \\\\[1ex]\n",
	    _("Obs"), vname,
	    _("prediction"), _("std. error"),
	    tmp);
    g_free(tmp);

    if (pmax < 4) {
	errpmax = pmax + 1;
    }

    for (t=fr->t1; t<=fr->t2; t++) {
	double xlo, xhi;

	if (na(fr->sderr[t])) {
	    xlo = xhi = NADBL;
	} else {
	    maxerr = tval * fr->sderr[t];
	    xlo = fr->fitted[t] - maxerr;
	    xhi = fr->fitted[t] + maxerr;
	}
	texprint_fcast_x(fr->actual[t], pmax, actual);
	texprint_fcast_x(fr->fitted[t], pmax, fitted);
	texprint_fcast_x(fr->sderr[t], errpmax, sderr);
	texprint_fcast_x(xlo, pmax, lo);
	texprint_fcast_x(xhi, pmax, hi);
	tex_print_obs_marker(t, pdinfo, prn);
	pprintf(prn, " & %s & %s & %s & %s & %s \\\\\n",
		actual, fitted, sderr, lo, hi);
    }

    pputs(prn, "\\end{longtable}\n\n");
    texprint_fcast_stats(fr, OPT_D, prn);
    pputs(prn, "\\end{center}\n\n");

}

#define FC_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
    "\\cellx800\\cellx2200\\cellx3600\n"

#define FCE_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262"	\
    "\\cellx800\\cellx2200\\cellx3600\\cellx5000"			\
    "\\cellx7800\n"

static void rtfprint_fcast_without_errs (const FITRESID *fr,
					 const DATASET *pdinfo,
					 PRN *prn)
{
    int t;

    pputs(prn, "{\\rtf1\\par\n\n");

    pputs(prn, "{" FC_ROW "\\intbl ");

    pprintf(prn,
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    _("Obs"), fr->depvar, _("prediction"));

    for (t=fr->t1; t<=fr->t2; t++) {
	rtf_print_obs_marker(t, pdinfo, prn);
	printf_rtf(fr->actual[t], prn, 0);
	printf_rtf(fr->fitted[t], prn, 0);
    }

    pputs(prn, "}}\n");
}

static void rtfprint_fcast_with_errs (const FITRESID *fr,
				      const DATASET *pdinfo,
				      PRN *prn)
{
    double maxerr, tval = 0;
    double conf = 100 * (1 - fr->alpha);
    gchar *tmp;
    int d, t;

    pputs(prn, "{\\rtf1\\par\n\\qc ");
    if (fr->asymp) {
	tval = normal_critval(fr->alpha / 2);
	pprintf(prn, _("For %g%% confidence intervals, z(%g) = %.2f"),
		conf, fr->alpha / 2, tval);
    } else {
	tval = student_critval(fr->df, fr->alpha / 2);
	pprintf(prn, _("For %g%% confidence intervals, t(%d, %g) = %.3f"),
		conf, fr->df, fr->alpha / 2, tval);
    }
    pputs(prn, "\\par\n\\par\n");

    tmp = g_strdup_printf(_("%g%% interval"), conf);

    pputs(prn, "{" FCE_ROW "\\intbl ");
    pprintf(prn,
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    _("Obs"), fr->depvar, _("prediction"),
	    _("std. error"),
	    tmp);
    g_free(tmp);

    d = get_gretl_digits();

    for (t=fr->t1; t<=fr->t2; t++) {
	rtf_print_obs_marker(t, pdinfo, prn);
	maxerr = tval * fr->sderr[t];
	printf_rtf(fr->actual[t], prn, 0);
	printf_rtf(fr->fitted[t], prn, 0);
	printf_rtf(fr->sderr[t], prn, 0);
	if (na(fr->sderr[t])) {
	    pputs(prn, "\\qc \\cell \\intbl \\row\n");
	} else {
	    maxerr = tval * fr->sderr[t];
	    pprintf(prn, "\\qc (%#.*g, %#.*g)\\cell \\intbl \\row\n",
		    d, fr->fitted[t] - maxerr,
		    d, fr->fitted[t] + maxerr);
	}
    }

    pputs(prn, "}}\n");
}

void special_print_forecast (const FITRESID *fr,
			     const DATASET *pdinfo,
			     PRN *prn)
{
    if (tex_format(prn)) {
	if (fr->sderr != NULL) {
	    texprint_fcast_with_errs(fr, pdinfo, prn);
	} else {
	    texprint_fcast_without_errs(fr, pdinfo, prn);
	}
    } else if (rtf_format(prn)) {
	if (fr->sderr != NULL) {
	    rtfprint_fcast_with_errs(fr, pdinfo, prn);
	} else {
	    rtfprint_fcast_without_errs(fr, pdinfo, prn);
	}
    }
}

static void
texprint_coeff_interval (const CoeffIntervals *cf, int i, PRN *prn)
{
    int odds = (cf->opt & OPT_O);
    double b = cf->coeff[i];
    double me = cf->maxerr[i];
    char vname[2*VNAMELEN];
    char tmp[32];

    tex_escape(vname, cf->names[i]);
    pprintf(prn, " %s & ", vname);

    if (isnan(b)) {
	pprintf(prn, "\\multicolumn{2}{c}{%s} & ", _("undefined"));
    } else {
	tex_rl_double(odds ? exp(b) : b, tmp);
	pprintf(prn, "%s & ", tmp);
    }

    if (cf->opt & OPT_E) {
	double se = me / cf->t;

	if (isnan(se)) {
	    pprintf(prn, "\\multicolumn{2}{c}{%s} & ", _("undefined"));
	} else {
	    tex_rl_double(odds ? exp(b) * se : se, tmp);
	    pprintf(prn, "%s & ", tmp);
	}
    }

    if (isnan(me)) {
	pprintf(prn, "\\multicolumn{4}{c}{%s}", _("undefined"));
    } else {
	char lo_str[32], hi_str[32];
	double lo = b - me;
	double hi = b + me;

	tex_rl_double(odds ? exp(lo) : lo, lo_str);
	tex_rl_double(odds ? exp(hi) : hi, hi_str);
	pprintf(prn, "%s & %s", lo_str, hi_str);
    }

    pputs(prn, "\\\\\n");
}

static void texprint_confints (const CoeffIntervals *cf, PRN *prn)
{
    const char *hds[] = {
	N_("coefficient"),
	N_("odds ratio")
    };
    const char *hd = hds[0];
    char pt = get_local_decpoint();
    double tail = cf->alpha / 2;
    gchar *cstr;
    int i;

    if (cf->asy) {
	pprintf(prn, "$z(%g) = %.3f$\n\n", tail, cf->t);
    } else {
	pprintf(prn, "$t(%d, %g) = %.3f$\n\n", cf->df, tail, cf->t);
    }

    if (cf->opt & OPT_O) {
	hd = hds[1];
    }

    if (cf->opt & OPT_E) {
	pprintf(prn, "\\begin{center}\n"
		"\\begin{tabular}{rr@{%c}lr@{%c}lr@{%c}lr@{%c}l}\n",
		pt, pt, pt, pt);
    } else {
	pprintf(prn, "\\begin{center}\n"
		"\\begin{tabular}{rr@{%c}lr@{%c}lr@{%c}l}\n",
		pt, pt, pt);
    }

    cstr = g_strdup_printf(_("%g\\%% confidence interval"), 100 * (1 - cf->alpha));

    if (cf->opt & OPT_E) {
	pprintf(prn, " & \\multicolumn{2}{c}{%s}%%\n"
		" & \\multicolumn{2}{c}{%s}%%\n"
		"  & \\multicolumn{4}{c}{%s}\\\\[1ex]\n",
		_(hd), _("std. error"), cstr);
    } else {
	pprintf(prn, " & \\multicolumn{2}{c}{%s}%%\n"
		"  & \\multicolumn{4}{c}{%s}\\\\[1ex]\n",
		_(hd), cstr);
    }

    g_free(cstr);

    for (i=0; i<cf->ncoeff; i++) {
	texprint_coeff_interval(cf, i, prn);
    }

    pputs(prn, "\\end{tabular}\n"
	  "\\end{center}\n");
}

static void
rtfprint_coeff_interval (const CoeffIntervals *cf, int i, PRN *prn)
{
    int odds = (cf->opt & OPT_O);
    double b = cf->coeff[i];
    double me = cf->maxerr[i];
    int d = get_gretl_digits();

    pprintf(prn, "\\qc %s\\cell", cf->names[i]);

    printf_rtf(odds ? exp(b) : b, prn, 0);

    if (cf->opt & OPT_E) {
	double se = me / cf->t;

	fprintf(stderr, "HERE se = %g\n", se);

	printf_rtf(odds ? exp(b) * se : se, prn, 0);
    }

    if (isnan(me)) {
	pprintf(prn, "\\qc %s\\cell ", _("undefined"));
    } else {
	double lo = b - me;
	double hi = b + me;

	pprintf(prn, "\\qc [%#.*g, %#.*g]\\cell ",
		d, odds ? exp(lo) : lo,
		d, odds ? exp(hi) : hi);
    }
    pputs(prn, " \\intbl \\row\n");
}

#define CF_ROW  "\\trowd \\trgaph60\\trleft-30\\trrh262" \
    "\\cellx2400\\cellx4000\\cellx7200\n"

#define CF_ROW2  "\\trowd \\trgaph60\\trleft-30\\trrh262" \
    "\\cellx2400\\cellx4000\\cellx5600\\cellx8800\n"

static void rtfprint_confints (const CoeffIntervals *cf, PRN *prn)
{
    int odds = (cf->opt & OPT_O);
    double tail = cf->alpha / 2;
    gchar *cstr;
    int i;

    if (cf->asy) {
	pprintf(prn, "{\\rtf1\\par\n\\qc z(%g) = %.3f\\par\n\\par\n",
		tail, cf->t);
    } else {
	pprintf(prn, "{\\rtf1\\par\n\\qc t(%d, %g) = %.3f\\par\n\\par\n",
		cf->df, tail, cf->t);
    }

    cstr = g_strdup_printf(_("%g\\%% confidence interval"), 100 * (1 - cf->alpha));

    if (cf->opt & OPT_E) {
	pputs(prn, "{" CF_ROW2 "\\intbl ");
	pprintf(prn,
		" \\qc \\cell"
		" \\qc %s\\cell"
		" \\qc %s\\cell"
		" \\qc %s\\cell"
		" \\intbl \\row\n",
		odds ? _("odds ratio") : _("coefficient"),
		_("std. error"),
		cstr);
    } else {
	pputs(prn, "{" CF_ROW "\\intbl ");
	pprintf(prn,
		" \\qc \\cell"
		" \\qc %s\\cell"
		" \\qc %s\\cell"
		" \\intbl \\row\n",
		odds ? _("odds ratio") :  _("coefficient"),
		cstr);
    }

    g_free(cstr);

    for (i=0; i<cf->ncoeff; i++) {
	rtfprint_coeff_interval(cf, i, prn);
    }

    pputs(prn, "}}\n");
}

void special_print_confints (const CoeffIntervals *cf, PRN *prn)
{
    if (tex_format(prn)) {
	texprint_confints(cf, prn);
    } else if (rtf_format(prn)) {
	rtfprint_confints(cf, prn);
    }
}

int scalars_to_prn (PRN *prn)
{
    GList *slist, *tail;
    char decpoint = get_data_export_decpoint();
    char delim = get_data_export_delimiter();
    const char *sname;
    double sval;
    user_var *u;

    if (delim == ',' && ',' == decpoint) {
	errbox(_("You can't use the same character for "
		 "the column delimiter and the decimal point"));
	return 1;
    }

    tail = slist = user_var_list_for_type(GRETL_TYPE_DOUBLE);

    if (decpoint != ',') {
	gretl_push_c_numeric_locale();
    }

    while (tail) {
	u = tail->data;
	sname = user_var_get_name(u);
	sval = user_var_get_scalar_value(u);
	if (na(sval)) {
	    pprintf(prn, "%s%cNA\n", sname, delim);
	} else {
	    pprintf(prn, "%s%c%.15g\n", sname, delim, sval);
	}
	tail = tail->next;
    }

    g_list_free(slist);

    if (decpoint != ',') {
	gretl_pop_c_numeric_locale();
    }

    return 0;
}

static int data_to_buf_as_rtf (const int *list, PRN *prn)
{
    int err;

    gretl_print_set_format(prn, GRETL_FORMAT_RTF);
    err = print_data_in_columns(list, NULL, dataset, OPT_NONE, prn);
    return err;
}

static int data_to_buf_as_csv (const int *list, gretlopt opt,
			       PRN *prn)
{
    int err;

    gretl_print_set_format(prn, GRETL_FORMAT_CSV);
    err = print_data_in_columns(list, NULL, dataset, opt, prn);
    return err;
}

static int real_series_to_clipboard (const int *list,
				     int format)
{
    PRN *prn = NULL;
    gretlopt opt = OPT_NONE;
    int err = 0;

    if (bufopen(&prn)) {
	return 1;
    }

    if (format == GRETL_FORMAT_XML) {
	err = gretl_write_gdt_to_prn(prn, list, dataset);
    } else {
	if (get_csv_exclude_obs()) {
	    opt = OPT_X;
	}
	err = data_to_buf_as_csv(list, opt, prn);
    }

    if (!err) {
	err = prn_to_clipboard(prn, format);
	if (err) {
	    fprintf(stderr, "prn_to_clipboard: err = %d\n", err);
	}
    } else {
	fprintf(stderr, "series_to_clipboard: err = %d\n", err);
    }

    gretl_print_destroy(prn);

    return err;
}

int selected_series_to_clipboard (void)
{
    int *list = main_window_selection_as_list();
    int err = 0;

    if (list != NULL) {
	const char *opts[] = {
	    N_("Delimited text"),
	    N_("XML (gdt)")
	};
	int resp = radio_dialog(_("Copy to clipboard"),
				_("Format:"),
				opts, 2, 0, 0, mdata->main);

	if (resp == 1) {
	    err = real_series_to_clipboard(list, GRETL_FORMAT_XML);
	} else if (resp == 0) {
	    resp = csv_options_dialog(COPY_CSV, GRETL_OBJ_DSET, NULL);
	    if (!canceled(resp)) {
		err = real_series_to_clipboard(list, GRETL_FORMAT_CSV);
	    }
	}
	free(list);
    }

    return err;
}

/* called from session.c: copy data to clipboard */

int csv_to_clipboard (GtkWidget *parent)
{
    gchar *liststr;
    int cancel, err = 0;

    data_export_selection_wrapper(COPY_CSV);
    liststr = get_selector_storelist();

    if (liststr != NULL) {
	int *list = command_list_from_string(liststr, &err);

	if (list != NULL) {
	    cancel = csv_options_dialog(COPY_CSV, GRETL_OBJ_DSET,
					parent);
	    if (!cancel) {
		err = real_series_to_clipboard(list, 0);
	    }
	    free(list);
	}
	g_free(liststr);
    }

    return err;
}

static void matrix_print_as_csv (const gretl_matrix *m, PRN *prn)
{
    char decpoint = get_data_export_decpoint();
    char delim = get_data_export_delimiter();
    char numstr[48];
    double x;
    int i, j;

    gretl_push_c_numeric_locale();

    for (i=0; i<m->rows; i++) {
	for (j=0; j<m->cols; j++) {
	    x = gretl_matrix_get(m, i, j);
	    sprintf(numstr, "%.*g", DBL_DIG, x);
	    if (decpoint != '.') {
		gretl_charsub(numstr, '.', decpoint);
	    }
	    pputs(prn, numstr);
	    if (j < m->cols - 1) {
		pputc(prn, delim);
	    }
	}
	pputc(prn, '\n');
    }

    gretl_pop_c_numeric_locale();
}

int matrix_to_clipboard_as_csv (const gretl_matrix *m,
				GtkWidget *parent)
{
    if (m != NULL) {
	int resp = csv_options_dialog(COPY_CSV, GRETL_OBJ_MATRIX,
				      parent);
	PRN *prn = NULL;

	if (canceled(resp)) {
	    return 0;
	} else if (bufopen(&prn)) {
	    return 1;
	} else {
	    matrix_print_as_csv(m, prn);
	    prn_to_clipboard(prn, GRETL_FORMAT_CSV);
	    gretl_print_destroy(prn);
	}
    }

    return 0;
}

int scalars_to_clipboard_as_csv (GtkWidget *parent)
{
    int err = 0;

    if (n_user_scalars() == 0) {
	warnbox(_("No scalar variables are currently defined"));
    } else {
	int resp = csv_options_dialog(COPY_CSV, GRETL_OBJ_SCALARS,
				      parent);
	PRN *prn = NULL;

	if (canceled(resp)) {
	    return 0;
	} else if (bufopen(&prn)) {
	    return 1;
	} else {
	    err = scalars_to_prn(prn);
	    if (!err) {
		prn_to_clipboard(prn, GRETL_FORMAT_CSV);
	    }
	    gretl_print_destroy(prn);
	}
    }

    return err;
}

#include "series_view.h"
#include "fileselect.h"

/* callback from "series view" window, for use when
   the delimited or RTF options are chosen
*/

int copy_vars_formatted (windata_t *vwin, int fmt, int action)
{
    char save_delim = get_data_export_delimiter();
    int *list = series_view_get_list(vwin);
    PRN *prn = NULL;
    int i, err = 0;

    if (list != NULL) {
	for (i=1; i<=list[0]; i++) {
	    if (list[i] >= dataset->v) {
		gui_errmsg(E_DATA);
		return E_DATA;
	    }
	}

	if (fmt == GRETL_FORMAT_CSV) {
	    set_data_export_delimiter(',');
	} else if (fmt == GRETL_FORMAT_TAB) {
	    fmt = GRETL_FORMAT_CSV;
	    set_data_export_delimiter('\t');
	}

	if (series_view_is_sorted(vwin)) {
	    prn = vwin_print_sorted_with_format(vwin, fmt);
	    if (prn == NULL) {
		err = 1;
	    }
	} else {
	    err = bufopen(&prn);
	    if (!err) {
		if (fmt == GRETL_FORMAT_RTF) {
		    err = data_to_buf_as_rtf(list, prn);
		} else {
		    err = data_to_buf_as_csv(list, OPT_NONE, prn);
		}
	    }
	}

	if (!err) {
	    if (action == W_COPY) {
		err = prn_to_clipboard(prn, fmt);
	    } else if (fmt == GRETL_FORMAT_RTF) {
		file_selector(SAVE_RTF, FSEL_DATA_PRN, prn);
	    } else {
		file_selector(EXPORT_CSV, FSEL_DATA_PRN, prn);
	    }
	}

	gretl_print_destroy(prn);
	free(list);
    }

    set_data_export_delimiter(save_delim);

    return err;
}

#ifdef G_OS_WIN32

static int get_latex_path (char *latex_path)
{
    int ret;
    char *p;

    ret = SearchPath(NULL, latex, NULL, MAXLEN, latex_path, &p);

    return (ret == 0);
}

#else

static int spawn_latex (char *texsrc)
{
    GError *error = NULL;
    gchar *errout = NULL, *sout = NULL;
    gchar *argv[] = {
	latex,
	"\\batchmode",
	"\\input",
	texsrc,
	NULL
    };
    int ok, status;
    int ret = LATEX_OK;

    ok = g_spawn_sync (gretl_dotdir(), /* working dir */
		       argv,
		       NULL,    /* envp */
		       G_SPAWN_SEARCH_PATH,
		       NULL,    /* child_setup */
		       NULL,    /* user_data */
		       &sout,   /* standard output */
		       &errout, /* standard error */
		       &status, /* exit status */
		       &error);

    if (!ok) {
	errbox(error->message);
	g_error_free(error);
	ret = LATEX_EXEC_FAILED;
    } else if (status != 0) {
	if (errout && *errout) {
	    errbox(errout);
	} else {
	    gchar *errmsg;

	    errmsg = g_strdup_printf("%s\n%s",
				     _("Failed to process TeX file"),
				     sout);
	    errbox(errmsg);
	    g_free(errmsg);
	}
	ret = LATEX_ERROR;
    } else if (errout && *errout) {
	fputs("spawn_latex: found stuff on stderr:\n", stderr);
	fputs(errout, stderr);
    }

    /* change above, 2008-08-22: before we flagged a LATEX_ERROR
       if we saw anything on standard error, regardless of the
       exit status
    */

    g_free(errout);
    g_free(sout);

    return ret;
}

#endif /* !G_OS_WIN32 */

int latex_compile (char *texshort)
{
#ifdef G_OS_WIN32
    static char latex_path[MAXLEN];
    gchar *tmp = NULL;
#endif
    int err = LATEX_OK;

#ifdef G_OS_WIN32
    if (*latex_path == '\0' && get_latex_path(latex_path)) {
	win_show_last_error();
	return LATEX_EXEC_FAILED;
    }

    tmp = g_strdup_printf("\"%s\" \\batchmode \\input %s", latex_path, texshort);
    if (win_run_sync(tmp, gretl_dotdir())) {
	err =LATEX_EXEC_FAILED;
    }
    g_free(tmp);
#else
    err = spawn_latex(texshort);
#endif /* G_OS_WIN32 */

    return err;
}

static int check_for_rerun (const char *texbase)
{
    char logfile[MAXLEN];
    char lline[512];
    FILE *fp;
    int ret = 0;

    gretl_path_compose(logfile, MAXLEN, texbase, ".log");
    fp = gretl_fopen(logfile, "r");

    if (fp != NULL) {
	while (fgets(lline, sizeof lline, fp)) {
	    if (strstr(lline, "Rerun LaTeX")) {
		ret = 1;
		break;
	    }
	}
	fclose(fp);
    }

    return ret;
}

static void view_or_save_latex (PRN *bprn, const char *fname, int saveit)
{
    char texfile[MAXLEN], texbase[MAXLEN], tmp[MAXLEN];
    int dot, err = LATEX_OK;
    char *texshort = NULL;
    const char *buf;
    PRN *fprn;

    *texfile = 0;

    if (fname != NULL) {
	strcpy(texfile, fname);
    } else {
	sprintf(texfile, "%swindow.tex", gretl_dotdir());
    }

    /* ensure we don't get stale output */
    remove(texfile);

    fprn = gretl_print_new_with_filename(texfile, &err);
    if (err) {
	gui_errmsg(err);
	return;
    }

    gretl_tex_preamble(fprn, prn_format(bprn));
    buf = gretl_print_get_buffer(bprn);
    pputs(fprn, buf);
    pputs(fprn, "\n\\end{document}\n");

    gretl_print_destroy(fprn);

    if (saveit) {
	return;
    }

    dot = gretl_dotpos(texfile);
    *texbase = 0;
    strncat(texbase, texfile, dot);

    texshort = strrslash(texbase);
    if (texshort == NULL) {
	errbox(_("Failed to process TeX file"));
	return;
    }
    texshort++;

    err = latex_compile(texshort);

    /* now maybe re-run latex (e.g. for longtable) */
    if (err == LATEX_OK) {
	if (check_for_rerun(texbase)) {
	    err = latex_compile(texshort);
	}
    }

    if (err == LATEX_OK) {
#if defined(G_OS_WIN32)
	sprintf(tmp, "%s.pdf", texbase);
	win32_open_file(tmp);
#elif defined(__APPLE__)
	sprintf(tmp, "%s.pdf", texbase);
	if (osx_open_file(tmp)) {
	    file_read_errbox(tmp);
	}
#else
	gretl_path_compose(tmp, MAXLEN, texbase, ".pdf");
	gretl_fork("viewpdf", tmp, NULL);
#endif
    }

    gretl_path_compose(tmp, MAXLEN, texbase, ".log");
    if (err == LATEX_ERROR) {
	view_file(tmp, 0, TMP_FILE, 78, 350, VIEW_FILE);
    } else {
	fprintf(stderr, "wrote '%s'\n", texfile);
	/* gretl_remove(texfile); */
	gretl_remove(tmp);
    }

    gretl_path_compose(tmp, MAXLEN, texbase, ".aux");
    gretl_remove(tmp);
}

void view_latex (PRN *prn)
{
    view_or_save_latex(prn, NULL, 0);
}

void save_latex (PRN *prn, const char *fname)
{
    if (prn != NULL) {
	view_or_save_latex(prn, fname, 1);
    } else {
	save_graph_page(fname);
    }
}
