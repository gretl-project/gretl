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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* texprint.c for gretl - LaTeX output from modeling */

#include "libgretl.h"

/* ......................................................... */

static void tex_print_float (double x, int tab, PRN *prn)
     /* prints a floating point number as a TeX math string.
	if tab != 0, print the sign in front, separated by a
	tab symbol (for equation-style regression printout).
     */
{
    char number[16];

    x = screen_zero(x);

    sprintf(number, "%#.*g", GRETL_DIGITS, x);

    if (!tab) {
	if (x < 0.) pprintf(prn, "$-$%s", number + 1);
	else pputs(prn, number);
    } else {
	if (x < 0.) pprintf(prn, "& $-$ & %s", number + 1);
	else pprintf(prn, "& $+$ & %s", number);
    }
}

/**
 * tex_escape:
 * @targ: target string (must be pre-allocated)
 * @src: source string.
 *
 * Copies from @src to @targ, escaping any characters in @src that are 
 * special to TeX (by inserting a leading backslash).
 * 
 * Returns: the transformed copy of the string
 */

char *tex_escape (char *targ, const char *src)
{
    while (*src) {
	if (*src == '$' || *src == '&' || *src == '_' || 
	    *src == '%' || *src == '#')
	    *targ++ = '\\';
	*targ++ = *src++;
    }
    *targ = '\0';
    return targ;
}

#define UPPER_F_LIMIT (pow(10, GRETL_DIGITS))
#define LOWER_F_LIMIT (pow(10, -4))

static void cut_extra_zero (char *numstr)
{
    int s = strspn(numstr, "-.,0");
    int p = (s == 0 && (strchr(numstr, '.') || strchr(numstr, ',')));

    numstr[s + p + GRETL_DIGITS] = '\0';
}

void tex_dcolumn_double (double xx, char *numstr)
{
    double a;

    xx = screen_zero(xx);
    a = fabs(xx);

    sprintf(numstr, "%#.*g", GRETL_DIGITS, xx);

    if (a != 0.0 && (a >= UPPER_F_LIMIT || a < LOWER_F_LIMIT)) {
	int expon;
	char *p, exponstr[8];

	p = strchr(numstr, 'e');
	expon = atoi(p + 2);
	strcpy(p, "\\mbox{e");
	sprintf(exponstr, "%s%02d}", (xx > 10)? "+" : "-", expon);
	strcat(numstr, exponstr);
    } else {
	cut_extra_zero(numstr);
    }
}

static void tex_make_cname (char *cname, const char *src)
{
    char *p;
    unsigned char c;

    if (src == NULL || strlen(src) == 0) return;

    p = strrchr(src, '_');
    if (p == NULL) {
	tex_escape(cname, src);
	return;
    }

    c = (unsigned char) *(p + 1);

    if (isdigit(c)) {
	int lag = atoi(++p);

	sprintf(cname, "$u_{t-%d}^2$", lag);
    } else {
	tex_escape(cname, src);
    }
}

static int tex_greek_param (char *targ, const char *src)
{
    *targ = 0;

    if (!strcmp(src, "alpha")) {
	strcpy(targ, "$\\alpha$");
    }
    else if (!strcmp(src, "beta")) {
	strcpy(targ, "$\\beta$");
    }
    else if (!strcmp(src, "gamma")) {
	strcpy(targ, "$\\gamma$");
    }
    else if (!strcmp(src, "delta")) {
	strcpy(targ, "$\\delta$");
    }
    else if (!strcmp(src, "pi")) {
	strcpy(targ, "$\\pi$");
    }
    else if (!strcmp(src, "lambda")) {
	strcpy(targ, "$\\lambda$");
    }

    return (*targ != 0);
}

int tex_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
		     int c, PRN *prn)
{
    char tmp[16], coeff[64], sderr[64], tratio[64], pval[64];

    if (isnan(pmod->coeff[c-2]) || na(pmod->coeff[c-2])) {
	sprintf(coeff, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
    } else {
	tex_dcolumn_double(pmod->coeff[c-2], coeff);
    }

    if (isnan(pmod->sderr[c-2]) || na(pmod->sderr[c-2])) {
	sprintf(sderr, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
	sprintf(tratio, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
	sprintf(pval, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
    } else {
	tex_dcolumn_double(pmod->sderr[c-2], sderr);
	sprintf(tratio, "%.4f", pmod->coeff[c-2] / pmod->sderr[c-2]);
	sprintf(pval, "%.4f", tprob(pmod->coeff[c-2] / pmod->sderr[c-2], 
				    pmod->dfd));
    }    

    *tmp = 0;
    if (pmod->aux == AUX_ARCH) {
	tex_make_cname(tmp, pdinfo->varname[pmod->list[c]]);
    } else if (pmod->ci == NLS) {
	if (!tex_greek_param(tmp, pmod->params[c-1])) {
	    tex_escape(tmp, pmod->params[c-1]);
	}
    } else {
	tex_escape(tmp, pdinfo->varname[pmod->list[c]]);
    }
	
    if (pmod->ci != LOGIT && pmod->ci != PROBIT) {
	pprintf(prn, "%s &\n"
		"  %s &\n"
		"    %s &\n"
		"      %s &\n"
		"        %s \\\\\n",  
		tmp,
		coeff,
		sderr,
		tratio,
		pval);	
    } else { /* LOGIT, PROBIT */
	char slope[32];

	if (pmod->list[c]) {
	    tex_dcolumn_double(pmod->slope[c-2], slope);
	}
	pprintf(prn, "%s &\n"
		"  %s &\n"
		"    %s &\n"
		"      %s &\n"
		"        %s \\\\\n",  
		tmp,
		coeff,
		sderr,
		tratio,
		(pmod->list[c])? slope : "");
    }

    return 0;
}

/* ......................................................... */

static int make_texfile (const PATHS *ppaths, int model_count,
			 int equation, char *texfile, PRN *prn)
{
    FILE *fp;

    if (*texfile == 0) {
	sprintf(texfile, "%s%s_%d.tex", ppaths->userdir,
		(equation)? "equation" : "model", model_count);
    }

    fp = fopen(texfile, "w");
    if (fp == NULL) return 1;

    gretl_print_attach_file(prn, fp);
    return 0;
}

/**
 * tex_print_equation:
 * @pmod:  pointer to gretl MODEL struct.
 * @pdinfo:  information regarding the data set.
 * @standalone: indicator variable.
 * @prn: gretl printing struct.
 *
 * Prints a gretl model in the form of a LaTeX equation, either as
 * a stand-alone document or as a fragment of LaTeX source for
 * insertion into a document.
 * 
 * Returns: 0 on successful completion.
 */

int tex_print_equation (const MODEL *pmod, const DATAINFO *pdinfo, 
			int standalone, PRN *prn)
{

    int i, start, ncoeff = pmod->ncoeff;
    double tstat, const_tstat = 0, const_coeff = 0;
    char tmp[16];

    if (standalone) {
	pputs(prn, "\\documentclass[11pt]{article}\n");

#ifdef ENABLE_NLS
	pputs(prn, "\\usepackage[latin1]{inputenc}\n\n");
#endif

	pputs(prn, "\\begin{document}\n\n"
		"\\thispagestyle{empty}\n\n");
    }

    pputs(prn, "\\begin{center}\n");

    if (pmod->ifc) {
	const_coeff = pmod->coeff[0];
	const_tstat = pmod->coeff[0] / pmod->sderr[0];
    }

    /* tabular header */
    pprintf(prn, "{\\setlength{\\tabcolsep}{.5ex}\n"
	    "\\renewcommand{\\arraystretch}{1}\n"
	    "\\begin{tabular}{rc"
	    "%s", (pmod->ifc)? "c" : "c@{\\,}l");
    start = (pmod->ifc)? 1 : 0;
    for (i=start; i<ncoeff; i++) pputs(prn, "cc@{\\,}l");
    pputs(prn, "}\n");

    /* dependent variable */
    *tmp = '\0';
    tex_escape(tmp, pdinfo->varname[pmod->list[1]]);
    pprintf(prn, "$\\widehat{\\rm %s}$ & = &\n", tmp);

    /* coefficients times indep vars */
    if (pmod->ifc) tex_print_float(const_coeff, 0, prn);
    start = (pmod->ifc)? 1 : 0;
    for (i=start; i<ncoeff; i++) {
	tex_print_float(pmod->coeff[i], (i > 0), prn);
	tex_escape(tmp, pdinfo->varname[pmod->list[i+2]]);
	pprintf(prn, " & %s ", tmp);
    }
    pputs(prn, "\\\\\n");

    /* t-stats in row beneath */
    if (pmod->ifc) {
	pprintf(prn, " & & {\\small $(%.3f)$} ", const_tstat);
    } 
    start = (pmod->ifc)? 1 : 0;
    for (i=start; i<ncoeff; i++) {
        tstat = pmod->coeff[i] / pmod->sderr[i];
	if (i == start) pprintf(prn, "& & \\small{$(%.3f)$} ", tstat);
	else pprintf(prn, "& & & \\small{$(%.3f)$} ", tstat);
    }
    pputs(prn, "\n\\end{tabular}}\n\n");

    /* additional info (R^2 etc) */
    pputs(prn, "\\vspace{.8ex}\n");

    if (pmod->ci == LAD) { 
	pprintf(prn, "$T = %d,\\, \\sum |\\hat{u}_t| = %g$\n",
		pmod->nobs, pmod->rho);
    } else {
	pprintf(prn, "$T$ = %d, $\\, \\bar{R}^2$ = %.3f, ",
		pmod->nobs, pmod->adjrsq);
	if (!na(pmod->fstt)) {
	    pprintf(prn, "$\\, F(%d,%d)$ = %.5g, ", 
		    pmod->dfn, pmod->dfd, pmod->fstt);
	}
	pprintf(prn, "$\\, \\hat{\\sigma}$ = %.4g", pmod->sigma);
	if (!floateq(pmod->rho_in, 0.0)) {
	    double r = pmod->rho_in;
	    char rstr[16];

	    if (r < 0.0) sprintf(rstr, "$-$%.4g", fabs(r));
	    else sprintf(rstr, "%.4g", r);
	    pprintf(prn, ", $\\, \\rho$ = %s", rstr);
	}
	pputs(prn, "\n");
    }

    pprintf(prn, "\n(%s)\n\\end{center}\n", 
	    I_("$t$-statistics in parentheses"));

    if (standalone) {
	pputs(prn, "\n\\end{document}\n");
    }

    return 0;
}

/**
 * tex_print_model:
 * @pmod:  pointer to gretl MODEL struct.
 * @pdinfo:  information regarding the data set.
 * @standalone: indicator variable.
 * @prn: gretl printing struct.
 *
 * Prints a gretl model in the form of a LaTeX table, either as
 * a stand-alone document or as a fragment of LaTeX source for
 * insertion into a document.
 * 
 * Returns: 0 on successful completion.
 */

int tex_print_model (const MODEL *pmod, const DATAINFO *pdinfo, 
		     int standalone, PRN *prn)
{
    if (standalone) prn->format = GRETL_PRINT_FORMAT_TEX_DOC;
    else prn->format = GRETL_PRINT_FORMAT_TEX;
    
    return printmodel (pmod, pdinfo, prn);
}

/**
 * tabprint:
 * @pmod: pointer to gretl MODEL struct.
 * @pdinfo: information regarding the data set.
 * @ppaths: struct containing information on paths.
 * @texfile: name of file to save.
 * @model_count: count of models estimated so far.
 * @oflag: option: standalone or not.
 *
 * Prints to file a gretl model in the form of a LaTeX table, either as
 * a stand-alone document or as a fragment of LaTeX source for
 * insertion into a document.
 * 
 * Returns: 0 on successful completion, 1 on error.
 */

int tabprint (const MODEL *pmod, const DATAINFO *pdinfo,
	      const PATHS *ppaths, char *texfile,
	      int model_count, int oflag)
{
    PRN prn;

    if (make_texfile(ppaths, model_count, 0, texfile, &prn))
	return 1;

    tex_print_model(pmod, pdinfo, oflag, &prn);
    if (prn.fp != NULL) fclose(prn.fp);
    return 0;
}

/**
 * eqnprint:
 * @pmod: pointer to gretl MODEL struct.
 * @pdinfo: information regarding the data set.
 * @ppaths: struct containing information on paths.
 * @texfile: name of file to save.
 * @model_count: count of models estimated so far.
 * @oflag: option: standalone or not.
 *
 * Prints to file a gretl model in the form of a LaTeX equation, either as
 * a stand-alone document or as a fragment of LaTeX source for
 * insertion into a document.
 * 
 * Returns: 0 on successful completion, 1 on error.
 */

int eqnprint (const MODEL *pmod, const DATAINFO *pdinfo,
	      const PATHS *ppaths, char *texfile,
	      int model_count, int oflag)
{
    PRN prn;

    if (make_texfile(ppaths, model_count, 1, texfile, &prn))
	return 1;

    tex_print_equation(pmod, pdinfo, oflag, &prn);
    if (prn.fp != NULL) fclose(prn.fp);
    return 0;
}
