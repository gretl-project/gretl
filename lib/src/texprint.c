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

    sprintf(number, "%#.*g", GRETL_DIGITS, x);

    if (!tab) {
	if (x < 0.) pprintf(prn, "$-$%s", number + 1);
	else pprintf(prn, "%s", number);
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

void tex_dcolumn_double (double xx, char *numstr)
{
    double a = fabs(xx);

    sprintf(numstr, "%#.*g", GRETL_DIGITS, xx);

    if (a >= UPPER_F_LIMIT || a < LOWER_F_LIMIT) {
	int expon;
	char *p, exponstr[8];

	p = strchr(numstr, 'e');
	expon = atoi(p + 2);
	strcpy(p, "\\mbox{e");
	sprintf(exponstr, "%c%02d}", (xx > 10)? '+' : '-', expon);
	strcat(numstr, exponstr);
    }
}

static void tex_make_cname (const char *orig, char *cname)
{
    char *p;
    unsigned char c;

    if (orig == NULL || strlen(orig) == 0) return;

    p = strrchr(orig, '_');
    if (p == NULL) {
	tex_escape(cname, orig);
	return;
    }

    c = (unsigned char) *(p + 1);

    if (isdigit(c)) {
	int lag = atoi(++p);

	sprintf(cname, "$u_{t-%d}^2$", lag);
    } else {
	tex_escape(cname, orig);
    }
}

int tex_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
		     int c, PRN *prn)
{
    char tmp[16], coeff[32], sderr[32];
    double t_ratio = pmod->coeff[c-1] / pmod->sderr[c-1];

    *tmp = 0;
    if (pmod->aux == AUX_ARCH) {
	tex_make_cname(pdinfo->varname[pmod->list[c]], tmp);
    } else {
	tex_escape(tmp, pdinfo->varname[pmod->list[c]]);
    }
	
    tex_dcolumn_double(pmod->coeff[c-1], coeff);
    tex_dcolumn_double(pmod->sderr[c-1], sderr);

    if (pmod->ci != LOGIT && pmod->ci != PROBIT) {
	pprintf(prn, "%s &\n"
		"  %s &\n"
		"    %s &\n"
		"      %.4f &\n"
		"        %.4f \\\\\n",  
		tmp,
		coeff,
		sderr,
		t_ratio,
		tprob(t_ratio, pmod->dfd));	
    } else { /* LOGIT, PROBIT */
	char slope[32];

	if (pmod->list[c]) {
	    tex_dcolumn_double(pmod->slope[c-1], slope);
	}
	pprintf(prn, "%s &\n"
		"  %s &\n"
		"    %s &\n"
		"      %.4f &\n"
		"        %s \\\\\n",  
		tmp,
		coeff,
		sderr,
		t_ratio,
		(pmod->list[c])? slope : "");
    }

    return 0;
}

/* ......................................................... */

static int make_texfile (const PATHS *ppaths, int model_count,
			 int equation, char *texfile, PRN *prn)
{
    prn->buf = NULL;

    sprintf(texfile, "%s%s_%d.tex", ppaths->userdir,
	    (equation)? "equation" : "model", model_count);

    prn->fp = fopen(texfile, "w");
    if (prn->fp == NULL) return 1;
    else return 0;
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

    int i, start, constneg = 0, ncoeff = pmod->list[0];
    double tstat, const_tstat = 0, const_coeff = 0;
    char tmp[16];

    if (standalone) {
	pprintf(prn, "\\documentclass[11pt]{article}\n\\begin{document}\n"
		"\\thispagestyle{empty}\n\n");
    }
    pprintf(prn, "\\begin{center}\n");

    if (pmod->ifc) {
	const_coeff = pmod->coeff[pmod->list[0]-1];
	const_tstat = pmod->coeff[pmod->list[0]-1]
	    / pmod->sderr[pmod->list[0]-1];
	if (const_coeff < 0.) constneg = 1;
	ncoeff--;
    }

    /* tabular header */
    pprintf(prn, "{\\setlength{\\tabcolsep}{.5ex}\n"
	    "\\renewcommand{\\arraystretch}{1}\n"
	    "\\begin{tabular}{rc"
	    "%s", (pmod->ifc)? "c" : "c@{\\,}l");
    start = (pmod->ifc)? 1 : 2;
    for (i=start; i<ncoeff; i++) pprintf(prn, "cc@{\\,}l");
    pprintf(prn, "}\n");

    /* dependent variable */
    tmp[0] = '\0';
    tex_escape(tmp, pdinfo->varname[pmod->list[1]]);
    pprintf(prn, "$\\widehat{\\rm %s}$ & = &\n", tmp);

    start++;
    /* coefficients times indep vars */
    if (pmod->ifc) tex_print_float(const_coeff, 0, prn);
    else {
	tex_escape(tmp, pdinfo->varname[pmod->list[2]]);
	tex_print_float(pmod->coeff[1], 0, prn);
	pprintf(prn, " & %s ", tmp);
    }
    for (i=start; i<=ncoeff; i++) {
	tex_print_float(pmod->coeff[i-1], 1, prn);
	tex_escape(tmp, pdinfo->varname[pmod->list[i]]);
	pprintf(prn, " & %s ", tmp);
    }
    pprintf(prn, "\\\\\n");

    /* t-stats in row beneath */
    if (pmod->ifc) {
	pprintf(prn, "& ");
	pprintf(prn, "& {\\small $(%.3f)$} ", const_tstat);
    } 
    for (i=2; i<=ncoeff; i++) {
        tstat = pmod->coeff[i-1]/pmod->sderr[i-1];
	if (i == 2) pprintf(prn, "& & \\small{$(%.3f)$} ", tstat);
	else pprintf(prn, "& & & \\small{$(%.3f)$} ", tstat);
    }
    pprintf(prn, "\n\\end{tabular}}\n\n");

    /* additional info (R^2 etc) */
    pprintf(prn, "\\vspace{.8ex}\n");

    if (na(pmod->fstt)) { /* LAD model */
	pprintf(prn, "$T = %d,\\, \\sum |\\hat{u}_t| = %g$\n",
		pmod->nobs, pmod->rho);
    } else {
	pprintf(prn, "$T = %d,\\, \\bar{R}^2 = %.3f,\\, F(%d,%d) = %g,\\, "
		"\\hat{\\sigma} = %g$\n",
		pmod->nobs, pmod->adjrsq, pmod->dfn, 
		pmod->dfd, pmod->fstt, pmod->sigma);
    }

    pprintf(prn, "\n(%s)\n\\end{center}\n", 
	    I_("$t$-statistics in parentheses"));

    if (standalone) 
	pprintf(prn, "\n\\end{document}\n");

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
