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

static void tex_modify_exponent (char *numstr)
{
    char *p = strchr(numstr, 'e');

    if (p != NULL) {
	int expon = atoi(p + 1);

	sprintf(p, "\\times 10^{%d}", expon);
    }
}

/* ......................................................... */

static void tex_print_float (double x, PRN *prn)
     /* prints a floating point number as a TeX math string.
	if tab != 0, print the sign in front, separated by a
	tab symbol (for equation-style regression printout).
     */
{
    char number[48];

    x = screen_zero(x);

    sprintf(number, "%#.*g", GRETL_DIGITS, x);
    tex_modify_exponent(number);
    pprintf(prn, "%s", (x < 0.0)? (number + 1) : number);
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
    char *p = targ;

    while (*src) {
	if (*src == '$' || *src == '&' || *src == '_' || 
	    *src == '%' || *src == '#')
	    *targ++ = '\\';
	*targ++ = *src++;
    }
    *targ = '\0';

    return p;
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

static void tex_arma_coeff_name (char *targ, const char *src,
				 int inmath)
{
    char vname[VNAMELEN], vnesc[16], texname[24];
    int lag;

    if (sscanf(src, "%[^(](-%d)", vname, &lag) == 2) {
	if (!strcmp(vname, "e")) {
	    if (!inmath) {
		strcpy(texname, "$\\varepsilon$");
	    } else {
		strcpy(texname, "\\varepsilon");
	    }
	} else if (!strcmp(vname, "y")) {
	    strcpy(texname, "y");
	} else {
	    tex_escape(vnesc, vname);
	    if (!inmath) {
		strcpy(texname, vnesc);
	    } else {
		sprintf(texname, "\\mbox{%s}", vnesc);
	    }
	}
	if (!inmath) {
	    sprintf(targ, "%s$_{t-%d}$", texname, lag);
	} else {
	    sprintf(targ, "%s_{t-%d}", texname, lag);
	}
    } else {
	tex_escape(vnesc, src);
	strcpy(targ, vnesc);
    }
}

int tex_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
		     int c, PRN *prn)
{
    char tmp[24], coeff[64], sderr[64], tratio[64], pval[64];
    int v = c - 2;

    if (isnan(pmod->coeff[v]) || na(pmod->coeff[v])) {
	sprintf(coeff, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
    } else {
	tex_dcolumn_double(pmod->coeff[v], coeff);
    }

    if (isnan(pmod->sderr[v]) || na(pmod->sderr[v])) {
	sprintf(sderr, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
	sprintf(tratio, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
	sprintf(pval, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
    } else {
	tex_dcolumn_double(pmod->sderr[v], sderr);
	sprintf(tratio, "%.4f", pmod->coeff[v] / pmod->sderr[v]);
	sprintf(pval, "%.4f", tprob(pmod->coeff[v] / pmod->sderr[v], 
				    pmod->dfd));
    }    

    *tmp = 0;
    if (pmod->aux == AUX_ARCH) {
	tex_make_cname(tmp, pdinfo->varname[pmod->list[c]]);
    } else if (pmod->ci == NLS) {
	if (!tex_greek_param(tmp, pmod->params[c-1])) {
	    tex_escape(tmp, pmod->params[c-1]);
	}
    } else if (pmod->ci == ARMA) {
	tex_arma_coeff_name(tmp, pmod->params[c-1], 0);
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
	double *slopes = gretl_model_get_data(pmod, "slopes");
	char slope[32];

	if (pmod->list[c]) {
	    tex_dcolumn_double(slopes[v], slope);
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
    double tstat;
    char tmp[48];
    int i;

    if (standalone) {
	pputs(prn, "\\documentclass[11pt]{article}\n");

#ifdef ENABLE_NLS
	pputs(prn, "\\usepackage[latin1]{inputenc}\n\n");
#endif
	pputs(prn, "\\usepackage{amsmath}\n\n");

	pputs(prn, "\\begin{document}\n\n"
		"\\thispagestyle{empty}\n\n");
    }

    /* initial setup */
    pputs(prn, "\\begin{gather}\n");

    /* dependent variable */
    *tmp = '\0';
    if (pmod->ci == ARMA) {
	tex_escape(tmp, pdinfo->varname[pmod->list[4]]);
    } else {
	tex_escape(tmp, pdinfo->varname[pmod->list[1]]);
    }
    pprintf(prn, "\\widehat{\\rm %s} = \n", tmp);

    /* coefficients times indep vars */
    for (i=0; i<pmod->ncoeff; i++) {
	tstat = pmod->coeff[i] / pmod->sderr[i];
	pprintf(prn, "%s\\underset{(%.3f)}{", 
		(pmod->coeff[i] < 0.0)? "-" :
		(i > 0)? "+" : "", tstat);
	tex_print_float(pmod->coeff[i], prn);
	pputc(prn, '}');
	if (i > 0 || pmod->ifc == 0) {
	    pputs(prn, "\\,");
	    if (pmod->ci == ARMA) {
		tex_arma_coeff_name(tmp, pmod->params[i+1], 1);
		pprintf(prn, "%s", tmp);
	    } else {
		tex_escape(tmp, pdinfo->varname[pmod->list[i+2]]);
		pprintf(prn, "\\mbox{%s}", tmp);
	    }
	}
	pputc(prn, '\n');
    }
    pputs(prn, " \\notag \\\\\n");

    /* additional info (R^2 etc) */
    if (pmod->ci == LAD) { 
	sprintf(tmp, "%g", pmod->rho);
	tex_modify_exponent(tmp);
	pprintf(prn, "T = %d \\quad \\sum |\\hat{u}_t| = %s",
		pmod->nobs, tmp);
    } else {
	pprintf(prn, "T = %d \\quad \\bar{R}^2 = %.4f ",
		pmod->nobs, pmod->adjrsq);
	if (!na(pmod->fstt)) {
	    sprintf(tmp, "%.5g", pmod->fstt);
	    tex_modify_exponent(tmp);
	    pprintf(prn, "\\quad F(%d,%d) = %s ", 
		    pmod->dfn, pmod->dfd, tmp);
	}
	sprintf(tmp, "%.5g", pmod->sigma);
	tex_modify_exponent(tmp);
	pprintf(prn, "\\quad \\hat{\\sigma} = %s", tmp);
	if (!na(gretl_model_get_double(pmod, "rho_in"))) {
	    double r = gretl_model_get_double(pmod, "rho_in");

	    sprintf(tmp, "%.5g", r);
	    tex_modify_exponent(tmp);
	    pprintf(prn, " \\quad \\rho = %s", tmp);
	}
    }

    pputs(prn, "\\notag \\\\\n");
    pprintf(prn, "\\centerline{(%s)} \\notag\n",
	  I_("$t$-statistics in parentheses"));
    pputs(prn, "\\end{gather}\n");

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
	      int model_count, unsigned long oflag)
{
    PRN prn;

    if (make_texfile(ppaths, model_count, 0, texfile, &prn))
	return 1;

    tex_print_model(pmod, pdinfo, (int) oflag, &prn);
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
	      int model_count, unsigned long oflag)
{
    PRN prn;

    if (make_texfile(ppaths, model_count, 1, texfile, &prn))
	return 1;

    tex_print_equation(pmod, pdinfo, (int) oflag, &prn);
    if (prn.fp != NULL) fclose(prn.fp);
    return 0;
}
