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
#include <stdarg.h>

#define PRECISION 6
/* #define TEX_ASCII 1 */

#ifdef TEX_ASCII

struct translation {
    unsigned code;		/* code being translated */
    const char *string;		/* translation string */
};

static struct translation const tex_translations [] = {
    {191, "?`"},		/* inverted question mark */
    {192, "\\`A"},		/* capital A with grave accent */
    {193, "\\'A"},		/* capital A with acute accent */
    {194, "\\^A"},		/* capital A with circumflex accent */
    {195, "\\~A"},		/* capital A with tilde */
    {196, "\\\"A"},		/* capital A diaeresis */
    {197, "\\AA{}"},		/* capital A with ring above */
    {198, "\\AE{}"},		/* capital diphthong A with E */
    {199, "\\c{C}"},		/* capital C with cedilla */
    {200, "\\`E"},		/* capital E with grave accent */
    {201, "\\'E"},		/* capital E with acute accent */
    {202, "\\^E"},		/* capital E with circumflex accent */
    {203, "\\\"E"},		/* capital E with diaeresis */
    {204, "\\`I"},		/* capital I with grave accent */
    {205, "\\'I"},		/* capital I with acute accent */
    {206, "\\^I"},		/* capital I with circumflex accent */
    {207, "\\\"I"},		/* capital I with diaeresis */
    {208, ""},    		/* capital Eth, Icelandic */
    {209, "\\~N"},		/* capital N with tilde */
    {210, "\\`O"},		/* capital O with grave accent */
    {211, "\\'O"},		/* capital O with acute accent */
    {212, "\\^O"},		/* capital O with circumflex accent */
    {213, "\\~O"},		/* capital O with tilde */
    {214, "\\\"O"},		/* capital O with diaeresis */
    {215, "\\times "},    	/* multiply sign */
    {216, "\\O{}"},		/* capital O with oblique stroke */
    {217, "\\`U"},		/* capital U with grave accent */
    {218, "\\'U"},		/* capital U with acute accent */
    {219, "\\^U"},		/* capital U with circumflex accent */
    {220, "\\\"U"},		/* capital U with diaeresis */
    {221, "\\'Y"},		/* capital Y with acute accent */
    {222, ""},    		/* capital THORN, Icelandic */
    {223, "\\ss{}"},		/* small german sharp s */
    {224, "\\`a"},		/* small a with grave accent */
    {225, "\\'a"},		/* small a with acute accent */
    {226, "\\^a"},		/* small a with circumflex accent */
    {227, "\\~a"},		/* small a with tilde */
    {228, "\\\"a"},		/* small a with diaeresis */
    {229, "\\aa{}"},		/* small a with ring above */
    {230, "\\ae{}"},		/* small diphthong a with e */
    {231, "\\c{c}"},		/* small c with cedilla */
    {232, "\\`e"},		/* small e with grave accent */
    {233, "\\'e"},		/* small e with acute accent */
    {234, "\\^e"},		/* small e with circumflex accent */
    {235, "\\\"e"},		/* small e with diaeresis */
    {236, "\\`{\\i}"},		/* small i with grave accent */
    {237, "\\'{\\i}"},		/* small i with acute accent */
    {238, "\\^{\\i}"},		/* small i with circumflex accent */
    {239, "\\\"{\\i}"},		/* small i with diaeresis */
    {240, ""},    		/* small eth, Icelandic */
    {241, "\\~n"},		/* small n with tilde */
    {242, "\\`o"},		/* small o with grave accent */
    {243, "\\'o"},		/* small o with acute accent */
    {244, "\\^o"},		/* small o with circumflex accent */
    {245, "\\~o"},		/* small o with tilde */
    {246, "\\\"o"},		/* small o with diaeresis */
    {247, ""},    		/* division sign */
    {248, "\\o{}"},		/* small o with oblique stroke */
    {249, "\\`u"},		/* small u with grave accent */
    {250, "\\'u"},		/* small u with acute accent */
    {251, "\\^u"},		/* small u with circumflex accent */
    {252, "\\\"u"},		/* small u with diaeresis */
    {253, "\\'y"},		/* small y with acute accent */
    {254, ""},    		/* small thorn, Icelandic */
    {255, "\\\"y"},		/* small y with diaeresis */
    {0, NULL}
};

static char *latin1_to_tex (const char *src)
{
    static char *dest;
    char add[2];
    size_t buflen, len = strlen(src);
    unsigned char c;

    buflen = 2 * len;
    if (buflen < 16) buflen = 16;
    dest = malloc(buflen);
    if (dest == NULL) return NULL;

    *dest = 0;
    add[1] = 0;

    while ((c = *src)) {
	if (c > 255) continue;
	if (buflen - strlen(dest) <= 8) {
	    buflen *= 2;
	    dest = realloc(dest, buflen);
	    if (dest == NULL) return NULL;
	}
	if (c >= 191) {
	    strcat(dest, tex_translations[c - 191].string);
	} else {
	    *add = c;
	    strcat(dest, add);
	}
	src++;
    }
    return dest;
}

#endif /* TEX_ASCII */

static int recode_pprintf (PRN *prn, const char *template, ...)
{
    va_list args;
    char tmp[1024]; /* watch out, can't use for big chunks */
    int err;
#ifdef TEX_ASCII
    char *trans;
#endif

    va_start(args, template);
    vsprintf(tmp, template, args);
    va_end(args);

#ifdef TEX_ASCII
    trans = latin1_to_tex(tmp);
    if (trans == NULL) return 1;
    err = pprintf(prn, "%s", trans);
    free(trans);    
#else
    err = pprintf(prn, "%s", tmp);
#endif

    return err;
}

/* ......................................................... */

static void tex_print_float (const double x, const int tab, PRN *prn)
     /* prints a floating point number as a TeX math string.
	if tab != 0, print the sign in front, separated by a
	tab symbol (for equation-style regression printout).
     */
{
    char number[16];

    sprintf(number, "%#.*g", PRECISION, x);

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

#define UPPER_F_LIMIT (pow(10, PRECISION))
#define LOWER_F_LIMIT (pow(10, -4))

static void dcol_float (double xx, char *numstr)
{
    double a = fabs(xx);

    sprintf(numstr, "%#.*g", PRECISION, xx);

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

static void tex_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			     const int c, PRN *prn)
{
    char tmp[16], coeff[32], stderr[32];
    double t_ratio = pmod->coeff[c-1] / pmod->sderr[c-1];
    
    tmp[0] = '\0';
    tex_escape(tmp, pdinfo->varname[pmod->list[c]]);

    dcol_float(pmod->coeff[c-1], coeff);
    dcol_float(pmod->sderr[c-1], stderr);

    pprintf(prn, "%s &\n"
	    "  %s &\n"
	    "    %s &\n"
	    "      %.4f &\n"
	    "        %.4f \\\\\n",  
	    tmp,
	    coeff,
	    stderr,
	    t_ratio,
	    tprob(t_ratio, pmod->dfd));	
}

/* ......................................................... */

static int make_texfile (const PATHS *ppaths, const int model_count,
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
			const int standalone, PRN *prn)
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
    pprintf(prn, "$T = %d,\\, \\bar{R}^2 = %.3f,\\, F(%d,%d) = %#g,\\, "
	    "\\hat{\\sigma} = %#g$\n",
	    pmod->nobs, pmod->adjrsq, pmod->dfn, 
	    pmod->dfd, pmod->fstt, pmod->sigma);

    recode_pprintf(prn, "\n(%s)\n\\end{center}\n", 
		   I_("$t$-statistics in parentheses"));

    if (standalone) 
	pprintf(prn, "\n\\end{document}\n");

    return 0;
}

/* stats printing functions */

static void tex_depvarstats (const MODEL *pmod, PRN *prn,
			     char *x1str, char *x2str)
{
    dcol_float(pmod->ybar, x1str);
    dcol_float(pmod->sdy, x2str);
    recode_pprintf(prn, "%s & %s \\\\\n %s & %s \\\\\n",
		   I_("Mean of dependent variable"), x1str,
		   I_("S.D. of dependent variable"), x2str);
} 

static void tex_essline (const MODEL *pmod, PRN *prn,
			 char *x1str, char *x2str)
{
    dcol_float(pmod->ess, x1str);
    dcol_float(pmod->sigma, x2str);
    recode_pprintf(prn, "%s & %s \\\\\n %s ($\\hat{\\sigma}$) & %s \\\\\n",
		   I_("Sum of squared residuals"), x1str,
		   I_("Standard error of residuals"), x2str);
} 

static void tex_rsqline (const MODEL *pmod, PRN *prn,
			 char *x1str, char *x2str)
{
    dcol_float(pmod->rsq, x1str);
    dcol_float(pmod->adjrsq, x2str);
    recode_pprintf(prn, "%s & %s \\\\\n %s & %s \\\\\n",
		   I_("Unadjusted $R^2$"), x1str, 
		   I_("Adjusted $\\bar{R}^2$"), x2str);
} 

static void tex_fline (const MODEL *pmod, PRN *prn,
		       char *x1str, char *x2str)
{
    dcol_float(pmod->fstt, x1str);
    dcol_float(fdist(pmod->fstt, pmod->dfn, pmod->dfd), x2str);
    recode_pprintf(prn, "%s (%d, %d) & %s \\\\\n %s & %s \\\\\n",
		   I_("F-statistic"), pmod->dfn, pmod->dfd, x1str,
		   I_("p-value for F()"), x2str);
} 

static void tex_dwline (const MODEL *pmod, PRN *prn,
			char *x1str, char *x2str)
{
    dcol_float(pmod->dw, x1str);
    dcol_float(pmod->rho, x2str);
    recode_pprintf(prn, "%s & %s \\\\\n %s ($\\hat{\\rho}$) & %s \n",
		   I_("Durbin--Watson statistic"), x1str, 
		   I_("First-order autocorrelation coeff."), x2str);
} 

static void tex_print_aicetc (const MODEL *pmod, PRN *prn)
{
    recode_pprintf(prn, 
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

/* ......................................................... */

static const char *tex_estimator_string (int ci)
{
    if (ci == OLS || ci == VAR) return I_("OLS");
    else if (ci == WLS) return I_("WLS"); 
    else if (ci == ARCH) return I_("WLS (ARCH)");
    else if (ci == CORC) return I_("Cochrane--Orcutt");
    else if (ci == HILU) return I_("Hildreth--Lu");
    else if (ci == TSLS) return I_("TSLS");
    else if (ci == HSK) return I_("Heteroskedasticity");
    else if (ci == AR) return I_("AR");
    else if (ci == HCCM) return I_("HCCM");
    else if (ci == PROBIT) return I_("Probit");
    else if (ci == LOGIT) return I_("Logit");
    else if (ci == POOLED) return I_("Pooled OLS");
    else return "";
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
		     const int standalone, PRN *prn)
{
    int i, ncoeff = pmod->list[0];
    int t1 = pmod->t1, t2 = pmod->t2;
    char tmp[16], x1str[32], x2str[32];
    char startdate[9], enddate[9], topstr[128];
    char pt = get_local_decpoint();

    modelprint_setup_obs(pmod, &t1, &t2);

    ncoeff = pmod->list[0];
    ntodate(startdate, t1, pdinfo);
    ntodate(enddate, t2, pdinfo);

    if (standalone) {
	pprintf(prn, "\\documentclass{article}\n"
		"\\usepackage{dcolumn}\n");
#ifndef TEX_ASCII
	pprintf(prn, "\\usepackage[T1]{fontenc}\n");
#endif
	pprintf(prn, "\\begin{document}\n\n"
		"\\thispagestyle{empty}\n");
    }

    pprintf(prn, "\\begin{center}\n");

    tex_escape(tmp, pdinfo->varname[pmod->list[1]]);

    sprintf(topstr, I_("Model %d: %s estimates using the %d observations %s--%s"),
	    pmod->ID, tex_estimator_string(pmod->ci), pmod->nobs,
	    startdate, enddate);
    
    recode_pprintf(prn, "\\textbf{%s}\\\\\n%s: %s", 
		   topstr, I_("Dependent variable"), tmp);

    if (pmod->ci == WLS || pmod->ci == ARCH) { 
	tex_escape(tmp, pdinfo->varname[pmod->nwt]);
	recode_pprintf(prn, "\\\\\n%s: %s\n\n", I_("Variable used as weight"), tmp);
    } else
	pprintf(prn, "\n\n");

    /* start table of coefficients */
    recode_pprintf(prn, "\\vspace{1em}\n\n"
		   "\\begin{tabular*}{\\textwidth}"
		   "{@{\\extracolsep{\\fill}}\n"
		   "l%% col 1: varname\n"
		   "  D{%c}{%c}{-1}%% col 2: coeff\n"
		   "    D{%c}{%c}{-1}%% col 3: stderr\n"
		   "      D{%c}{%c}{-1}%% col 4: t-stat\n"
		   "        D{%c}{%c}{4}}%% col 5: p-value\n"
		   "Variable &\n"
		   "  \\multicolumn{1}{c}{%s} &\n"
		   "    \\multicolumn{1}{c}{%s} &\n"
		   "      \\multicolumn{1}{c}{%s} &\n"
		   "        \\multicolumn{1}{c}{%s} \\\\[1ex]\n",
		   pt, pt, pt, pt, pt, pt, pt, pt,
		   I_("Coefficient"), I_("Std.\\ Error"), 
		   I_("$t$-statistic"), I_("p-value"));

    if (pmod->ifc) {
	tex_print_coeff(pdinfo, pmod, ncoeff, prn);
	ncoeff--;
    }
    for (i=2; i<=ncoeff; i++) tex_print_coeff(pdinfo, pmod, i, prn);

    pprintf(prn, "\\end{tabular*}\n\n");

    if (pmod->aux == AUX_ARCH || pmod->aux == AUX_ADF)
	goto stop_tex;

    if (pmod->ci == CORC || pmod->ci == HILU)
	recode_pprintf(prn, "\\vspace{1em}\n%s ($\\rho=%g$)\n\n", 
		       I_("Statistics based on quasi-differenced data"), 
		       pmod->rhot[1]);

    pprintf(prn, 
	    "\\vspace{1em}\n\n"
	    "\\begin{tabular}{lD{%c}{%c}{-1}}\n",
	    pt, pt, pt, pt);

    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG ||
	pmod->aux == AUX_WHITE || pmod->aux == AUX_AR) {
	tex_rsqline(pmod, prn, x1str, x2str);
    } else {
	if (pmod->ci != CORC && pmod->ci != HILU) {
	    tex_depvarstats(pmod, prn, x1str, x2str);
	}
	tex_essline(pmod, prn, x1str, x2str);
	tex_rsqline(pmod, prn, x1str, x2str);
	tex_fline(pmod, prn, x1str, x2str);
	tex_dwline(pmod, prn, x1str, x2str);
    }

    pprintf(prn, "\\end{tabular}\n\n");
    
    tex_print_aicetc(pmod, prn);

 stop_tex:
    pprintf(prn, "\n\\end{center}\n");

    if (standalone) 
	pprintf(prn, "\n\\end{document}\n");

    return 0;
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
	      const int model_count, int oflag)
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
	      const int model_count, int oflag)
{
    PRN prn;

    if (make_texfile(ppaths, model_count, 1, texfile, &prn))
	return 1;

    tex_print_equation(pmod, pdinfo, oflag, &prn);
    if (prn.fp != NULL) fclose(prn.fp);
    return 0;
}
