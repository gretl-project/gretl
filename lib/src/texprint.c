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
#include "johansen.h"
#include "texprint.h"

static char colspec[4][8];
static int use_custom;

int tex_using_custom_tabular (void)
{
    return use_custom;
}

const char *tex_column_format (int i)
{
    if (i >= 0 && i < 4) {
	return colspec[i];
    } else {
	return "";
    }
}

#define tex_screen_zero(x)  ((fabs(x) > 1.0e-17)? x : 0.0)

static void tex_modify_exponent (char *numstr)
{
    char *p = strchr(numstr, 'e');

    if (p != NULL) {
	int expon = atoi(p + 1);

	sprintf(p, "\\mbox{e%s%d}", (expon > 0)? "+" : "", expon);
    }
}

/* prints a floating point number as a TeX math string.  if tab != 0,
   print the sign in front, separated by a tab symbol (for
   equation-style regression printout).
*/

static void tex_print_float (double x, PRN *prn)
{
    char number[48];

    x = tex_screen_zero(x);

    sprintf(number, "%#.*g", GRETL_DIGITS, x);
    tex_modify_exponent(number);
    pprintf(prn, "%s", (x < 0.0)? (number + 1) : number);
}

static void tex_print_signed_float (double x, PRN *prn)
{
    char number[48];

    x = tex_screen_zero(x);

    sprintf(number, "%#.*g", GRETL_DIGITS, x);
    tex_modify_exponent(number);

    if (x < 0.0) {
	pprintf(prn, "$-$%s", number + 1);
    } else {
	pputs(prn, number);
    }
}

static void tex_print_math_float (double x, PRN *prn)
{
    char number[48];

    x = tex_screen_zero(x);

    sprintf(number, "%#.*g", GRETL_DIGITS, x);
    tex_modify_exponent(number);

    pprintf(prn, "$%s$", number);
}

/**
 * tex_float_string:
 * @x: value to be printed.
 * @targ: target string.
 *
 * Prints the value @x into @targ with @prec digits of precision.
 * The result is intended for use in a non-math environment: if
 * @x is negative, the minus sign is put into math mode.
 * 
 * Returns: @targ.
 */

char *tex_float_string (double x, int prec, char *targ)
{
    int offset = 0;

    *targ = '\0';

    x = tex_screen_zero(x);
    if (x < 0.0) {
	strcat(targ, "$-$");
	x = fabs(x);
	offset = 3;
    }

    sprintf(targ + offset, "%#.*g", prec, x);
    tex_modify_exponent(targ);

    return targ;
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

void tex_dcolumn_double (double x, char *numstr)
{
    double a;

    if (na(x)) {
	strcpy(numstr, "\\multicolumn{1}{c}{}");
	return;
    }

    x = screen_zero(x);
    a = fabs(x);

    sprintf(numstr, "%#.*g", GRETL_DIGITS, x);

    if (a != 0.0 && (a >= UPPER_F_LIMIT || a < LOWER_F_LIMIT)) {
	int expon;
	char *p, exponstr[8];

	p = strchr(numstr, 'e');
	expon = atoi(p + 2);
	strcpy(p, "\\mbox{e");
	sprintf(exponstr, "%s%02d}", (x > 10)? "+" : "-", expon);
	strcat(numstr, exponstr);
    } else {
	cut_extra_zero(numstr);
    }
}

static void tex_make_cname (char *cname, const char *src)
{
    char *p;
    unsigned char c;

    if (src == NULL || *src == '\0') return;

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
    } else if (!strcmp(src, "beta")) {
	strcpy(targ, "$\\beta$");
    } else if (!strcmp(src, "gamma")) {
	strcpy(targ, "$\\gamma$");
    } else if (!strcmp(src, "delta")) {
	strcpy(targ, "$\\delta$");
    } else if (!strcmp(src, "pi")) {
	strcpy(targ, "$\\pi$");
    } else if (!strcmp(src, "lambda")) {
	strcpy(targ, "$\\lambda$");
    }

    return (*targ != 0);
}

static void tex_arbond_coeff_name (char *targ, const char *src,
				   int inmath)
{
    char vname[VNAMELEN], vnesc[16];
    int lag;

    if (sscanf(src, "D%10[^(](%d)", vname, &lag) == 2) {
	tex_escape(vnesc, vname);
	if (!inmath) {
	    sprintf(targ, "$\\Delta \\mbox{\\rm %s}_{%d}$", vnesc, lag);
	} else {
	    sprintf(targ, "\\Delta \\mbox{\\rm %s}_{%d}", vnesc, lag);
	}
    } else {
	tex_escape(vnesc, src);
	if (inmath) {
	    sprintf(targ, "\\mbox{%s}", vnesc);
	} else {
	    strcpy(targ, vnesc);
	}
    }
}

static void tex_garch_coeff_name (char *targ, const char *src,
				  int inmath)
{
    char vname[VNAMELEN], vnesc[16];
    int lag;

    if (sscanf(src, "%15[^(](%d)", vname, &lag) == 2) {
	/* e.g. "alpha(0)" */
	if (!inmath) {
	    sprintf(targ, "$\\%s_{%d}$", vname, lag);
	} else {
	    sprintf(targ, "\\%s_{%d}", vname, lag);
	}
    } else {
	/* regular variable name */
	tex_escape(vnesc, src);
	if (inmath) {
	    sprintf(targ, "\\mbox{%s}", vnesc);
	} else {
	    strcpy(targ, vnesc);
	}
    }
}

static void tex_mp_coeff_name (char *targ, const char *src,
			       int inmath)
{
    char vname[VNAMELEN], vnesc[24];
    int power;

    tex_escape(vnesc, src);

    if (sscanf(vnesc, "%15[^^]^%d", vname, &power) == 2) {
	/* variable raised to some power */
	if (!inmath) {
	    sprintf(targ, "%s$^{%d}$", vname, power);
	} else {
	    sprintf(targ, "\\mbox{%s}^%d", vname, power);
	}
    } else {
	/* regular variable name */
	if (inmath) {
	    sprintf(targ, "\\mbox{%s}", vnesc);
	} else {
	    strcpy(targ, vnesc);
	}
    }
}

static void tex_arma_coeff_name (char *targ, const char *src,
				 int inmath)
{
    char vname[VNAMELEN], vnesc[32], texname[32];
    int i;

    if (sscanf(src, "phi_%d", &i)) {
	if (!inmath) {
	    sprintf(targ, "$\\phi_{%d}$", i);
	} else {
	    sprintf(targ, "\\phi_{%d}", i);
	}
    } else if (sscanf(src, "Phi_%d", &i)) {
	if (!inmath) {
	    sprintf(targ, "$\\Phi_{%d}$", i);
	} else {
	    sprintf(targ, "\\Phi_{%d}", i);
	}
    } else if (sscanf(src, "theta_%d", &i)) {
	if (!inmath) {
	    sprintf(targ, "$\\theta_{%d}$", i);
	} else {
	    sprintf(targ, "\\theta_{%d}", i);
	}
    } else if (sscanf(src, "Theta_%d", &i)) {
	if (!inmath) {
	    sprintf(targ, "$\\Theta_{%d}$", i);
	} else {
	    sprintf(targ, "\\Theta_{%d}", i);
	}
    } else if (sscanf(src, "%15[^(](-%d)", vname, &i) == 2) {
	if (!strcmp(vname, "y")) {
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
	    sprintf(targ, "%s$_{t-%d}$", texname, i);
	} else {
	    sprintf(targ, "%s_{t-%d}", texname, i);
	}
    } else {
	tex_escape(vnesc, src);
	strcpy(targ, vnesc);
    }
}

static void tex_lagname (char *s, const DATAINFO *pdinfo, int v)
{
    const char *lbl = VARLABEL(pdinfo, v);
    int gotit = 0;

    if (strlen(lbl) > 2) {
	char myvar[32], tmp[VNAMELEN];
	int lag;

	lbl += 2;
	if (sscanf(lbl, "%15[^(](t - %d)", tmp, &lag) == 2) {
	    tex_escape(myvar, tmp);
	    sprintf(s, "%s$_{t-%d}$", myvar, lag);
	    gotit = 1;
	}
    }
	
    if (!gotit) {
	tex_escape(s, pdinfo->varname[v]); 
    }
}

static void tex_vecm_varname (char *s, const DATAINFO *pdinfo, int v)
{
    const char *lbl = VARLABEL(pdinfo, v);
    int cvnum;
    int gotit = 0;

    if (sscanf(pdinfo->varname[v], "EC%d", &cvnum)) {
	sprintf(s, "EC%d$_{t-1}$", cvnum);
	gotit = 1;
    } else if (strlen(lbl) > 2) {
	char myvar[32], tmp[VNAMELEN];
	int lag;

	lbl += 2;
	if (sscanf(lbl, "d_%15[^(](t - %d)", tmp, &lag) == 2) {
	    tex_escape(myvar, tmp);
	    sprintf(s, "$\\Delta$%s$_{t-%d}$", myvar, lag);
	    gotit = 1;
	}
    }

    if (!gotit) {
	tex_escape(s, pdinfo->varname[v]); 
    }
}

static int tex_print_coeff_custom (const char *pname, const MODEL *pmod, 
				   int i, PRN *prn)
{
    double bi = pmod->coeff[i];
    double se = pmod->sderr[i];
    char fmt[12];
    double x;

    pprintf(prn, "%s & ", pname);

    if (colspec[0][0]) {
	/* coefficient */
	if (xna(bi)) {
	    pprintf(prn, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
	} else {
	    sprintf(fmt, "$%s$", colspec[0]);
	    pprintf(prn, fmt, pmod->coeff[i]);
	}
    }

    if (!colspec[1][0] && !colspec[2][0] && !colspec[3][0]) {
	pputs(prn, " \\\\\n");
	return 0;
    }

    if (colspec[1][0]) {
	if (colspec[0][0]) {
	    pputs(prn, " & ");
	}
	/* standard error */
	if (isnan(pmod->sderr[i]) || na(pmod->sderr[i])) {
	    pprintf(prn, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
	} else {
	    pprintf(prn, colspec[1], pmod->sderr[i]);
	}
    }

    if (!colspec[2][0] && !colspec[3][0]) {
	pputs(prn, " \\\\\n");
	return 0;
    }

    if (colspec[2][0]) {
	if (colspec[0][0] || colspec[1][0]) {
	    pputs(prn, " & ");
	}
	/* t-ratio */
	if (xna(bi) || xna(se)) {
	    pprintf(prn, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
	} else {
	    x = bi / se;
	    sprintf(fmt, "$%s$", colspec[2]);
	    pprintf(prn, fmt, x);
	}
    } 

    if (colspec[3][0]) {
	if (colspec[0][0] || colspec[1][0] || colspec[2][0]) {
	    pputs(prn, " & ");
	}
	/* p-value */
	if (xna(bi) || xna(se)) {
	    pprintf(prn, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
	} else {
	    x = coeff_pval(pmod, bi / se, pmod->dfd);
	    pprintf(prn, colspec[3], x);
	}
    }  

    pputs(prn, " \\\\\n");

    return 0;
}

int tex_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
		     int i, PRN *prn)
{
    char tmp[32], coeff[64], sderr[64], tratio[64], pval[64];
    double bi = pmod->coeff[i];
    double se = pmod->sderr[i];
    int j = i + 2;

    *tmp = 0;
    if (pmod->aux == AUX_ARCH) {
	tex_make_cname(tmp, pdinfo->varname[pmod->list[i+2]]);
    } else if (pmod->ci == NLS) {
	if (!tex_greek_param(tmp, pmod->params[i])) {
	    tex_escape(tmp, pmod->params[i]);
	}
    } else if (pmod->ci == ARMA) {
	tex_arma_coeff_name(tmp, pmod->params[i], 0);
    } else if (pmod->ci == GARCH) {
	tex_garch_coeff_name(tmp, pmod->params[i], 0);
    } else if (pmod->ci == VAR) {
	tex_lagname(tmp, pdinfo, pmod->list[j]);
    } else if (pmod->aux == AUX_VECM) {
	tex_vecm_varname(tmp, pdinfo, pmod->list[j]);
    } else if (pmod->ci == MPOLS && pmod->params != NULL) {
	tex_mp_coeff_name(tmp, pmod->params[i], 0);
    } else if ((pmod->ci == PROBIT || pmod->ci == LOGIT) &&
	       pmod->params != NULL) {
	tex_escape(tmp, pmod->params[i]);
    } else if (pmod->ci == PANEL) {
	tex_escape(tmp, pmod->params[i]);
    } else if (pmod->ci == ARBOND) {
	tex_arbond_coeff_name(tmp, pmod->params[i], 0);
    } else {
	tex_escape(tmp, pdinfo->varname[pmod->list[j]]);
    }

    if (use_custom) {
	return tex_print_coeff_custom(tmp, pmod, i, prn);
    }

    if (xna(bi)) {
	sprintf(coeff, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
    } else {
	tex_dcolumn_double(bi, coeff);
    }

    if (xna(se)) {
	sprintf(sderr, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
	sprintf(tratio, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
	sprintf(pval, "\\multicolumn{1}{c}{\\rm %s}", I_("undefined"));
    } else {
	tex_dcolumn_double(se, sderr);
	sprintf(tratio, "%.4f", bi / se);
	sprintf(pval, "%.4f", coeff_pval(pmod, bi / se, pmod->dfd));
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
    } else { 
	/* LOGIT, PROBIT */
	double *slopes = gretl_model_get_data(pmod, "slopes");
	char slope[32];

	if (pmod->list[j]) {
	    tex_dcolumn_double(slopes[i], slope);
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
		(pmod->list[j])? slope : "");
    }

    return 0;
}

void tex_custom_coeff_table_start (const char *col1, const char *col2,
				   PRN *prn)
{
    int i, ncols = 0;

    for (i=0; i<4; i++) {
	if (colspec[i][0]) ncols++;
    }

    pputs(prn, "\\vspace{1em}\n\n"
	  "\\begin{tabular}{l");

    for (i=0; i<ncols; i++) {
	pputs(prn, "r");
    }

    pputs(prn, "}\n");

    pprintf(prn, "\\multicolumn{1}{c}{%s} &\n", I_(col1));

    if (colspec[0][0]) {
	pprintf(prn, "\\multicolumn{1}{c}{%s}", I_(col2));
    }

    if (!colspec[1][0] && !colspec[2][0] && !colspec[3][0]) {
	pputs(prn, " \\\\\n");
	return;
    }

    if (colspec[1][0]) {
	if (colspec[0][0]) {
	    pputs(prn, " &\n");
	}
	pprintf(prn, "\\multicolumn{1}{c}{%s}", I_("Std.\\ Error"));
    }

    if (!colspec[2][0] && !colspec[3][0]) {
	pputs(prn, " \\\\\n");
	return;
    }

    if (colspec[2][0]) {
	if (colspec[0][0] || colspec[1][0]) {
	    pputs(prn, " &\n");
	}
	pprintf(prn, "\\multicolumn{1}{c}{%s}", I_("$t$-statistic"));
    }

    if (colspec[3][0]) {
	if (colspec[0][0] || colspec[1][0] || colspec[2][0]) {
	    pputs(prn, " &\n");
	}
	pprintf(prn, "\\multicolumn{1}{c}{%s}", I_("p-value"));
    }

    pputs(prn, " \\\\\n");
}

void tex_coeff_table_start (const char *col1, const char *col2,
			    int binary, PRN *prn)
{

    char pt;

    if (use_custom) {
	tex_custom_coeff_table_start(col1, col2, prn);
	return;
    }

    pt = get_local_decpoint();

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
	    (binary)? I_("Slope"): I_("p-value"),
	    (binary)? "$^*$" : "");
}

void tex_coeff_table_end (PRN *prn)
{
    if (use_custom) {
	pputs(prn, "\\end{tabular}\n\n");
    } else {
	pputs(prn, "\\end{tabular*}\n\n");
    }
}

void tex_print_VECM_omega (GRETL_VAR *vecm, const DATAINFO *pdinfo, PRN *prn)
{
    char vname[48];
    const int *list = vecm->jinfo->list;
    double x;
    int i, j;

    pprintf(prn, "%s\n\n", I_("Cross-equation covariance matrix"));
    pputs(prn, "\\vspace{1em}\n");

    pputs(prn, "\\begin{tabular}{");
    pputs(prn, "l");
    for (i=0; i<vecm->neqns; i++) {
	pputs(prn, "r");
    }
    pputs(prn, "}\n & ");

    for (i=0; i<vecm->neqns; i++) {
	tex_escape(vname, pdinfo->varname[list[i+1]]);
	pprintf(prn, "$\\Delta$%s ", vname);
	if (i == vecm->neqns - 1) {
	    pputs(prn, "\\\\\n");
	} else {
	    pputs(prn, "& ");
	}
    }
    pputc(prn, '\n');

    for (i=0; i<vecm->neqns; i++) {
	tex_escape(vname, pdinfo->varname[list[i+1]]);
	pprintf(prn, "$\\Delta$%s & ", vname);
	for (j=0; j<vecm->neqns; j++) {
	    x = gretl_matrix_get(vecm->S, i, j);
	    tex_print_math_float(x, prn);
	    if (j == vecm->neqns - 1) {
		pputs(prn, "\\\\\n");
	    } else {
		pputs(prn, " & ");
	    }
	}
    }

    pputs(prn, "\\end{tabular}\n\n");
    pputs(prn, "\\vspace{1em}\n");

    pputs(prn, "\\noindent\n");
    pprintf(prn, "%s = ", I_("determinant"));
    tex_print_math_float(exp(vecm->ldet), prn);
    pputs(prn, "\\\\\n");
}

void tex_print_VECM_coint_eqns (GRETL_VAR *vecm, const DATAINFO *pdinfo, PRN *prn)
{
    char s[32];
    JohansenInfo *jv = vecm->jinfo;
    int rows = gretl_matrix_rows(jv->Beta);
    int i, j;
    double x;

    pputs(prn, "\\noindent\n");
    pputs(prn, _("Cointegrating vectors"));
    if (jv->Bse != NULL) {
	pprintf(prn, " (%s)\n", _("standard errors in parentheses"));
    } 

    pputs(prn, "\n\\vspace{1em}\n");

    pputs(prn, "\\begin{tabular}{");
    pputs(prn, "l");
    for (i=0; i<jv->rank; i++) {
	pputs(prn, "r");
    }
    pputs(prn, "}\n");

    for (i=0; i<rows; i++) {
	if (i < jv->list[0]) {
	    tex_escape(s, pdinfo->varname[jv->list[i+1]]);
	    pprintf(prn, "%s$_{t-1}$ & ", s);
	} else if (jv->code == J_REST_CONST) {
	    pputs(prn, "const & ");
	} else if (jv->code == J_REST_TREND) {
	    pputs(prn, "trend & ");
	}

	/* coefficients */
	for (j=0; j<jv->rank; j++) {
	    x = gretl_matrix_get(jv->Beta, i, j);
	    if (jv->Bse == NULL) {
		x /= gretl_matrix_get(jv->Beta, j, j);
	    }
	    tex_print_signed_float(x, prn);
	    if (j == jv->rank - 1) {
		pputs(prn, "\\\\\n");
	    } else {
		pputs(prn, "& ");
	    }	    
	}

	if (jv->Bse != NULL) {
	    /* standard errors */
	    pputs(prn, " & ");
	    for (j=0; j<jv->rank; j++) {
		if (i < jv->rank) {
		    x = 0.0;
		} else {
		    x = gretl_matrix_get(jv->Bse, i - jv->rank, j);
		}
		pputc(prn, '(');
		tex_print_float(x, prn);
		pputc(prn, ')');
		if (j == jv->rank - 1) {
		    pputs(prn, "\\\\\n");
		} else {
		    pputs(prn, "& ");
		}
	    }
	}
    }

    pputs(prn, "\\end{tabular}\n\n\\vspace{1em}\n");
}

void tex_print_VAR_ll_stats (GRETL_VAR *var, PRN *prn)
{
    pprintf(prn, "\\noindent\n%s = ", I_("Log-likelihood"));
    tex_print_math_float(var->ll, prn);
    pputs(prn, "\\par\n");

    pprintf(prn, "\\noindent\n%s = ", I_("Determinant of covariance matrix"));
    tex_print_math_float(exp(var->ldet), prn);
    pputs(prn, "\\par\n");

    pprintf(prn, "\\noindent\n%s $= %.4f$ \\par\n", I_("AIC"), var->AIC);
    pprintf(prn, "\\noindent\n%s $= %.4f$ \\par\n", I_("BIC"), var->BIC);
    pprintf(prn, "\\noindent\n%s $= %.4f$ \\par\n", I_("HQC"), var->HQC);
}

static PRN *make_tex_prn (int ID, char *texfile,
			  int eqn, int doc)
{
    PrnFormat fmt = GRETL_FORMAT_TEX;
    PRN *prn;

    if (*texfile == '\0') {
	sprintf(texfile, "%s%s_%d.tex", gretl_user_dir(),
		(eqn)? "equation" : "model", ID);
    }

    prn = gretl_print_new_with_filename(texfile);
    if (prn != NULL) {
	if (eqn) {
	    fmt |= GRETL_FORMAT_EQN;
	}
	if (doc) {
	    fmt |= GRETL_FORMAT_DOC;
	}
	gretl_print_set_format(prn, fmt);
    }

    return prn;
}

/* mechanism for customizing gretl's tex preamble */

static char tex_preamble_file[MAXLEN];

#ifdef ENABLE_NLS
static const char *get_gretltex_local (void)
{
    static char localtex[32] = {0};
    char *lang = getenv("LANG");

    if (lang != NULL) {
	char lstr[3] = {0};

	strncat(lstr, lang, 2);
	sprintf(localtex, "gretlpre_%s.tex", lstr);
    }

    return localtex;
}
#endif

void set_gretl_tex_preamble (void)
{
    FILE *fp;
    const char *gretltex = "gretlpre.tex";
#ifdef ENABLE_NLS
    const char *localtex = get_gretltex_local();

    /* first choice: localized preamble file */
    sprintf(tex_preamble_file, "%s%s", gretl_user_dir(), localtex);
    fp = gretl_fopen(tex_preamble_file, "r");
    if (fp == NULL) {
	tex_preamble_file[0] = '\0';
    } else {
	fclose(fp);
	return;
    }    
#endif

    /* preamble file on disk */
    sprintf(tex_preamble_file, "%s%s", gretl_user_dir(), gretltex);
    fp = gretl_fopen(tex_preamble_file, "r");
    if (fp == NULL) {
	tex_preamble_file[0] = '\0';
    } else {
	fclose(fp);
    }
}

static int tex_use_utf;

void set_tex_use_utf (int s)
{
    tex_use_utf = s;
}

void gretl_tex_preamble (PRN *prn, int ams)
{
    FILE *fp = NULL;
    int userfile = 0;

    if (tex_preamble_file[0] != '\0') {
	fp = gretl_fopen(tex_preamble_file, "r");
	if (fp != NULL) {
	    char line[128];

	    while (fgets(line, sizeof line, fp)) {
		pputs(prn, line);
	    }
	    userfile = 1;
	    fclose(fp);
	}
    }

    if (!userfile) {
	pputs(prn, "\\documentclass[11pt]{article}\n");

#ifdef ENABLE_NLS
	if (tex_use_utf) {
	    pputs(prn, "\\usepackage{ucs}\n");
	    pputs(prn, "\\usepackage[utf8x]{inputenc}\n\n");
	} else {
	    pputs(prn, "\\usepackage[latin1]{inputenc}\n\n");
	}
#endif
	if (ams) {
	    pputs(prn, "\\usepackage{amsmath}\n\n");
	} else {
	    pputs(prn, "\\usepackage{dcolumn,longtable}\n\n");
	}

	pputs(prn, "\\begin{document}\n\n"
	      "\\thispagestyle{empty}\n\n");
    }
}

#define MAXCOEFF 4

/**
 * tex_print_equation:
 * @pmod:  pointer to gretl MODEL struct.
 * @pdinfo:  information regarding the data set.
 * @opt: can include %OPT_S for a standalone document, and
 * %OPT_T to print t-ratios rather than standard errors.
 * @prn: gretl printing struct.
 *
 * Prints to @prn a gretl model in the form of a LaTeX equation, either as
 * a stand-alone document or as a fragment of LaTeX source for
 * insertion into a document.
 * 
 * Returns: 0 on successful completion.
 */

int tex_print_equation (const MODEL *pmod, const DATAINFO *pdinfo, 
			gretlopt opt, PRN *prn)
{
    double x;
    char tmp[48];
    int i, nc = pmod->ncoeff;
    int split = 0, offvar = 0;
    int cchars = 0, ccount = 0;

    if (pmod->ci == POISSON) {
	offvar = gretl_model_get_int(pmod, "offset_var");
	if (offvar > 0) {
	    nc++;
	}
    }

    split = (nc > MAXCOEFF);

    if (opt & OPT_S) {
	gretl_tex_preamble(prn, 1);
    } else{
	pputs(prn, "%%% the following needs the amsmath LaTeX package\n\n");
    }

    /* initial setup */
    pputs(prn, "\\begin{gather}\n");
    if (split) {
	pputs(prn, "\\begin{split}\n");
    }

    /* dependent variable */
    *tmp = '\0';
    if (pmod->depvar != NULL) {
	tex_escape(tmp, pmod->depvar);
    } else if (pmod->ci == POISSON) {
	char vname[32];

	tex_escape(vname, pdinfo->varname[pmod->list[1]]);
	sprintf(tmp, "log(%s)", vname);
    } else {
	i = gretl_model_get_depvar(pmod);
	tex_escape(tmp, pdinfo->varname[i]);
    }

    if (pmod->ci == ARBOND) {
	pprintf(prn, "\\widehat{\\Delta \\rm %s} %s= \n", tmp, (split? "&" : ""));
    } else {
	pprintf(prn, "\\widehat{\\rm %s} %s= \n", tmp, (split? "&" : ""));
    }

    if (pmod->ci == GARCH) {
	nc -= (1 + pmod->list[1] + pmod->list[2]);
    } else if (pmod->ci == PANEL) {
	nc = pmod->list[0] - 1;
    }

    /* coefficients times indep vars */
    for (i=0; i<nc; i++) {
	if (offvar > 0 && i == nc - 1) {
	    pputc(prn, '+');
	    tex_print_float(1.0, prn);
	} else {
	    if (opt & OPT_T) {
		x = pmod->coeff[i] / pmod->sderr[i];
		pprintf(prn, "%s\\underset{(%.3f)}{", 
			(pmod->coeff[i] < 0.0)? "-" :
			(i > 0)? "+" : "", x);
	    } else {
		pprintf(prn, "%s\\underset{(%.5g)}{", 
			(pmod->coeff[i] < 0.0)? "-" :
			(i > 0)? "+" : "", pmod->sderr[i]);
	    }
	    tex_print_float(pmod->coeff[i], prn);
	    pputc(prn, '}');
	}
	if (i > 0 || pmod->ifc == 0) {
	    pputs(prn, "\\,");
	    if (pmod->ci == ARMA) {
		cchars += strlen(pmod->params[i]);
		tex_arma_coeff_name(tmp, pmod->params[i], 1);
		pputs(prn, tmp);
	    } else if (pmod->ci == GARCH) {
		cchars += strlen(pmod->params[i]);
		tex_garch_coeff_name(tmp, pmod->params[i], 1);
		pputs(prn, tmp);
	    } else if (pmod->ci == PANEL) {
		cchars += strlen(pmod->params[i]);
		tex_escape(tmp, pmod->params[i]);
		pprintf(prn, "\\mbox{%s}", tmp);
	    } else if (pmod->ci == ARBOND) {
		if (strcmp(pmod->params[i], "const")) {
		    cchars += strlen(pmod->params[i]);
		    tex_arbond_coeff_name(tmp, pmod->params[i], 1);
		    pputs(prn, tmp);
		}
	    } else if (pmod->ci == MPOLS && pmod->params != NULL) {
		cchars += strlen(pmod->params[i]);
		tex_mp_coeff_name(tmp, pmod->params[i], 1);
		pputs(prn, tmp);
	    } else if ((pmod->ci == PROBIT || pmod->ci == LOGIT) &&
		       pmod->params != NULL) {
		cchars += strlen(pmod->params[i]);
		tex_escape(tmp, pmod->params[i]);
		pprintf(prn, "\\mbox{%s}", tmp);
	    } else if (offvar > 0 && i == nc - 1) {
		cchars += strlen(pdinfo->varname[offvar]);
		tex_escape(tmp, pdinfo->varname[offvar]);
		pprintf(prn, "\\mbox{%s}", tmp);
	    } else {
		cchars += strlen(pdinfo->varname[pmod->list[i+2]]);
		tex_escape(tmp, pdinfo->varname[pmod->list[i+2]]);
		pprintf(prn, "\\mbox{%s}", tmp);
	    }
	}
	ccount++;
	if (split && (cchars > 30 || ccount > 3)) {
	    pputs(prn, "\\\\\n& ");
	    cchars = ccount = 0;
	} else {
	    pputc(prn, '\n');
	}
    }

    if (split) {
	pputs(prn, "\\end{split}\n");
    }

    pputs(prn, " \\notag \\\\\n");

    if (pmod->ci == GARCH) {
	int q = pmod->list[1];
	int p = pmod->list[2];
	int r = pmod->list[0] - 4;

	if (opt & OPT_T) {
	    x = pmod->coeff[r] / pmod->sderr[r];
	    pprintf(prn, "\\hat{\\sigma}^2_t = \\underset{(%.3f)}{%g} ", 
		    x, pmod->coeff[r]);
	} else {
	    pprintf(prn, "\\hat{\\sigma}^2_t = \\underset{(%.5g)}{%g} ", 
		    pmod->sderr[r], pmod->coeff[r]);
	}	    

	for (i=1; i<=q; i++) {
	    if (opt & OPT_T) {
		x = pmod->coeff[r+i] / pmod->sderr[r+i];
		pprintf(prn, "%s\\underset{(%.3f)}{", 
			(pmod->coeff[r+i] < 0.0)? "-" : "+", x);
	    } else {
		pprintf(prn, "%s\\underset{(%.5g)}{", 
			(pmod->coeff[r+i] < 0.0)? "-" : "+", 
			pmod->sderr[r+i]);
	    }		
	    tex_print_float(pmod->coeff[r+i], prn);
	    pputs(prn, "}\\,");
	    pprintf(prn, "\\varepsilon^2_{t-%d}", i);
	}

	for (i=1; i<=p; i++) {
	    if (opt & OPT_T) {
		x = pmod->coeff[q+r+i] / pmod->sderr[q+r+i];
		pprintf(prn, "%s\\underset{(%.3f)}{", 
			(pmod->coeff[q+r+i] < 0.0)? "-" : "+", x);
	    } else {
		pprintf(prn, "%s\\underset{(%.5g)}{", 
			(pmod->coeff[q+r+i] < 0.0)? "-" : "+", 
			pmod->sderr[q+r+i]);
	    }		
	    tex_print_float(pmod->coeff[q+r+i], prn);
	    pputs(prn, "}\\,");
	    pprintf(prn, "\\sigma^2_{t-%d}", i);
	}

	pputs(prn, "\\notag \\\\\n");
    }

    pprintf(prn, "T = %d ", pmod->nobs);

    /* additional info (R^2 etc) */
    if (pmod->ci == LAD) { 
	sprintf(tmp, "%g", pmod->rho);
	tex_modify_exponent(tmp);
	pprintf(prn, "\\quad \\sum |\\hat{u}_t| = %s", tmp);
    } else {
	if (!na(pmod->adjrsq)) {
	    pprintf(prn, "\\quad \\bar{R}^2 = %.4f ", pmod->adjrsq);
	} else if (!na(pmod->lnL)) {
	    pprintf(prn, "\\quad \\mbox{ln}L = %.4f ", pmod->lnL);
	}

	if (pmod->ci != LOGIT && pmod->ci != PROBIT && !na(pmod->fstt)) {
	    sprintf(tmp, "%.5g", pmod->fstt);
	    tex_modify_exponent(tmp);
	    pprintf(prn, "\\quad F(%d,%d) = %s ", 
		    pmod->dfn, pmod->dfd, tmp);
	}

	if (!na(pmod->sigma)) {
	    sprintf(tmp, "%.5g", pmod->sigma);
	    tex_modify_exponent(tmp);
	    pprintf(prn, "\\quad \\hat{\\sigma} = %s", tmp);
	}

	if (!na(gretl_model_get_double(pmod, "rho_in"))) {
	    double r = gretl_model_get_double(pmod, "rho_in");

	    sprintf(tmp, "%.5g", r);
	    tex_modify_exponent(tmp);
	    pprintf(prn, " \\quad \\rho = %s", tmp);
	}
    }

    pputs(prn, "\\notag \\\\\n");
    pprintf(prn, "\\centerline{(%s)} \\notag\n",
	    (opt & OPT_T)? I_("$t$-statistics in parentheses") :
	    I_("standard errors in parentheses"));
    pputs(prn, "\\end{gather}\n");

    if (opt & OPT_S) {
	pputs(prn, "\n\\end{document}\n");
    }

    return 0;
}

/**
 * tex_print_model:
 * @pmod:  pointer to gretl MODEL struct.
 * @pdinfo: information regarding the data set.
 * @opt: may include %OPT_T in case of printing a model
 * in equation format, to use t-ratios instead of standard
 * errors.
 * @prn: gretl printing struct.
 *
 * Prints to @prn a gretl model in the form of either a LaTeX 
 * table or an equation, and either as a stand-alone document or 
 * as a fragment of LaTeX source for insertion into a document.
 *
 * The options are read from the format field of @prn --
 * gretl_print_set_format().
 * 
 * Returns: 0 on successful completion.
 */

int tex_print_model (MODEL *pmod, const DATAINFO *pdinfo, 
		     gretlopt opt, PRN *prn)
{
    int ret;

    if (tex_doc_format(prn)) {
	opt |= OPT_S;
    }

    if (tex_eqn_format(prn)) { 
	ret = tex_print_equation(pmod, pdinfo, opt, prn);
    } else {
	ret = printmodel(pmod, pdinfo, OPT_NONE, prn);
    }

    return ret;
}

/**
 * texprint:
 * @pmod: pointer to model.
 * @pdinfo: information regarding the data set.
 * @texfile: name of file to save.
 * @opt: if opt & %OPT_O, complete doc, else fragment;
 * if opt & %OPT_E print as equation, otherwise use tabular
 * format; if opt & %OPT_T show t-ratios rather than standard
 * errors when printing in equation format.
 *
 * Prints to file a gretl model in the form of a LaTeX table or
 * equation, either as a stand-alone document or as a fragment 
 * of LaTeX source for insertion into a document.
 * 
 * Returns: 0 on successful completion, 1 on error.
 */

int texprint (MODEL *pmod, const DATAINFO *pdinfo, char *texfile, 
	      gretlopt opt)
{
    PRN *prn;
    int eqn = (opt & OPT_E);
    int doc = (opt & OPT_O);
    int err = 0;

    prn = make_tex_prn(pmod->ID, texfile, eqn, doc);
    if (prn == NULL) {
	err = 1;
    } else {
	err = tex_print_model(pmod, pdinfo, opt, prn);
	gretl_print_destroy(prn);
    }

    return err;
}

/**
 * tex_print_obs_marker:
 * @t: observation number.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Print a string (label, date or obs number) representing the given @t.
 */

void tex_print_obs_marker (int t, const DATAINFO *pdinfo, PRN *prn)
{
    if (pdinfo->markers) { 
	pprintf(prn, "\\texttt{%s} ", pdinfo->S[t]); 
    } else {
	char tmp[OBSLEN]; 

	ntodate(tmp, t, pdinfo);
	pprintf(prn, "%8s ", tmp);
    }
}

static int check_colspec (const char *s)
{
    const char *ok = "eEfgG";
    int w = 0, p = 0;
    char c = 0;
    int err = 1;

    /* blank is OK */
    if (*s == '\0') {
	return 0;
    }

    if (*s != '%') {
	return 1;
    }

    s++;

    if (*s == '#') { /* OK */
	s++;
    }

    if (sscanf(s, "%d.%d%c", &w, &p, &c) == 3) {
	if (w != 0 && p > 0 && strchr(ok, c)) {
	    err = 0;
	}
    } else if (sscanf(s, "%d%c", &w, &c) == 2) {
	if (w != 0 && strchr(ok, c)) {
	    err = 0;
	}
    } else if (sscanf(s, ".%d%c", &p, &c) == 2) {
	if (p > 0 && strchr(ok, c)) {
	    err = 0;
	}
    } else if (sscanf(s, "%c", &c) == 1) {
	if (strchr(ok, c)) {
	    err = 0;
	}
    } 

    return err;
}

/**
 * set_tex_param_format:
 * @s: stylized format string.
 *
 * Sets the format with which parameters will be printed, when
 * producing TeX tabular output.
 */

void set_tex_param_format (const char *s)
{
    const char *p = s;
    int i, n = 0;
    int err = 0;

    if (s == NULL) {
	use_custom = 0;
	return;
    }

    for (i=0; i<4; i++) {
	colspec[i][0] = '\0';
    }

    i = 0;

    while (i < 4) {
	if (*s == '|' || *s == '\0') {
	    if (n > 7) {
		n = 7;
	    }
	    strncat(colspec[i], p, n);
	    fprintf(stderr, "spec %d = '%s'\n", i, colspec[i]);
	    err = check_colspec(colspec[i]);
	    if (err || *s == '\0') {
		break;
	    }
	    p = s + 1;
	    i++;
	    n = 0;
	} else {
	    n++;
	}
	s++;
    }

    if (!err) {
	/* all columns can't be blank */
	n = 0;
	for (i=0; i<4; i++) {
	    if (colspec[i][0] != '\0') n++;
	}
	if (n == 0) {
	    err = 1;
	}
    }

    if (err) {
	for (i=0; i<4; i++) {
	    colspec[i][0] = '\0';
	}
	use_custom = 0;
    } else {
	use_custom = 1;
    }
}
