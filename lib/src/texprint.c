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

/* texprint.c for gretl - LaTeX output from modeling */

#include "libgretl.h"
#include "var.h"
#include "johansen.h"
#include "texprint.h"

static char colspec[4][8];
static int use_custom;
static int use_pdf;

#define tex_screen_zero(x)  ((fabs(x) > 1.0e-17)? x : 0.0)

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

void set_tex_use_pdf (const char *prog)
{
    const char *p = strrslash(prog);
    char test[4];

    /* looking for "pdflatex", possibly preceded by
       an absolute path */

    p = (p == NULL)? prog : p + 1;

    *test = '\0';
    strncat(test, p, 3);
    gretl_lower(test);
    use_pdf = !strcmp(test, "pdf");
}

int get_tex_use_pdf (void)
{
    return use_pdf;
}

static const char *tex_greek_var (const char *s)
{
    if (!strcmp(s, "alpha")) {
	return "\\alpha";
    } else if (!strcmp(s, "beta")) {
	return "\\beta";
    } else if (!strcmp(s, "gamma")) {
	return "\\gamma";
    } else if (!strcmp(s, "delta")) {
	return "\\delta";
    } else if (!strcmp(s, "epsilon")) {
	return "\\epsilon";
    } else if (!strcmp(s, "chi")) {
	return "\\chi";
    } else if (!strcmp(s, "pi")) {
	return "\\pi";
    } else if (!strcmp(s, "phi")) {
	return "\\phi";
    } else if (!strcmp(s, "psi")) {
	return "\\psi";
    } else if (!strcmp(s, "lambda")) {
	return "\\lambda";
    }

    return NULL;
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

    if (src == NULL) {
	fprintf(stderr, "tex_escape: src is NULL\n");
	*p = '\0';
	return p;
    }

    while (*src) {
	if (*src == '$' || *src == '&' || *src == '_' ||
	    *src == '%' || *src == '#')
	    *targ++ = '\\';
	*targ++ = *src++;
    }

    *targ = '\0';

    return p;
}

/**
 * tex_escape_new:
 * @src: source string.
 *
 * Returns: a copy of @src in which any characters that require
 * escaping are preceded by a backslash.
 */

char *tex_escape_new (const char *src)
{
    const char *s = src;
    char *ret = NULL;
    int i, len = 0;

    if (src == NULL) {
	fprintf(stderr, "tex_escape: src is NULL\n");
	return gretl_strdup("");
    }

    while (*s) {
	if (*s == '$' || *s == '&' || *s == '_' ||
	    *s == '%' || *s == '#') {
	    len++;
	}
	len++;
    }

    ret = calloc(1, len + 1);
    s = src;
    i = 0;

    while (*s) {
	if (*s == '$' || *s == '&' || *s == '_' ||
	    *s == '%' || *s == '#') {
	    ret[i++] = '\\';
	}
	ret[i++] = *s;
    }

    return ret;
}

static int tex_math_pname (char *targ, const char *s)
{
    char base[16], op[2], mod[8];
    int n;

    n = sscanf(s, "%15[^_^]%1[_^]%7s", base, op, mod);

    if (n == 3 && (*mod == '{' || isdigit(*mod))) {
	const char *tgreek = tex_greek_var(base);
	const char *tbase = (tgreek != NULL)? tgreek : base;

	if (*mod == '{') {
	    sprintf(targ, "$%s%s%s$", tbase, op, mod);
	} else {
	    sprintf(targ, "$%s%s{%s}$", tbase, op, mod);
	}
	return 1;
    }

    return 0;
}

/**
 * tex_escape_special:
 * @targ: target string (must be pre-allocated)
 * @src: source string.
 *
 * Copies from @src to @targ, escaping characters in @src that are
 * special to TeX (by inserting a leading backslash).  Unlike
 * tex_escape(), this function does not mess with '$' in the
 * source string, and it attempts to handle greek letters
 * correctly.
 *
 * Returns: the transformed copy of the string.
 */

char *tex_escape_special (char *targ, const char *src)
{
    const char *tgreek;
    char *p = targ;

    if (strchr(src, '$')) {
	/* don't mess with it */
	strcpy(targ, src);
	return targ;
    }

    tgreek = tex_greek_var(src);

    if (tgreek != NULL) {
	sprintf(targ, "$%s$", tgreek);
    } else if (tex_math_pname(targ, src)) {
	; /* handled */
    } else {
	/* regular escape routine */
	while (*src) {
	    if (*src == '&' || *src == '_' ||
		*src == '%' || *src == '#') {
		*p++ = '\\';
	    }
	    *p++ = *src++;
	}

	*p = '\0';
    }

    return targ;
}

/* Print the floating point number @x into the string @s, using the C
   format "%*.f", with the digits following the decimal point given by
   @dig.  This is intended for use with the LaTeX tabular column
   format "r@{.}l", so the decimal point is replaced by '&'.  In
   addition, it is presumed that we're _not_ in TeX math mode, so the
   leading minus sign, if present, is "mathized" as "$-$".

   NADBL is handled by printing a blank "r@{.}l" value.
*/

char *tex_rl_float (double x, char *s, int dig)
{
    char *p;

    if (na(x)) {
	return strcpy(s, "\\multicolumn{2}{c}{}");
    }

    x = screen_zero(x);

    if (x < 0) {
	sprintf(s, "$-$%.*f", dig, -x);
    } else {
	sprintf(s, "%.*f", dig, x);
    }

    p = strchr(s, '.');

    if (p == NULL) {
	p = strchr(s, ',');
    }

    if (p != NULL) {
	*p = '&';
    } else {
	strcat(s, "&");
    }

    return s;
}

/* When a floating point value is printed as TeX in the C format "%g",
   process the exponent part.  It is assumed we're not in TeX math
   mode, so in a negative exponent such as "e-06" the minus will
   appear as a hyphen, which does not look very good.  Using a true
   minus sign in the exponent doesn't look very good either (too
   wide); we compromise by writing the minus as an en dash.
*/

char *tex_modify_exponent (char *s)
{
    char *p = strchr(s, 'e');

    if (p != NULL) {
	int minus = (*(p+1) == '-');
	char tmp[16];

	sprintf(tmp, "\\textrm{e%s%s}", (minus)? "--" : "+", p + 2);
	strcpy(p, tmp);
    }

    return s;
}

static char *tex_rl_double_dig (double x, char *s, int d)
{
    char *p;

    if (na(x)) {
	return strcpy(s, "\\multicolumn{2}{c}{}");
    }

    x = screen_zero(x);

    if (x < 0) {
	sprintf(s, "$-$%#.*g", d, -x);
    } else {
	sprintf(s, "%#.*g", d, x);
    }

    if (strchr(s, 'e') != NULL) {
	tex_modify_exponent(s);
    }

    p = strchr(s, '.');

    if (p == NULL) {
	p = strchr(s, ',');
    }

    if (p != NULL) {
	*p = '&';
    } else {
	strcat(s, "&");
    }

    return s;
}

/* Print the floating point number @x into the string @s, using the C
   format "%#*.g", with GRETL_DIGITS of precision.  This is intended
   for use with the LaTeX tabular column format "r@{.}l", so the
   decimal point is replaced by '&'.  In addition, it is presumed that
   we're _not_ in TeX math mode, so the leading minus sign, if
   present, is "mathized" as "$-$".

   NADBL is handled by printing a blank "r@{.}l" value.
*/

char *tex_rl_double (double x, char *s)
{
    return tex_rl_double_dig(x, s, get_gretl_digits());
}

/* Print the floating point number @x into the string @s, using the C
   format "%#*.g", with @dig digits of precision.  It is presumed
   that we're _not_ in TeX math mode, so the leading minus sign, if
   present, is "mathized" as "$-$".

   NADBL is handled by printing a space.
*/

char *tex_sprint_double_digits (double x, char *s, int dig)
{
    if (na(x)) {
	return strcpy(s, " ");
    }

    x = screen_zero(x);

    if (x < 0.0) {
	sprintf(s, "$-$%#.*g", dig, -x);
    } else {
	sprintf(s, "%#.*g", dig, x);
    }

    if (strchr(s, 'e') != NULL) {
	tex_modify_exponent(s);
    }

    return s;
}

/* Basically as above, but for use when in TeX math mode */

static char *tex_sprint_math_double_digits (double x, char *s, int dig)
{
    if (na(x)) {
	return strcpy(s, "\\mbox{NA}");
    }

    x = screen_zero(x);

    sprintf(s, "%#.*g", dig, x);

    if (strchr(s, 'e') != NULL) {
	tex_modify_exponent(s);
    }

    return s;
}

/* Print the floating point number @x into the string @s, using the C
   format "%#*.g", with GRETL_DIGITS of precision.  It is presumed
   that we're _not_ in TeX math mode, so the leading minus sign, if
   present, is "mathized" as "$-$".

   NADBL is handled by printing a space.
*/

char *tex_sprint_double (double x, char *s)
{
    int d = get_gretl_digits();

    if (na(x)) {
	return strcpy(s, " ");
    }

    x = screen_zero(x);

    if (x < 0.0) {
	sprintf(s, "$-$%#.*g", d, -x);
    } else {
	sprintf(s, "%#.*g", d, x);
    }

    if (strchr(s, 'e') != NULL) {
	tex_modify_exponent(s);
    }

    return s;
}

/* Print the floating point number @x directly to @prn, using the C
   format "%#*.g", with GRETL_DIGITS of precision.  It is presumed
   that we're _not_ in TeX math mode, so the leading minus sign, if
   present, is "mathized" as "$-$".

   NADBL is handled by printing a space.
*/

void tex_print_double (double x, PRN *prn)
{
    char s[32];

    tex_sprint_double(x, s);
    pputs(prn, s);
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
    const char *tgreek = tex_greek_var(src);

    if (tgreek != NULL) {
	sprintf(targ, "$%s$", tgreek);
	return 1;
    } else {
	*targ = '\0';
	return 0;
    }
}

static void tex_garch_coeff_name (char *targ, const char *src,
				  int inmath)
{
    char fmt[16], vname[VNAMELEN], vnesc[16];
    int lag;

    sprintf(fmt, "%%%d[^(](%%d)", VNAMELEN - 1);

    if (sscanf(src, fmt, vname, &lag) == 2) {
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
    char fmt[12], vname[VNAMELEN], vnesc[24];
    int power;

    tex_escape(vnesc, src);

    sprintf(fmt, "%%%d[^^]^%%d", VNAMELEN - 1);

    if (sscanf(vnesc, fmt, vname, &power) == 2) {
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
    char vname[VNAMELEN], vnesc[32], texname[48];
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
    } else if (strstr(src, "(-") != NULL) {
	char fmt[16];

	sprintf(fmt, "%%%d[^(](-%%d)", VNAMELEN - 1);

	if (sscanf(src, fmt, vname, &i) == 2) {
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
    } else {
	tex_escape(vnesc, src);
	strcpy(targ, vnesc);
    }
}

static void tex_VAR_varname (char *s, const MODEL *pmod,
			     const DATASET *dset, int v)
{
    char tmp[32], base[12];
    int lag;

    gretl_model_get_param_name(pmod, dset, v, tmp);

    if (sscanf(tmp, "%11[^_]_%d", base, &lag) == 2) {
	sprintf(s, "%s$_{t-%d}$", base, lag);
    } else {
	tex_escape(s, tmp);
    }
}

static void tex_VECM_varname (char *s, const MODEL *pmod,
			      const DATASET *dset, int v)
{
    char tmp[32], base[12];
    int lag;

    gretl_model_get_param_name(pmod, dset, v, tmp);

    if (sscanf(tmp, "d_%11[^_]_%d", base, &lag) == 2) {
	sprintf(s, "$\\Delta$%s$_{t-%d}$", base, lag);
    } else {
	tex_escape(s, tmp);
    }
}

static int tex_print_coeff_custom (const model_coeff *mc, PRN *prn)
{
    char fmt[12];

    pprintf(prn, "%s & ", mc->name);

    if (colspec[0][0]) {
	/* coefficient */
	if (na(mc->b)) {
	    pprintf(prn, "\\multicolumn{1}{c}{\\rm %s}", _("undefined"));
	} else {
	    sprintf(fmt, "$%s$", colspec[0]);
	    pprintf(prn, fmt, mc->b);
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
	if (na(mc->se)) {
	    pprintf(prn, "\\multicolumn{1}{c}{\\rm %s}", _("undefined"));
	} else {
	    pprintf(prn, colspec[1], mc->se);
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
	if (na(mc->tval)) {
	    if (mc->show_tval) {
		pprintf(prn, "\\multicolumn{1}{c}{\\rm %s}", _("undefined"));
	    } else {
		pprintf(prn, "\\multicolumn{1}{c}{}");
	    }
	} else {
	    sprintf(fmt, "$%s$", colspec[2]);
	    pprintf(prn, fmt, mc->tval);
	}
    }

    if (colspec[3][0]) {
	if (colspec[0][0] || colspec[1][0] || colspec[2][0]) {
	    pputs(prn, " & ");
	}
	/* p-value (or perhaps slope) */
	if (mc->show_pval) {
	    if (na(mc->pval)) {
		pprintf(prn, "\\multicolumn{1}{c}{\\rm %s}", _("undefined"));
	    } else {
		pprintf(prn, colspec[3], mc->pval);
	    }
	} else if (!na(mc->slope)) {
	    pprintf(prn, colspec[3], mc->slope);
	} else {
	    pprintf(prn, "\\multicolumn{1}{c}{}");
	}
    }

    pputs(prn, " \\\\\n");

    return 0;
}

void make_tex_coeff_name (const MODEL *pmod, const DATASET *dset, int i,
			  char *name)
{
    char pname[VNAMELEN];

    gretl_model_get_param_name(pmod, dset, i, pname);

    if (pmod->aux == AUX_ARCH) {
	tex_make_cname(name, pname);
    } else if (pmod->ci == NLS) {
	if (!tex_greek_param(name, pname)) {
	    tex_escape(name, pname);
	}
    } else if (pmod->ci == ARMA) {
	tex_arma_coeff_name(name, pname, 0);
    } else if (pmod->ci == GARCH) {
	tex_garch_coeff_name(name, pname, 0);
    } else if (pmod->ci == VAR) {
	tex_VAR_varname(name, pmod, dset, i);
    } else if (pmod->aux == AUX_VECM) {
	tex_VECM_varname(name, pmod, dset, i);
    } else if (pmod->ci == MPOLS) {
	tex_mp_coeff_name(name, pname, 0);
    } else {
	tex_escape(name, pname);
    }
}

static char *tex_multi_double (double x, char *numstr)
{
    char *p;

    if (na(x)) {
	strcpy(numstr, " ");
    } else if (x < 0) {
	sprintf(numstr, "$-$%.15E", -x);
    } else {
	sprintf(numstr, "%.15E", x);
    }

    if ((p = strstr(numstr, "E-")) != NULL) {
	char tmp[8];

	sprintf(tmp, "E--%s", p + 2);
	strcpy(p, tmp);
    }

    return numstr;
}

void tex_print_coeff (const model_coeff *mc, PRN *prn)
{
    char col1[64], col2[64], col3[64], col4[64];
    int ncols = 4;

    if (mc->multi) {
	tex_multi_double(mc->b, col1);
	tex_multi_double(mc->se, col2);
	pprintf(prn, "%s & %s & %s \\\\\n", mc->name, col1, col2);
	return;
    }

    if (use_custom) {
	tex_print_coeff_custom(mc, prn);
	return;
    }

    if (na(mc->b)) {
	sprintf(col1, "\\multicolumn{2}{c}{\\rm %s}", _("undefined"));
    } else {
	tex_rl_double(mc->b, col1);
    }

    if (!na(mc->lo) && !na(mc->hi)) {
	tex_rl_double(mc->lo, col2);
	tex_rl_double(mc->hi, col3);
	ncols = 3;
    } else {
	if (na(mc->se)) {
	    sprintf(col2, "\\multicolumn{2}{c}{\\rm %s}", _("undefined"));
	} else {
	    tex_rl_double(mc->se, col2);
	}

	if (na(mc->tval)) {
	    if (mc->show_tval) {
		sprintf(col3, "\\multicolumn{2}{c}{\\rm %s}", _("undefined"));
	    } else {
		strcpy(col3, "\\multicolumn{2}{c}{}");
	    }
	} else {
	    /* note: was tex_rl_float(mc->tval, col3, 4) */
	    tex_rl_double_dig(mc->tval, col3, 4);
	}
    }

    *col4 = '\0';

    if (!mc->show_pval && na(mc->slope)) {
	strcpy(col4, "\\multicolumn{2}{c}{}");
    } else if (!na(mc->slope)) {
	tex_rl_double(mc->slope, col4);
    } else if (mc->show_pval) {
	if (!na(mc->pval)) {
	    tex_rl_float(mc->pval, col4, 4);
	}
    }

    pprintf(prn, "%s &\n"
	    "  %s &\n"
	    "    %s &\n",
	    mc->name,
	    col1,
	    col2);

    if (ncols == 4) {
	pprintf(prn,
		"      %s &\n"
		"        %s \\\\\n",
		col3,
		col4);
    } else {
	pprintf(prn,
		"      %s \\\\\n",
		col3);
    }
}

static int
tex_custom_coeff_table_start (const char **cols, gretlopt opt, PRN *prn)
{
    int i, ncols = 0;

    for (i=0; i<4; i++) {
	if (colspec[i][0]) {
	    ncols++;
	}
    }

    if (!(opt & OPT_U)) {
	/* not a user-defined model */
	pputs(prn, "\\vspace{1em}\n\n");
    }

    pputs(prn, "\\begin{tabular}{l");

    for (i=0; i<ncols; i++) {
	pputc(prn, 'r');
    }

    pputs(prn, "}\n");

    pprintf(prn, "\\multicolumn{1}{c}{%s} &\n", _(cols[0]));

    if (colspec[0][0]) {
	pprintf(prn, "\\multicolumn{1}{c}{%s}", _(cols[1]));
    }

    if (!colspec[1][0] && !colspec[2][0] && !colspec[3][0]) {
	pputs(prn, " \\\\\n");
	return ncols;
    }

    if (colspec[1][0]) {
	if (colspec[0][0]) {
	    pputs(prn, " &\n");
	}
	pprintf(prn, "\\multicolumn{1}{c}{%s}", _(cols[2]));
    }

    if (!colspec[2][0] && !colspec[3][0]) {
	pputs(prn, " \\\\\n");
	return ncols;
    }

    if (colspec[2][0]) {
	if (colspec[0][0] || colspec[1][0]) {
	    pputs(prn, " &\n");
	}
	pprintf(prn, "\\multicolumn{1}{c}{%s}", _(cols[3]));
    }

    if (colspec[3][0]) {
	if (colspec[0][0] || colspec[1][0] || colspec[2][0]) {
	    pputs(prn, " &\n");
	}
	pprintf(prn, "\\multicolumn{1}{c}{%s}", _(cols[4]));
    }

    pputs(prn, " \\\\\n");

    return ncols;
}

/* returns the number of columns in the coeff table */

int tex_coeff_table_start (const char **cols, gretlopt opt, PRN *prn)
{
    char pt = get_local_decpoint();
    int i, mcols, binary = (opt & OPT_B);
    int ncols = 1;

    if (use_custom) {
	return tex_custom_coeff_table_start(cols, opt, prn);
    }

    if (!(opt & OPT_U)) {
	/* not a user-defined model */
	pputs(prn, "\\vspace{1em}\n\n");
    }

    pputs(prn, "\\begin{tabular}{l");

    for (i=1; cols[i] != NULL; i++) {
	if (opt & OPT_M) {
	    pputc(prn, 'r');
	} else {
	    pprintf(prn, "r@{%c}l", pt);
	}
	ncols += 2;
    }

    pprintf(prn, "}\n%s &\n", _(cols[0]));

    mcols = (opt & OPT_M)? 1 : 2;

    for (i=1; cols[i] != NULL; i++) {
	bufspace(i, prn);
	pprintf(prn, "\\multicolumn{%d}{c}{%s%s} %s\n", mcols, _(cols[i]),
		(cols[i+1] == NULL && binary)? "$^*$" : "",
		(cols[i+1] == NULL)? "\\\\[1ex]" : "&");
    }

    return ncols;
}

void tex_coeff_table_end (PRN *prn)
{
    pputs(prn, "\\end{tabular}\n\n");
}

void tex_print_VECM_omega (GRETL_VAR *vecm, const DATASET *dset, PRN *prn)
{
    char vname[48];
    const int *list = vecm->ylist;
    double x;
    int i, j;

    pprintf(prn, "%s\n\n", _("Cross-equation covariance matrix"));
    pputs(prn, "\\vspace{1em}\n");

    pputs(prn, "\\begin{tabular}{");
    pputs(prn, "l");
    for (i=0; i<vecm->neqns; i++) {
	pputs(prn, "r");
    }
    pputs(prn, "}\n & ");

    for (i=0; i<vecm->neqns; i++) {
	tex_escape(vname, dset->varname[list[i+1]]);
	pprintf(prn, "$\\Delta$%s ", vname);
	if (i == vecm->neqns - 1) {
	    pputs(prn, "\\\\\n");
	} else {
	    pputs(prn, "& ");
	}
    }
    pputc(prn, '\n');

    for (i=0; i<vecm->neqns; i++) {
	tex_escape(vname, dset->varname[list[i+1]]);
	pprintf(prn, "$\\Delta$%s & ", vname);
	for (j=0; j<vecm->neqns; j++) {
	    x = gretl_matrix_get(vecm->S, i, j);
	    tex_print_double(x, prn);
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
    pprintf(prn, "%s = ", _("determinant"));
    tex_print_double(exp(vecm->ldet), prn);
    pputs(prn, "\\\\\n");
}

static void tex_beta_vname (char *s,
			    const GRETL_VAR *v,
			    const DATASET *dset,
			    int i, PRN *prn)
{
    if (i < v->neqns) {
	tex_escape(s, dset->varname[v->ylist[i+1]]);
	pprintf(prn, "%s$_{t-1}$ & ", s);
    } else if (auto_restr(v) && i == v->neqns) {
	pprintf(prn, "%s & ", (jcode(v) == J_REST_CONST)? "const" : "trend");
    } else if (v->rlist != NULL) {
	int k = i - v->ylist[0] - auto_restr(v) + 1;

	tex_escape(s, dset->varname[v->rlist[k]]);
	pprintf(prn, "%s$_{t-1}$ & ", s);
    }
}

void tex_print_VECM_coint_eqns (GRETL_VAR *vecm, const DATASET *dset, PRN *prn)
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
    } else {
	pputc(prn, '\n');
    }

    pputs(prn, "\n\\vspace{1em}\n");

    pputs(prn, "\\begin{tabular}{");
    pputs(prn, "l");
    for (i=0; i<jv->rank; i++) {
	pputs(prn, "r");
    }
    pputs(prn, "}\n");

    for (i=0; i<rows; i++) {
	tex_beta_vname(s, vecm, dset, i, prn);

	/* coefficients */
	for (j=0; j<jv->rank; j++) {
	    x = gretl_matrix_get(jv->Beta, i, j);
	    if (jv->Bse == NULL) {
		x /= gretl_matrix_get(jv->Beta, j, j);
	    }
	    tex_print_double(x, prn);
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
		x = gretl_matrix_get(jv->Bse, i, j);
		pputc(prn, '(');
		tex_print_double(x, prn);
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
    pputc(prn, '\n');

    rows = gretl_matrix_rows(jv->Alpha);

    pputs(prn, "\\noindent\n");
    pprintf(prn, _("Adjustment vectors"));
    if (jv->Ase != NULL) {
	pprintf(prn, " (%s)\n", _("standard errors in parentheses"));
    } else {
	pputc(prn, '\n');
    }

    pputs(prn, "\n\\vspace{1em}\n");

    pputs(prn, "\\begin{tabular}{");
    pputs(prn, "l");
    for (i=0; i<jv->rank; i++) {
	pputs(prn, "r");
    }
    pputs(prn, "}\n");

    for (i=0; i<rows; i++) {
	tex_beta_vname(s, vecm, dset, i, prn);

	/* coefficients */
	for (j=0; j<jv->rank; j++) {
	    x = gretl_matrix_get(jv->Alpha, i, j);
	    if (jv->Ase == NULL) {
		x /= gretl_matrix_get(jv->Alpha, j, j);
	    }
	    tex_print_double(x, prn);
	    if (j == jv->rank - 1) {
		pputs(prn, "\\\\\n");
	    } else {
		pputs(prn, "& ");
	    }
	}

	if (jv->Ase != NULL) {
	    /* standard errors */
	    pputs(prn, " & ");
	    for (j=0; j<jv->rank; j++) {
		x = gretl_matrix_get(jv->Ase, i, j);
		pputc(prn, '(');
		tex_print_double(x, prn);
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
    pputc(prn, '\n');

}

void tex_print_VAR_ll_stats (GRETL_VAR *var, PRN *prn)
{
    pprintf(prn, "\\noindent\n%s = ", _("Log-likelihood"));
    tex_print_double(var->ll, prn);
    pputs(prn, "\\par\n");

    pprintf(prn, "\\noindent\n%s = ", _("Determinant of covariance matrix"));
    tex_print_double(exp(var->ldet), prn);
    pputs(prn, "\\par\n");

    pprintf(prn, "\\noindent\n%s $= %.4f$ \\par\n", _("AIC"), var->AIC);
    pprintf(prn, "\\noindent\n%s $= %.4f$ \\par\n", _("BIC"), var->BIC);
    pprintf(prn, "\\noindent\n%s $= %.4f$ \\par\n", _("HQC"), var->HQC);
}

static PRN *make_tex_prn (const char *fname,
			  int eqn, int doc,
			  int *err)
{
    PrnFormat fmt = GRETL_FORMAT_TEX;
    PRN *prn;

    gretl_maybe_switch_dir(fname);
    prn = gretl_print_new_with_filename(fname, err);

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

static const char *set_tex_locale_filename (char *local)
{
    char *lang = getenv("LANG");

    *local = '\0';

    if (lang != NULL) {
	char lstr[3] = {0};

	strncat(lstr, lang, 2);
	sprintf(local, "gretlpre_%s.tex", lstr);
    }

    return local;
}

static int find_gretlpre (const char *path, const char *localname)
{
    char test[MAXLEN];
    int err, gotit = 0;

    if (*localname != '\0') {
	/* localized preamble file? */
	sprintf(test, "%s%s", path, localname);
	err = gretl_test_fopen(test, "r");
	if (!err) {
	    strcpy(tex_preamble_file, test);
	    gotit = 1;
	}
    }

    if (!gotit) {
	/* regular preamble file? */
	sprintf(test, "%sgretlpre.tex", path);
	err = gretl_test_fopen(test, "r");
	if (!err) {
	    strcpy(tex_preamble_file, test);
	    gotit = 1;
	}
    }

    return gotit;
}

void set_gretl_tex_preamble (void)
{
    const char *path = gretl_workdir();
    char localname[16];
    int gotit;

    set_tex_locale_filename(localname);
    gotit = find_gretlpre(path, localname);

    if (!gotit) {
	path = maybe_get_default_workdir();
	if (path != NULL) {
	    gotit = find_gretlpre(path, localname);
	}
    }

    if (!gotit) {
#ifdef __APPLE__
	path = gretl_app_support_dir();
#else
	path = gretl_dotdir();
#endif
	gotit = find_gretlpre(path, localname);
    }

    gretl_error_clear();
}

static void landscape_modify_line (char *line)
{
    char *p, *rem;

    if (strstr(line, "landscape")) {
	return;
    }

    p = strstr(line, "documentclass");

    if (p != NULL) {
	if (*(p + 13) == '[') {
	    p = strchr(p, ']');
	    if (p != NULL) {
		rem = gretl_strdup(p);
		if (rem != NULL) {
		    sprintf(p, ",landscape%s", rem);
		    free(rem);
		}
	    }
	} else {
	    p += 13;
	    rem = gretl_strdup(p);
	    if (rem != NULL) {
		sprintf(p, "[landscape]%s", rem);
		free(rem);
	    }
	}
    }
}

void gretl_tex_preamble (PRN *prn, int fmt)
{
    char *lang = getenv("LANG");
    FILE *fp = NULL;
    int userfile = 0;

    if (*tex_preamble_file != '\0') {
	fp = gretl_fopen(tex_preamble_file, "r");
	if (fp != NULL) {
	    char line[256];

	    /* FIXME model table: longtable and geom packages */

	    while (fgets(line, sizeof line, fp)) {
		if (strstr(line, "documentclass") &&
		    (fmt & GRETL_FORMAT_LANDSCAPE)) {
		    landscape_modify_line(line);
		}
		pputs(prn, line);
	    }
	    userfile = 1;
	    fclose(fp);
	    fprintf(stderr, "gretltex: using preamble file\n %s\n",
		    tex_preamble_file);
	}
    }

    if (!userfile) {
	const char *paper = in_usa()? "letterpaper" : "a4paper";
	const char *driver = use_pdf ? "pdftex" : "dvips";
	const char *margin = "";

	if (fmt & GRETL_FORMAT_MODELTAB) {
	    margin = "margin=2cm,";
	}

	pputs(prn, "\\documentclass");

	if (fmt & GRETL_FORMAT_MODELTAB) {
	    if (fmt & GRETL_FORMAT_LANDSCAPE) {
		pputs(prn, "[landscape]");
	    }
	} else if (fmt & GRETL_FORMAT_LANDSCAPE) {
	    pputs(prn, "[11pt,landscape]");
	} else {
	    pputs(prn, "[11pt]");
	}

	pputs(prn, "{article}\n");

#ifdef ENABLE_NLS
	pputs(prn, "\\usepackage[utf8]{inputenc}\n");
#endif

	if (lang != NULL && !strncmp(lang, "ru", 2)) {
	    pputs(prn, "\\usepackage[russian]{babel}\n");
	}

	pprintf(prn, "\\usepackage[%s,%s%s]{geometry}\n", paper, margin,
		driver);

	if (fmt & GRETL_FORMAT_EQN) {
	    pputs(prn, "\\usepackage{amsmath}\n");
	} else {
	    pputs(prn, "\\usepackage{longtable}\n");
	}

	pputs(prn, "\n\\begin{document}\n\n"
	      "\\thispagestyle{empty}\n\n");
    }
}

/* For use when printing a model in equation style: print the value
   unsigned, since the sign will be handled separately.
*/

static void tex_print_unsigned_double (double x, PRN *prn)
{
    tex_print_double(fabs(x), prn);
}

#define MAXCOEFF 4

/**
 * tex_print_equation:
 * @pmod:  pointer to gretl MODEL struct.
 * @dset:  information regarding the data set.
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

int tex_print_equation (const MODEL *pmod, const DATASET *dset,
			gretlopt opt, PRN *prn)
{
    double x;
    char tmp[48], vname[32];
    int i, nc = pmod->ncoeff;
    int split = 0, offvar = 0;
    int cchars = 0, ccount = 0;
    int se_digits;
    int sderr_ok = 1;

    if (pmod->ci == HECKIT) {
	return E_NOTIMP;
    }

    if (COUNT_MODEL(pmod->ci)) {
	offvar = gretl_model_get_int(pmod, "offset_var");
	if (offvar > 0) {
	    nc++;
	}
    }

    split = (nc > MAXCOEFF);

    if (opt & OPT_S) {
	gretl_tex_preamble(prn, GRETL_FORMAT_EQN);
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
    } else {
	i = gretl_model_get_depvar(pmod);
	tex_escape(tmp, dset->varname[i]);
    }

    if (0 /* FIXME should this apply for DPANEL? (was ARBOND-specific) */) {
	pprintf(prn, "\\widehat{\\Delta \\rm %s} %s= \n", tmp, (split? "&" : ""));
    } else {
	pprintf(prn, "\\widehat{\\rm %s} %s= \n", tmp, (split? "&" : ""));
    }

    if (pmod->ci == GARCH) {
	nc -= (1 + pmod->list[1] + pmod->list[2]);
    } else if (pmod->ci == PANEL) {
	nc = pmod->list[0] - 1;
    }

    se_digits = get_gretl_digits();
    se_digits = se_digits > 5 ? 5 : se_digits;

    /* coefficients times indep vars */
    for (i=0; i<nc; i++) {
	if (offvar > 0 && i == nc - 1) {
	    pputc(prn, '+');
	    tex_print_double(1.0, prn);
	} else {
	    if (na(pmod->sderr[i])) {
		sderr_ok = 0;
		pprintf(prn, "%s{", (pmod->coeff[i] < 0.0)? "-" :
			(i > 0)? "+" : "");
	    } else if (sderr_ok) {
		if (opt & OPT_T) {
		    /* t-ratios */
		    x = pmod->coeff[i] / pmod->sderr[i];
		    pprintf(prn, "%s\\underset{(%.3f)}{",
			    (pmod->coeff[i] < 0.0)? "-" :
			    (i > 0)? "+" : "", x);
		} else {
		    /* standard errors */
		    tex_sprint_math_double_digits(pmod->sderr[i], tmp, se_digits);
		    pprintf(prn, "%s\\underset{(%s)}{",
			    (pmod->coeff[i] < 0.0)? "-" :
			    (i > 0)? "+" : "", tmp);
		}
	    }
	    tex_print_unsigned_double(pmod->coeff[i], prn);
	    pputc(prn, '}');
	}

	if (i > 0 || pmod->ifc == 0) {
	    /* regular coefficient, not const */
	    if (offvar > 0 && i == nc - 1) {
		strcpy(vname, dset->varname[offvar]);
	    } else {
		gretl_model_get_param_name(pmod, dset, i, vname);
	    }
	    cchars += strlen(vname);

	    pputs(prn, "\\,");

	    if (pmod->ci == ARMA) {
		tex_arma_coeff_name(tmp, vname, 1);
		pputs(prn, tmp);
	    } else if (pmod->ci == GARCH) {
		tex_garch_coeff_name(tmp, vname, 1);
		pputs(prn, tmp);
	    } else if (pmod->ci == MPOLS) {
		tex_mp_coeff_name(tmp, vname, 1);
		pputs(prn, tmp);
	    } else {
		tex_escape(tmp, vname);
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
	    tex_sprint_math_double_digits(pmod->sderr[r], tmp, se_digits);
	    pprintf(prn, "\\hat{\\sigma}^2_t = \\underset{(%s)}{%g} ", /* FIXME? */
		    tmp, pmod->coeff[r]);
	}

	for (i=1; i<=q; i++) {
	    if (opt & OPT_T) {
		x = pmod->coeff[r+i] / pmod->sderr[r+i];
		pprintf(prn, "%s\\underset{(%.3f)}{",
			(pmod->coeff[r+i] < 0.0)? "-" : "+", x);
	    } else {
		tex_sprint_math_double_digits(pmod->sderr[r+i], tmp, se_digits);
		pprintf(prn, "%s\\underset{(%s)}{",
			(pmod->coeff[r+i] < 0.0)? "-" : "+", tmp);
	    }
	    tex_print_unsigned_double(pmod->coeff[r+i], prn);
	    pputs(prn, "}\\,");
	    pprintf(prn, "\\varepsilon^2_{t-%d}", i);
	}

	for (i=1; i<=p; i++) {
	    if (opt & OPT_T) {
		x = pmod->coeff[q+r+i] / pmod->sderr[q+r+i];
		pprintf(prn, "%s\\underset{(%.3f)}{",
			(pmod->coeff[q+r+i] < 0.0)? "-" : "+", x);
	    } else {
		tex_sprint_math_double_digits(pmod->sderr[q+r+i], tmp, se_digits);
		pprintf(prn, "%s\\underset{(%s)}{",
			(pmod->coeff[q+r+i] < 0.0)? "-" : "+", tmp);
	    }
	    tex_print_unsigned_double(pmod->coeff[q+r+i], prn);
	    pputs(prn, "}\\,");
	    pprintf(prn, "\\sigma^2_{t-%d}", i);
	}

	pputs(prn, "\\notag \\\\\n");
    }

    pprintf(prn, "T = %d ", pmod->nobs);

    /* additional info (R^2 etc) */
    if (pmod->ci == LAD) {
	x = gretl_model_get_double(pmod, "ladsum");
	if (!na(x)) {
	    tex_sprint_math_double_digits(x, tmp, 5);
	    pprintf(prn, "\\quad \\sum |\\hat{u}_t| = %s ", tmp);
	}
    } else {
	if (!na(pmod->adjrsq)) {
	    pprintf(prn, "\\quad \\bar{R}^2 = %.4f ", pmod->adjrsq);
	} else if (!na(pmod->lnL)) {
	    pprintf(prn, "\\quad \\mbox{ln}L = %.4f ", pmod->lnL);
	}

	if (pmod->ci != LOGIT && pmod->ci != PROBIT && !na(pmod->fstt)) {
	    tex_sprint_math_double_digits(pmod->fstt, tmp, 5);
	    pprintf(prn, "\\quad F(%d,%d) = %s ",
		    pmod->dfn, pmod->dfd, tmp);
	}

	if (!na(pmod->sigma)) {
	    tex_sprint_math_double_digits(pmod->sigma, tmp, 5);
	    pprintf(prn, "\\quad \\hat{\\sigma} = %s ", tmp);
	}

	if (!na(gretl_model_get_double(pmod, "rho_gls"))) {
	    x = gretl_model_get_double(pmod, "rho_gls");
	    tex_sprint_math_double_digits(x, tmp, 5);
	    pprintf(prn, " \\quad \\rho = %s", tmp);
	}
    }

    pputs(prn, "\\notag \\\\\n");

    if (sderr_ok) {
	pprintf(prn, "\\centerline{(%s)} \\notag\n",
		(opt & OPT_T)? _("$t$-statistics in parentheses") :
		_("standard errors in parentheses"));
    } else {
	pputs(prn, "\\notag\n");
    }

    pputs(prn, "\\end{gather}\n");

    if (opt & OPT_S) {
	pputs(prn, "\n\\end{document}\n");
    }

    return 0;
}

/**
 * tex_print_model:
 * @pmod:  pointer to gretl MODEL struct.
 * @dset: information regarding the data set.
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

int tex_print_model (MODEL *pmod, const DATASET *dset,
		     gretlopt opt, PRN *prn)
{
    int err = 0;

    if (RQ_SPECIAL_MODEL(pmod)) {
	return E_NOTIMP;
    }

    if (tex_doc_format(prn)) {
	opt |= OPT_S;
    }

    if (tex_eqn_format(prn)) {
	err = tex_print_equation(pmod, dset, opt, prn);
    } else {
	if (opt & OPT_T) {
	    /* --format option */
	    const char *s = get_optval_string(TABPRINT, OPT_T);

	    err = set_tex_param_format(s);
	}
	if (!err) {
	    err = printmodel(pmod, dset, OPT_NONE, prn);
	}
    }

    return err;
}

/**
 * texprint:
 * @pmod: pointer to model.
 * @dset: information regarding the data set.
 * @fname: name of file to save.
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

int texprint (MODEL *pmod, const DATASET *dset,
	      const char *fname, gretlopt opt)
{
    PRN *prn;
    int eqn = (opt & OPT_E);
    int doc = (opt & OPT_O);
    int err = 0;

    if (RQ_SPECIAL_MODEL(pmod)) {
	return E_NOTIMP;
    }

    prn = make_tex_prn(fname, eqn, doc, &err);

    if (!err) {
	err = tex_print_model(pmod, dset, opt, prn);
	gretl_print_destroy(prn);
    }

    return err;
}

static void out_crlf (const char *buf, FILE *fp)
{
    const char *p = buf;

#ifdef G_OS_WIN32
    /* fputs on Windows should take care of CR, LF */
    fputs(p, fp);
#else
    while (*p) {
	if (*p == '\n') {
	    fputs("\r\n", fp);
	} else {
	    fputc(*p, fp);
	}
	p++;
    }
#endif
}

int rtfprint (MODEL *pmod, const DATASET *dset,
	      const char *fname, gretlopt opt)
{
    const char *buf = NULL;
    char *trbuf = NULL;
    PRN *prn;
    int err = 0;

    if (RQ_SPECIAL_MODEL(pmod)) {
	return E_NOTIMP;
    }

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);

    if (!err) {
	/* print the model to buffer first */
	gretl_print_set_format(prn, GRETL_FORMAT_RTF);
	err = printmodel(pmod, dset, opt, prn);
    }

    if (!err) {
	/* recode if necessary */
	buf = gretl_print_get_buffer(prn);
	if (!gretl_is_ascii(buf)) {
	    trbuf = utf8_to_rtf(buf);
	    if (trbuf == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (!err) {
	/* now send to file, converting LF to CR + LF
	   as we go, if required
	*/
	FILE *fp;

	gretl_maybe_switch_dir(fname);
	fp = gretl_fopen(fname, "w");

	if (fp == NULL) {
	    err = E_FOPEN;
	} else if (trbuf != NULL) {
	    out_crlf(trbuf, fp);
	} else {
	    out_crlf(buf, fp);
	}

	if (fp != NULL) {
	    fclose(fp);
	}
    }

    if (trbuf != NULL) {
	free(trbuf);
    }

    if (prn != NULL) {
	gretl_print_destroy(prn);
    }

    return err;
}

static void out_native (const char *buf, FILE *fp)
{
#ifdef G_OS_WIN32
    if (!gretl_is_ascii(buf)) {
	/* Windows: if the text is UTF-8, prepend
	   the UTF-8 BOM */
	fputc(0xEF, fp);
	fputc(0xBB, fp);
	fputc(0xBF, fp);
    }
#endif
    fputs(buf, fp);
}

int csvprint (MODEL *pmod, const DATASET *dset,
	      const char *fname, gretlopt opt)
{
    PRN *prn;
    int err = 0;

    if (RQ_SPECIAL_MODEL(pmod)) {
	return E_NOTIMP;
    }

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);

    if (!err) {
	/* print to buffer first */
	gretl_print_set_format(prn, GRETL_FORMAT_CSV);
	err = printmodel(pmod, dset, opt, prn);
    }

    if (!err) {
	/* then send to file */
	const char *buf = gretl_print_get_buffer(prn);
	FILE *fp;

	gretl_maybe_switch_dir(fname);
	fp = gretl_fopen(fname, "w");

	if (fp == NULL) {
	    err = E_FOPEN;
	} else {
	    out_native(buf, fp);
	}

	if (fp != NULL) {
	    fclose(fp);
	}
    }

    if (prn != NULL) {
	gretl_print_destroy(prn);
    }

    return err;
}

/**
 * tex_print_obs_marker:
 * @t: observation number.
 * @dset: data information struct.
 * @prn: gretl printing struct.
 *
 * Print a string (label, date or obs number) representing the given @t.
 */

void tex_print_obs_marker (int t, const DATASET *dset, PRN *prn)
{
    if (dset->markers) {
	pprintf(prn, "\\texttt{%s} ", dset->S[t]);
    } else {
	char tmp[OBSLEN];

	ntolabel(tmp, t, dset);
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
 *
 * Returns: 0 on success, non-zero code on error.
 */

int set_tex_param_format (const char *s)
{
    const char *p = s;
    int i, n = 0;
    int err = 0;

    if (s == NULL || !strcmp(s, "default")) {
	use_custom = 0;
	return 0;
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
#if 0
	    fprintf(stderr, "spec %d = '%s'\n", i, colspec[i]);
#endif
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
	    err = E_ARGS;
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

    return err;
}
