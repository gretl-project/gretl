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

/* support for C-like (s)printf and sscanf */

#include "libgretl.h"
#include "libset.h"
#include "gretl_func.h"
#include "gretl_string_table.h"
#include "gretl_scalar.h"
#include "matrix_extra.h"

#include <errno.h>

#define PSDEBUG 0

/* For the point of "x *= (1 + 0x1p-52)" below, see
   http://sourceware.org/bugzilla/show_bug.cgi?id=4943 ,
   contribution from Eric Postpischil, 2007-10-05.
*/

static void printf_series (int v, const double **Z,
			   const DATAINFO *pdinfo,
			   const char *fmt, 
			   int wid, int prec, 
			   int wstar, int pstar,
			   PRN *prn)
{
    char label[OBSLEN];
    double x;
    int n, t;

    n = max_obs_label_length(pdinfo) + 1;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	get_obs_string(label, t, pdinfo);
	pprintf(prn, "%*s ", n, label);
	x = Z[v][t];
	if (na(x)) {
	    pputc(prn, '\n');
	    continue;
	}
	x *= (1 + 0x1p-52);
	if (wstar && pstar) {
	    pprintf(prn, fmt, wid, prec, x);
	} else if (wstar || pstar) {
	    wid = (wstar)? wid : prec;
	    pprintf(prn, fmt, wid, x);
	} else {
	    pprintf(prn, fmt, x);
	}
	pputc(prn, '\n');
    }
}

static int printf_escape (int c, PRN *prn)
{
    switch (c) {
    case 'n':
	pputc(prn, '\n');
	break;
    case 't':
	pputc(prn, '\t');
	break;
    case 'v':
	pputc(prn, '\v');
	break;
    case '\\':
	pputc(prn, '\\');
	break;
    case '"':
	pputc(prn, '"');
	break;
    default:
	/* treat as literal backslash */
	pputc(prn, '\\');
	pputc(prn, c);
    }

    return 0;
}

static char *printf_get_string (char *s, double ***pZ,
				DATAINFO *pdinfo, int t, 
				int *err)
{
    const char *p = NULL;
    const char *q = NULL;
    char *ret = NULL;
    int len = 0;

    /* special, not yet in genr */

    if (!strncmp(s, "marker", 6) && pdinfo->S != NULL) {
	/* observation label */
	int offset = 0;

	p = pdinfo->S[t];
	len = strlen(p);
	q = s + 6;
	while (isspace(*q)) q++;
	if (*q == '+') {
	    q++;
	    offset = atoi(q);
	    len -= offset;
	}
	if (len >= 0) {
	    ret = gretl_strndup(p + offset, len);
	}
    } else {
	ret = generate_string(s, pZ, pdinfo, err);
    }

    if (ret == NULL) {
	ret = gretl_strdup("");
    }

    if (ret == NULL) {
	*err = E_ALLOC;
    }

    return ret;
}

static double printf_get_scalar (char *s, double ***pZ,
				 DATAINFO *pdinfo, int t, 
				 int *err)
{
    double x = NADBL;
    int v;

#if PSDEBUG
    fprintf(stderr, "printf_get_scalar: looking at '%s'\n", s);
#endif

    if (numeric_string(s)) {
	x = dot_atof(s);
    } else if (gretl_is_scalar(s)) {
	x = gretl_scalar_get_value(s);
    } else {
	v = series_index(pdinfo, s);

	if (v < pdinfo->v) {
	    char genstr[32];

	    sprintf(genstr, "%s[%d]", s, t + 1);
	    x = generate_scalar(genstr, pZ, pdinfo, err);
	} else {
	    x = generate_scalar(s, pZ, pdinfo, err);
	}
    }

#if PSDEBUG
    fprintf(stderr, "printf_get_scalar: returning %g\n", x);
#endif

    return x;
}

/* dup argv (s) up to the next free comma */

static char *get_next_arg (const char *s, int *len, int *err)
{
    const char *p;
    char *arg = NULL;
    int par = 0, br = 0;
    int quoted = 0;
    int n = 0;

    if (s == NULL) {
	*err = E_PARSE;
	return NULL;
    }

    *len = strspn(s, ", ");
    s += *len;
    p = s;
    
    while (*p) {
	if (*p == '"' && (p == s || *(p-1) != '\\')) {
	    quoted = !quoted;
	}
	if (!quoted) {
	    if (*p == '(') par++;
	    else if (*p == ')') par--;
	    else if (*p == '[') br++;
	    else if (*p == ']') br--;
	}
	if (!quoted && !par && !br && *p == ',') {
	    break;
	}
	p++;
    }

    n = p - s;
    *len += n;

    if (n > 0) {
	arg = gretl_strndup(s, n);
	if (arg == NULL) {
	    *err = E_ALLOC;
	}
    } else {
	*err = E_PARSE;
    }

#if PSDEBUG
    fprintf(stderr, "get_next_arg: got '%s'\n", arg);
#endif

    return arg;
}

/* dup format string up to the end of the current conversion:
   e.g. "%10.4f", "%6g", "%.8g", %3s", "%.*s", "%#12.4g"
*/

static char *
get_printf_format_chunk (const char *s, int *fc, 
			 int *len, int *wstar, int *pstar,
			 int *err)
{
    const char *cnvchars = "eEfgGduxs";
    const char *numchars = "0123456789";
    char *chunk = NULL;
    const char *p = s;
    int n;

    p++; /* move past % */

    /* '#', if present, must come first */
    if (*p == '#') {
	p++;
    }

    /* zero padding? */
    if (*p == '0') {
	p++;
    }

    /* left justification? */
    if (*p == '-') {
	p++;
    }

    /* leading space before positive number? */
    if (*p == ' ') {
	p++;
    }

    /* always print sign? */
    if (*p == '+') {
	p++;
    }
 
    /* optional width? */
    n = strspn(p, numchars);
    if (n == 0 && *p == '*') {
	/* variable version */
	*wstar = 1;
	p++;
    } else if (n > 3) {
	*err = E_PARSE;
	return NULL;
    } else {
	p += n;
    }

    /* optional dot plus precision? */
    if (*p == '.') {
	p++;
	n = strspn(p, numchars);
	if (n == 0 && *p == '*') {
	    *pstar = 1;
	    p++;
	} else if (n > 3) {
	    *err = E_PARSE;
	    return NULL;
	} else {
	    p += n;
	}
    }

    /* now we should have a conversion character */
    if (*p == '\0' || strchr(cnvchars, *p) == NULL) {
	fprintf(stderr, "bad conversion '%c'\n", *p);
	*err = E_PARSE;
	return NULL;
    }

    *fc = *p++;
    *len = n = p - s;
    
    if (n > 0) {
	chunk = gretl_strndup(s, n);
	if (chunk == NULL) {
	    *err = E_ALLOC;
	}
    } else {
	*err = E_PARSE;
    }

#if PSDEBUG
    fprintf(stderr, "get_printf_format_chunk: got '%s'\n", chunk);
#endif

    return chunk;
}

/* extract the next conversion spec from *pfmt, find and evaluate the
   corresponding elements in *pargs, and print the result */

static int print_arg (char **pfmt, char **pargs, 
		      double ***pZ, DATAINFO *pdinfo,
		      int t, PRN *prn)
{
    const char *intconv = "dxul";
    char *fmt = NULL;
    char *arg = NULL;
    char *str = NULL;
    gretl_matrix *m = NULL;
    int series_v = -1;
    double x = NADBL;
    int flen = 0, alen = 0;
    int wstar = 0, pstar = 0;
    int wid = 0, prec = 0;
    int fc = 0;
    int err = 0;

#if PSDEBUG
    fprintf(stderr, "print_arg: *pfmt='%s', *pargs='%s'\n",
	    *pfmt, *pargs);
#endif

    /* select current conversion format */
    fmt = get_printf_format_chunk(*pfmt, &fc, &flen, &wstar, &pstar, &err);
    if (err) {
	return err;
    }

    *pfmt += flen;

    if (wstar) {
	/* evaluate field width specifier */
	arg = get_next_arg(*pargs, &alen, &err);
	if (!err) {
	    x = printf_get_scalar(arg, pZ, pdinfo, t, &err);
	    if (!err && (na(x) || fabs(x) > 255)) {
		err = E_DATA;
	    }
	}
	free(arg);
	if (err) {
	    goto bailout;
	}
	*pargs += alen;
	wid = x;
    }

    if (pstar) {
	/* evaluate precision specifier */
	arg = get_next_arg(*pargs, &alen, &err);
	if (!err) {
	    x = printf_get_scalar(arg, pZ, pdinfo, t, &err);
	    if (!err && (na(x) || fabs(x) > 255)) {
		err = E_DATA;
	    }
	}
	free(arg);
	if (err) {
	    goto bailout;
	}
	*pargs += alen;
	prec = x;
    }

    /* get next substantive arg */
    arg = get_next_arg(*pargs, &alen, &err);
    if (!err) {
	int v;

	if (fc == 's') {
	    str = printf_get_string(arg, pZ, pdinfo, t, &err);
	} else if ((m = get_matrix_by_name(arg)) != NULL) {
	    ; /* OK, we'll print the matrix */
	} else if ((v = series_index(pdinfo, arg)) < pdinfo->v) {
	    series_v = v; /* OK, we'll print the series */
	} else {
	    x = printf_get_scalar(arg, pZ, pdinfo, t, &err);
	    if (!err && na(x)) {
		fc = fmt[flen - 1] = 's';
		str = gretl_strdup("NA");
	    }
	}
	*pargs += alen;
    }
    free(arg);

    if (err) {
	goto bailout;
    } 

    /* do the actual printing */

    if (m != NULL) {
	/* printing a matrix */
	if (!wstar) wid = -1;
	if (!pstar) prec = -1;
	gretl_matrix_print_with_format(m, fmt, wid, prec, prn);
    } else if (series_v > 0) {
	/* printing a series */
	printf_series(series_v, (const double **)*pZ,
		      pdinfo, fmt, wid, prec, 
		      wstar, pstar, prn);
    } else if (fc == 's') {
	/* printing a string */
	if (wstar && pstar) {
	    pprintf(prn, fmt, wid, prec, str);
	} else if (wstar || pstar) {
	    wid = (wstar)? wid : prec;
	    pprintf(prn, fmt, wid, str);
	} else {
	    pprintf(prn, fmt, str);
	}
    } else if (strchr(intconv, fc)) {
	/* printing a scalar as int */
	if (wstar && pstar) {
	    pprintf(prn, fmt, wid, prec, (int) x);
	} else if (wstar || pstar) {
	    wid = (wstar)? wid : prec;
	    pprintf(prn, fmt, wid, (int) x);
	} else {
	    pprintf(prn, fmt, (int) x);
	}
    } else {
	/* printing a scalar value */
	x *= (1 + 0x1p-52);
	if (wstar && pstar) {
	    pprintf(prn, fmt, wid, prec, x);
	} else if (wstar || pstar) {
	    wid = (wstar)? wid : prec;
	    pprintf(prn, fmt, wid, x);
	} else {
	    pprintf(prn, fmt, x);
	}
    }

 bailout:

    if (err) {
	pputc(prn, '\n');
    }

    free(fmt);
    free(str);
    
    return err;
}

/* split line into format and args, copying both parts */

static int split_printf_line (const char *s, char *targ, int *sp,
			      char **format, char **args)
{
    const char *p;
    int n, err = 0;

    *sp = 0;

    if (!strncmp(s, "printf ", 7)) {
	s += 7;
    } else if (!strncmp(s, "sprintf ", 8)) {
	s += 8;
	*sp = 1;
    }

    if (*sp) {
	/* need a target name */
	s += strspn(s, " ");
	err = extract_varname(targ, s, &n);
	if (!err && n == 0) {
	    err = E_PARSE;
	}
	if (err) {
	    return err;
	}
	s += n;
	/* allow comma after target */
	s += strspn(s, " ");
	if (*s == ',') {
	    s++;
	}
    }

    s += strspn(s, " ");
    if (*s != '"' || *(s+1) == '\0') {
	return E_PARSE;
    }

    s++;
    p = s;

    n = 0;
    while (*s) {
	if (*s == '"' && *(s-1) != '\\') {
	    break;
	}
	n++;
	s++;
    }

    if (n == 0) {
	/* empty format string */
	return 0;
    }

    *format = gretl_strndup(p, n);
    if (*format == NULL) {
	return E_ALLOC;
    }

    s++;
    s += strspn(s, " ");

    if (*s != ',') {
	/* empty args */
	*args = NULL;
	return 0;
    }

    s++;
    s += strspn(s, " ");

    *args = gretl_strdup(s);
    if (*args == NULL) {
	return E_ALLOC;
    }

    return 0;
}

static int real_do_printf (const char *line, double ***pZ, 
			   DATAINFO *pdinfo, PRN *inprn, 
			   int t)
{
    PRN *prn = inprn;
    char targ[VNAMELEN];
    char *format = NULL;
    char *args = NULL;
    int sp, err = 0;

    gretl_error_clear();

    *targ = '\0';

    if (t < 0) {
	t = pdinfo->t1;
    }

    err = split_printf_line(line, targ, &sp, &format, &args);
    if (err) {
	return err;
    }

    if (sp) {
	/* "sprintf" */
	if (*targ == '\0') {
	    return E_PARSE;
	}
	prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
	if (err) {
	    return err;
	}
    } 

#if PSDEBUG
    fprintf(stderr, "do_printf: line = '%s'\n", line);
    fprintf(stderr, " targ = '%s'\n", targ);
    fprintf(stderr, " format = '%s'\n", format);
    fprintf(stderr, " args = '%s'\n", args);
#endif

    if (format != NULL) {
	char *p = format;
	char *q = args;

	while (*p && !err) {
	    if (*p == '%' && *(p+1) == '%') {
		pputc(prn, '%');
		p += 2;
	    } else if (*p == '%') {
		err = print_arg(&p, &q, pZ, pdinfo, t, prn);
	    } else if (*p == '\\') {
		err = printf_escape(*(p+1), prn);
		p += 2;
	    } else {
		pputc(prn, *p);
		p++;
	    }
	}

	if (q != NULL && *q != '\0') {
	    pprintf(prn, "\nunmatched argument '%s'", q);
	    err = E_PARSE;
	}
    }

    if (sp) {
	if (!err) {
	   err = save_named_string(targ, gretl_print_get_buffer(prn),
				   inprn);
	}
	gretl_print_destroy(prn);
    } else if (err) {
	pputc(prn, '\n');
    }

    free(format);
    free(args);

    return err;
}

/**
 * do_printf:
 * @line: command line.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @prn: printing struct.
 *
 * Implements a somewhat limited version of C's printf()
 * for use in gretl scripts.
 *
 * Returns: 0 on success, non-zero on error.
 */

int do_printf (const char *line, double ***pZ, 
	       DATAINFO *pdinfo, PRN *prn)
{
    return real_do_printf(line, pZ, pdinfo, prn, -1);
}

/* below: sscanf apparatus */

static int 
sscanf_target_var (const char *vname, DATAINFO *pdinfo, int *err)
{
    if (gretl_is_scalar(vname)) {
	return 0;
    } else {
	strcpy(gretl_errmsg, _("sscanf: numerical target must be scalar"));
	*err = E_DATA;
	return 0;
    }
}

struct bracket_scan {
    int reverse;
    char *chrs;
};

static int 
parse_bracket_scanner (const char **pfmt, struct bracket_scan *bscan)
{
    const char *p, *s = *pfmt;
    int k = 0, n = 0, closed = 0;

    s++;
    k++;

    if (*s == '^') {
	bscan->reverse = 1;
	s++;
	k++;
    } 
    if (*s == ']') {
	s++;
	k++;
    }
    
    p = s;
    while (*p) {
	if (*p == ']') {
	    closed = 1;
	    break;
	}
	p++;
	n++;
    }

    if (!closed || n == 0) {
	return E_PARSE;
    }

    bscan->chrs = gretl_strndup(s, n);
    if (bscan->chrs == NULL) {
	return E_ALLOC;
    }

#if PSDEBUG
    fprintf(stderr, "bracket_scan: reverse = %d, chars = '%s'\n",
	    bscan->reverse, bscan->chrs);
#endif

    *pfmt += (k + n + 1);

    return 0;
}

#define WIDTH_UNSPEC -1

static char *
get_scanf_format_chunk (const char *s, int *fc, int *len, 
			int *width, int *supp, 
			struct bracket_scan *bscan,
			int *err)
{
    const char *cnvchars = "dfgsm";
    const char *numchars = "0123456789";
    char *chunk = NULL;
    const char *p = s;
    int n;

    p++; /* % */

    if (*p == '*') {
	/* suppress assignment */
	*supp = 1;
	p++;
    } 

    /* optional max width? */
    n = strspn(p, numchars);
    if (n > 0) {
	*width = atoi(p);
	p += n;
    } else {
	*width = WIDTH_UNSPEC;
    }

    /* now we should have a conversion spec */
    if (*p == '[') {
	*err = parse_bracket_scanner(&p, bscan);
	if (!*err) {
	    *fc = 's';
	}
    } else {
	if (*p == 'l' && *(p+1) == 'f') {
	    p++;
	}
	if (*p == '\0' || strchr(cnvchars, *p) == NULL) {
	    fprintf(stderr, "bad conversion '%c'\n", *p);
	    *err = E_PARSE;
	    return NULL;
	}
	*fc = *p++;
    }

    *len = n = p - s;
    
    if (n > 0) {
	chunk = gretl_strndup(s, n);
	if (chunk == NULL) {
	    *err = E_ALLOC;
	}
    } else {
	*err = E_PARSE;
    }

#if PSDEBUG
    fprintf(stderr, "get_scanf_format_chunk: got '%s'\n", chunk);
#endif

    return chunk;
}

static int 
scan_string (const char *targ, char **psrc, int width, 
	     struct bracket_scan *bscan, int *ns)
{
    int n, err = 0;

    *psrc += strspn(*psrc, " \t\n");

    if (bscan->chrs != NULL) {
	if (bscan->reverse) {
	    n = strcspn(*psrc, bscan->chrs);
	} else {
	    n = strspn(*psrc, bscan->chrs);
	}
    } else {	    
	n = strcspn(*psrc, " \t\n");
    }

    if (width != WIDTH_UNSPEC && n > width) {
	n = width;
    }

    if (targ != NULL) {
	if (!is_user_string(targ)) {
	    sprintf(gretl_errmsg, _("%s: not a string variable"), targ);
	    err = E_UNKVAR;
	} else if (n > 0) {
	    char *conv = gretl_strndup(*psrc, n);

	    if (conv == NULL) {
		err = E_ALLOC;
	    } else {
		err = save_named_string(targ, conv, NULL);
		free(conv);
		*psrc += n;
	    }
	    if (!err) {
		*ns += 1;
	    }
	}
    } else {
	*psrc += n;
    }

    return err;
}

static int scan_scalar (const char *targ, char **psrc,
			int fc, int width, double **Z, 
			DATAINFO *pdinfo, int *ns)
{
    char *endp = NULL;
    int v = 0;
    long k = 0;
    double x = 0;
    int err = 0;

    if (targ != NULL) {
	v = sscanf_target_var(targ, pdinfo, &err);
	if (err) {
	    return err;
	}
    }

    errno = 0;

    if (fc == 'd') {
	k = strtol(*psrc, &endp, 10);
    } else {
	x = strtod(*psrc, &endp);
    }

    if (width != WIDTH_UNSPEC) {
	int eaten = endp - *psrc;

	if (eaten > width) {
	    char *tmp = gretl_strndup(*psrc, width);

	    if (tmp == NULL) {
		return E_ALLOC;
	    }
	    if (fc == 'd') {
		k = strtol(tmp, &endp, 10);
	    } else {
		x = strtod(tmp, &endp);
	    }
	    endp = *psrc + width;
	    free(tmp);
	}
    }

    if (!errno && endp != *psrc) {
	if (v == 0) {
	    gretl_scalar_set_value(targ, (fc == 'd')? k : x);
	} else if (v > 0) {
	    Z[v][0] = (fc == 'd')? k : x;
	}
	*ns += 1;
    }

    *psrc = endp;

    return err;
}

static int scan_matrix (const char *targ, char **psrc, int rows, int *ns)
{
    char *src = *psrc;
    char *endp = NULL;
    gretl_matrix *m = NULL;
    double x;
    int r = 0, c = 0, c0 = 0;
    int i, j, err = 0;

    errno = 0;

    while (*src) {
	strtod(src, &endp);
	if (endp == src || errno == ERANGE) {
	    break;
	} 
	c++;
	src = endp;
	src += strspn(src, " \t");
	if (*src == '\r' && *(src+1) == '\n') {
	    src++;
	}
	if (*src == '\n' || *src == '\r' || *src == '\0') {
	    if (r == 0) {
		c0 = c;
	    } else if (c != c0) {
		break;
	    }
	    r++;
	    if (r == rows) {
		break;
	    }
	    c = 0;
	}
    }

    if (targ == NULL) {
	*psrc = src;
	return 0;
    }

    if (get_matrix_by_name(targ) == NULL) {
	return E_UNKVAR;
    }

    if (r > 0 && c0 > 0) {
	m = gretl_matrix_alloc(r, c0);
	if (m == NULL) {
	    return E_ALLOC;
	}	

	src = *psrc;
	for (i=0; i<r; i++) {
	    for (j=0; j<c0; j++) {
		x = strtod(src, &endp);
		gretl_matrix_set(m, i, j, x);
		src = endp;
		src += strspn(src, " \t\n\r");
	    }
	}
 
	err = user_matrix_replace_matrix_by_name(targ, m);
	if (err) {
	    gretl_matrix_free(m);
	} else {
	    *ns += 1;
	}
    }

    *psrc = src;

    return err;
}

static int scan_arg (char **psrc, char **pfmt, char **pargs, 
		     double **Z, DATAINFO *pdinfo,
		     PRN *prn, int *ns)
{
    char *fmt = NULL;
    char *arg = NULL;
    char *str = NULL;
    char *targ = NULL;
    struct bracket_scan bscan;
    int flen = 0, alen = 0;
    int wid = 0, supp = 0;
    int fc = 0;
    int err = 0;

#if PSDEBUG
    fprintf(stderr, "scan_arg: *psrc='%s', *pfmt='%s', *pargs='%s'\n",
	    *psrc, *pfmt, *pargs);
#endif

    bscan.reverse = 0;
    bscan.chrs = NULL;

    /* select current conversion format */
    fmt = get_scanf_format_chunk(*pfmt, &fc, &flen, &wid, &supp, 
				 &bscan, &err);
    if (err) {
	return err;
    }

    *pfmt += flen;

    if (!supp) {
	/* get next substantive arg */
	arg = get_next_arg(*pargs, &alen, &err);
	if (!err) {
	    targ = arg;
	    if (*targ == '&') {
		targ++;
	    }
	    if (*targ == '\0') {
		err = E_PARSE;
	    } else {
		*pargs += alen;
	    }
	}
    }

    /* do the actual scanning */

    if (!err && fc == 's') {
	err = scan_string(targ, psrc, wid, &bscan, ns);
    } else if (!err && fc == 'm') {
	err = scan_matrix(targ, psrc, wid, ns);
    } else if (!err) {
	err = scan_scalar(targ, psrc, fc, wid, Z, pdinfo, ns);
    }

    free(fmt);
    free(arg);
    free(str);
    free(bscan.chrs);
    
    return err;
}

static int get_literal_length (const char *s)
{
    int n = 0;

    while (*s) {
	if (*s == '"' && *(s-1) != '\\') {
	    break;
	}
	n++;
	s++;
    }

    return n;
}

/* split line into source string, format and arguments: 
   acceptable sources are string literals or string variables
*/

static int 
split_scanf_line (const char *s, char **src, char **format, 
		  char **args, int *literal)
{
    const char *p;
    char *srcname;
    int gotcom = 0;
    int n = 0;

    if (!strncmp(s, "sscanf ", 7)) {
	s += 7;
    } 

    s += strspn(s, " ");

    /* string literal or string variable as source */
    if (*s == '"') {
	*literal = 1;
	s++;
	n = get_literal_length(s);
    } else {
	if (*s == '@') s++;
	n = gretl_namechar_spn(s);
    } 

    if (n == 0) {
	return E_PARSE;
    }

    srcname = gretl_strndup(s, n);
    if (srcname == NULL) {
	return E_ALLOC;
    }

    if (*literal) {
	*src = srcname;
	s++;
    } else {
	*src = get_string_by_name(srcname);
	free(srcname);
	if (*src == NULL) {
	    return E_UNKVAR;
	}
    }

    s += n;
    if (*s == ',') {
	gotcom = 1;
	s++;
    }

    /* skip to next field */
    s += strspn(s, " ");
    if (!gotcom && *s == ',') {
	s++;
	s += strspn(s, " ");
    }

    /* get format string */
    if (*s != '"' || *(s+1) == '\0') {
	return E_PARSE;
    }

    s++;
    p = s;

    n = 0;
    while (*s) {
	if (*s == '"' && *(s-1) != '\\') {
	    break;
	}
	n++;
	s++;
    }

    if (n == 0) {
	/* empty format string */
	return 0; 
    }

    *format = gretl_strndup(p, n);
    if (*format == NULL) {
	return E_ALLOC;
    }

    s++;
    s += strspn(s, " ");

    if (*s != ',') {
	/* empty args */
	*args = NULL;
	return 0;
    }

    s++;
    s += strspn(s, " ");

    *args = gretl_strdup(s);
    if (*args == NULL) {
	return E_ALLOC;
    }

    return 0;
}

static int nscan;

int n_scanned_items (void)
{
    return nscan;
}

/**
 * do_sscanf:
 * @line: command line.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @prn: printing struct.
 *
 * Implements a somewhat limited version of C's sscanf()
 * for use in gretl scripts.
 *
 * Returns: 0 on success, non-zero on error.
 */

int do_sscanf (const char *line, double ***pZ, 
	       DATAINFO *pdinfo, PRN *prn)
{
    char *r, *p, *q;
    char *src = NULL;
    char *format = NULL;
    char *args = NULL;
    int literal = 0;
    int err = 0;

    nscan = 0;
    gretl_error_clear();

    err = split_scanf_line(line, &src, &format, &args, &literal);
    if (err) {
	return err;
    }

#if PSDEBUG
    fprintf(stderr, "do_sscanf: line = '%s'\n", line);
    fprintf(stderr, " src = '%s'\n", src);
    fprintf(stderr, " format = '%s'\n", format);
    fprintf(stderr, " args = '%s'\n", args);
#endif

    r = src;
    p = format;
    q = args;

    while (*r && *p && !err) {
	if (isspace(*p)) {
	    while (isspace(*r)) r++;
	    while (isspace(*p)) p++;
	}
	if (*r == '\t' && *p == '\\' && *(p+1) == 't') {
	    r++;
	    p += 2;
	} else if (*r == '\n' && *p == '\\' && *(p+1) == 'n') {
	    r++;
	    p += 2;
	} else if (*r == *p) {
	    r++;
	    p++;
	} else if (*p == '%') {
	    err = scan_arg(&r, &p, &q, *pZ, pdinfo, prn, &nscan);
	} else {
	    break;
	}
    }

    if (literal) {
	free(src);
    }

    free(format);
    free(args);

    if (!err && gretl_messages_on()) {
	pprintf(prn, "Number of items successfully scanned = %d\n", 
		nscan);
    }

    return err;
}

/* The originating command is of form:

     genr markers = format, arg1,...

   We're assuming that we're just getting the part after
   the equals sign here, in the variable s.
*/

int generate_obs_markers (const char *s, double ***pZ, DATAINFO *pdinfo)
{
    PRN *prn;
    int t, err = 0;

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);

    if (err) {
	return err;
    }

    if (pdinfo->S == NULL) {
	err = dataset_allocate_obs_markers(pdinfo);
    }

    if (!err) {
	const char *buf;

	for (t=0; t<pdinfo->n && !err; t++) {
	    gretl_print_reset_buffer(prn);
	    err = real_do_printf(s, pZ, pdinfo, prn, t);
	    if (!err) {
		buf = gretl_print_get_buffer(prn);
		pdinfo->S[t][0] = '\0';
		strncat(pdinfo->S[t], buf, OBSLEN - 1);
	    }
	}
    }

    gretl_print_destroy(prn);
	
    return err;
}
