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
#include "gretl_func.h"
#include "gretl_string_table.h"

#include <errno.h>

#define PSDEBUG 0

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

/* various string argument variants, optionally followed
   by "+ offset" */

static char *printf_get_string (const char *s, const double **Z,
				const DATAINFO *pdinfo, 
				int t, int *err)
{
    char *ret = NULL;
    char tmpstr[16], strarg[16];
    const char *p = NULL, *q = NULL;
    int v, offset = 0;
    int len = 0;

#if PSDEBUG
    fprintf(stderr, "printf_get_string: looking at '%s'\n", s);
#endif

    if (*s == '@') {
	char sname[VNAMELEN] = {0};
	int n = gretl_varchar_spn(s + 1);

	len = n;
	if (len >= VNAMELEN) {
	    len = VNAMELEN - 1;
	}
	strncat(sname, s + 1, len);
	p = get_named_string(sname);
	if (p != NULL) {
	    len = strlen(p);
	    q = s + 1 + n;
	}
    } else if (*s == '"') {
	/* literal string */
	p = strrchr(s + 1, '"');
	if (p != NULL) {
	    q = p + 1;
	    len = p - s - 1;
	    p = s + 1;
	}
    } else if (sscanf(s, "varname(%d)", &v)) {
	/* name of variable identified by number */
	if (v >= 0 && v < pdinfo->v) {
	    p = pdinfo->varname[v];
	    len = strlen(p);
	    q = strchr(s, ')') + 1;
	}
    } else if (sscanf(s, "argname(%15[^)])", strarg)) {
	/* name of function argument */
	char *aname = gretl_func_get_arg_name(strarg);

	if (aname != NULL) {
	    strcpy(tmpstr, aname);
	    p = tmpstr;
	    len = strlen(p);
	    q = strchr(s, ')') + 1;
	    free(aname);
	}
    } else if (sscanf(s, "date(%15[^)])", strarg)) {
	/* date string */
	t = -1;
	if (isdigit(*strarg)) {
	    t = atoi(strarg);
	} else {
	    v = varindex(pdinfo, strarg);
	    if (v < pdinfo->v) {
		t = Z[v][0];
	    }
	}
	if (t > 0 && t <= pdinfo->n) {
	    ntodate(tmpstr, t - 1, pdinfo);
	    p = tmpstr;
	    len = strlen(p);
	    q = strchr(s, ')') + 1;
	}
    } else if (!strncmp(s, "marker", 6) && pdinfo->S != NULL) {
	/* observation label */
	p = pdinfo->S[t];
	len = strlen(p);
	q = s + 6;
    }

    if (p != NULL) {
	if (q != NULL) {
	    while (isspace(*q)) q++;
	    if (*q == '+') {
		q++;
		offset = atoi(q);
		len -= offset;
	    }
	}
	if (len >= 0) {
	    ret = gretl_strndup(p + offset, len);
	}
    }

    if (ret == NULL) {
	ret = gretl_strdup("NA");
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
	return atof(s);
    }

    v = varindex(pdinfo, s);

    if (v < pdinfo->v && var_is_series(pdinfo, v)) {
	char genstr[32];

	sprintf(genstr, "%s[%d]", s, t + 1);
	x = generate_scalar(genstr, pZ, pdinfo, err);
    } else {
	x = generate_scalar(s, pZ, pdinfo, err);
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
    const char *cnvchars = "eEfgGxduxs";
    const char *numchars = "0123456789";
    char *chunk = NULL;
    const char *p = s;
    int n;

    p++; /* % */

    /* '#', if present, must come first */
    if (*p == '#') {
	p++;
    }

    /* 'width' could be < 0 */
    if (*p == '-') {
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
	if (fc == 's') {
	    str = printf_get_string(arg, (const double **) *pZ, 
				    pdinfo, t, &err);
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

    if (fc == 's') {
	if (wstar && pstar) {
	    pprintf(prn, fmt, wid, prec, str);
	} else if (wstar || pstar) {
	    wid = (wstar)? wid : prec;
	    pprintf(prn, fmt, wid, str);
	} else {
	    pprintf(prn, fmt, str);
	}
    } else if (strchr(intconv, fc)) {
	if (wstar && pstar) {
	    pprintf(prn, fmt, wid, prec, (int) x);
	} else if (wstar || pstar) {
	    wid = (wstar)? wid : prec;
	    pprintf(prn, fmt, wid, (int) x);
	} else {
	    pprintf(prn, fmt, (int) x);
	}
    } else {
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
    int n;

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
	n = gretl_varchar_spn(s);
	if (n == 0 || n >= VNAMELEN) {
	    return E_PARSE;
	} else {
	    *targ = '\0';
	    strncat(targ, s, n);
	    s += n;
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
    char *p, *q;
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
	prn = gretl_print_new(GRETL_PRINT_BUFFER);
	if (prn == NULL) {
	    return E_ALLOC;
	}
    } 

#if PSDEBUG
    fprintf(stderr, "do_printf: line = '%s'\n", line);
    fprintf(stderr, " targ = '%s'\n", targ);
    fprintf(stderr, " format = '%s'\n", format);
    fprintf(stderr, " args = '%s'\n", args);
#endif

    p = format;
    q = args;

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
 * Implement a somewhat limited version of C's printf
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
    int v = varindex(pdinfo, vname);

    if (v <= 0 || v >= pdinfo->v) {
	*err = E_UNKVAR;
    } else if (!var_is_scalar(pdinfo, v)) {
	strcpy(gretl_errmsg, _("sscanf: numerical target must be scalar"));
	*err = E_DATA;
    }

    return v;
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
    fprintf(stderr, "get_printf_format_chunk: got '%s'\n", chunk);
#endif

    return chunk;
}

static int 
scan_string (const char *targ, const char **psrc, int width, 
	     struct bracket_scan *bscan)
{
    int n, err = 0;

    if (bscan->chrs != NULL) {
	if (bscan->reverse) {
	    n = strcspn(*psrc, bscan->chrs);
	} else {
	    n = strspn(*psrc, bscan->chrs);
	}
    } else {	    
	*psrc += strspn(*psrc, " \t\n");
	n = strcspn(*psrc, " \t\n");
    }

    if (width != WIDTH_UNSPEC && n > width) {
	n = width;
    }

    if (targ != NULL) {
	if (!is_user_string(targ)) {
	    sprintf(gretl_errmsg, _("%s: not a string variable"), targ);
	    err = E_UNKVAR;
	} else {
	    char *conv = gretl_strndup(*psrc, n);

	    if (conv == NULL) {
		err = E_ALLOC;
	    } else {
		err = save_named_string(targ, conv, NULL);
		free(conv);
		*psrc += n;
	    }
	}
    } else {
	*psrc += n;
    }

    return err;
}

static int scan_scalar (const char *targ, const char **psrc,
			int fc, double **Z, DATAINFO *pdinfo)
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

    if (errno == ERANGE) {
	sprintf(gretl_errmsg, _("'%s' -- number out of range!"), *psrc);
	err = E_DATA;
    } else if (endp == *psrc) {
	sprintf(gretl_errmsg, _("'%s' -- no numeric conversion performed!"), *psrc);
	err = E_PARSE;
    } else {
	if (v > 0) {
	    Z[v][0] = (fc == 'd')? k : x;
	}
	*psrc = endp;
    }

    return err;
}

static int scan_matrix (const char *targ, const char **psrc)
{
    const char *src = *psrc;
    char *endp = NULL;
    gretl_matrix *m = NULL;
    double x;
    int cmax = 0, cmin = 999999;
    int r = 0, c = 0;
    int i, j, err = 0;

    errno = 0;

    while (*src && !err) {
	strtod(src, &endp);
	if (endp == src) {
	    break;
	} else if (errno == ERANGE) {
	    err = E_DATA;
	    break;
	}
	c++;
	src = endp;
	src += strspn(src, " \t");
	if (*src == '\r' && *(src+1) == '\n') {
	    src++;
	}
	if (*src == '\n' || *src == '\r' || *src == '\0') {
	    r++;
	    if (c > cmax) cmax = c;
	    if (c < cmin) cmin = c;
	    c = 0;
	}
    }

    if (!err && (r == 0 || cmax == 0 || cmin != cmax)) {
	err = E_DATA;
    }

    if (err) return err;

    if (targ == NULL) {
	*psrc = src;
	return 0;
    }

    if (get_matrix_by_name(targ) == NULL) {
	return E_UNKVAR;
    }

    m = gretl_matrix_alloc(r, cmax);
    if (m == NULL) {
	return E_ALLOC;
    }

    src = *psrc;

    for (i=0; i<r; i++) {
	for (j=0; j<cmax; j++) {
	    x = strtod(src, &endp);
	    gretl_matrix_set(m, i, j, x);
	    src = endp;
	    src += strspn(src, " \t\n\r");
	}
    }

    err = user_matrix_replace_matrix_by_name(targ, m);
    if (err) {
	gretl_matrix_free(m);
    }

    *psrc = src;

    return err;
}

static int scan_arg (const char **psrc, char **pfmt, char **pargs, 
		     double **Z, DATAINFO *pdinfo,
		     PRN *prn)
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
	err = scan_string(targ, psrc, wid, &bscan);
    } else if (!err && fc == 'm') {
	err = scan_matrix(targ, psrc);
    } else if (!err) {
	err = scan_scalar(targ, psrc, fc, Z, pdinfo);
    }

    free(fmt);
    free(arg);
    free(str);
    free(bscan.chrs);
    
    return err;
}

/* split line into source string, format and arguments: for now
   we only accept named string variables as sources 
*/

static int 
split_scanf_line (const char *s, const char **src, char **format, char **args)
{
    const char *p;
    char *srcname;
    int gotcom = 0;
    int n;

    if (!strncmp(s, "sscanf ", 7)) {
	s += 7;
    } 

    s += strspn(s, " ");
    if (*s != '@') {
	return E_PARSE;
    }

    s++;
    n = gretl_varchar_spn(s);
    if (n == 0) {
	return E_PARSE;
    }

    srcname = gretl_strndup(s, n);
    if (srcname == NULL) {
	return E_ALLOC;
    } else {
	*src = get_named_string(srcname);
	free(srcname);
	if (*src == NULL) {
	    return E_UNKVAR;
	}
	s += n;
	if (*s == ',') {
	    gotcom = 1;
	    s++;
	}
    }

    /* skip to next field */
    s += strspn(s, " ");
    if (!gotcom && *s == ',') {
	s++;
	s += strspn(s, " ");
    }

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

/**
 * do_sscanf:
 * @line: command line.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @prn: printing struct.
 *
 * Implement a somewhat limited version of C's sscanf
 * for use in gretl scripts.
 *
 * Returns: 0 on success, non-zero on error.
 */

int do_sscanf (const char *line, double ***pZ, 
	       DATAINFO *pdinfo, PRN *prn)
{
    const char *r, *src = NULL;
    char *p, *q;
    char *format = NULL;
    char *args = NULL;
    int err = 0;

    gretl_error_clear();

    err = split_scanf_line(line, &src, &format, &args);
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
	while (isspace(*r)) r++;
	while (isspace(*p)) p++;
	if (*r == *p) {
	    r++;
	    p++;
	} else if (*p == '%') {
	    err = scan_arg(&r, &p, &q, *pZ, pdinfo, prn);
	} else {
	    sprintf(gretl_errmsg, _("Unmatched element in format: '%s'"), p);
	    err = E_DATA;
	}
    }

    if (!err && *p && !*r) {
	sprintf(gretl_errmsg, _("Unmatched element in format: '%s'"), p);
	err = E_DATA;
    }

    free(format);
    free(args);

    return err;
}

/* The originating command is of form:

     genr markers = format, arg1,...

   We're assuming that we're just getting the part after
   the equals sign here, in the variable s.
*/

int generate_obs_markers (const char *s, double ***pZ, DATAINFO *pdinfo)
{
    PRN *prn = gretl_print_new(GRETL_PRINT_BUFFER);
    int t, err = 0;

    if (prn == NULL) {
	return E_ALLOC;
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
