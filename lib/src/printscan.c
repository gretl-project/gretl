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
#include "matrix_extra.h"
#include "uservar.h"

#include <errno.h>

#define PSDEBUG 0

/* For the point of "x *= (1 + 0x1p-52)" below, see
   http://sourceware.org/bugzilla/show_bug.cgi?id=4943 ,
   contribution from Eric Postpischil, 2007-10-05.
*/

static void printf_series (int v, const DATASET *dset,
			   const char *fmt, 
			   int wid, int prec, 
			   int wstar, int pstar,
			   PRN *prn)
{
    int obslen = max_obs_marker_length(dset);
    double x;
    int t;

    for (t=dset->t1; t<=dset->t2; t++) {
	print_obs_marker(t, dset, obslen, prn);
	x = dset->Z[v][t];
	if (!na(x)) {
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

static char *printf_get_string (char *s, DATASET *dset, 
				int t, int *err)
{
    char *ret = NULL;

    if (!strncmp(s, "marker", 6) && dset != NULL && dset->S != NULL) {
	/* special, not in regular genr */
	const char *p = dset->S[t];
	const char *q = s + 6;
	int len = strlen(p);
	int offset = 0;

	while (isspace(*q)) q++;

	if (*q == '+') {
	    offset = atoi(++q);
	    len -= offset;
	}
	if (len >= 0) {
	    ret = gretl_strndup(p + offset, len);
	    if (ret == NULL) {
		*err = E_ALLOC;
	    }
	}
    } else {
	ret = generate_string(s, dset, err);
    }

    if (ret == NULL && !*err) {
	*err = E_DATA;
    }

    return ret;
}

/* here we're trying to get a field width or precision
   specifier */

static int printf_get_int (const char *s, DATASET *dset, 
			   int *err)
{
    double x;
    int ret = 0;

#if PSDEBUG
    fprintf(stderr, "printf_get_int: looking at '%s'\n", s);
#endif

    if (numeric_string(s)) {
	x = dot_atof(s);
    } else if (gretl_is_scalar(s)) {
	x = gretl_scalar_get_value(s, NULL);
    } else {
	x = generate_scalar(s, dset, err);
    }

    if (!*err && (xna(x) || fabs(x) > 255)) {
	*err = E_DATA;
    } 

    if (!*err) {
	ret = (int) x;
    }

#if PSDEBUG
    fprintf(stderr, "printf_get_int: returning %g\n", x);
#endif

    return ret;
}

/* gen_arg_val: when this function is called we have already
   determined that the argument string @s is not associated with a
   string format specifier, is not a literal numerical value, and is
   not the name of an existing scalar, matrix or series.  So now we
   attempt to generate a value (of type unknown).
*/

static int gen_arg_val (const char *s, DATASET *dset, 
			double *px, gretl_matrix **pm, 
			int *pv, int *scalar)
{
    const char *name = "$ptmp";
    char *genstr;
    int err = 0;

    genstr = gretl_strdup_printf("%s=%s", name, s); 
    if (genstr == NULL) {
	return E_ALLOC;
    }

    err = generate(genstr, dset, OPT_P, NULL);

    if (!err) {
	int type = genr_get_last_output_type();

	if (type == GRETL_TYPE_DOUBLE) {
	    *scalar = 1;
	    *px = gretl_scalar_get_value(name, NULL);
	    user_var_delete_by_name(name, NULL);
	} else if (type == GRETL_TYPE_MATRIX) {
	    gretl_matrix *m = get_matrix_by_name(name);

	    *pm = gretl_matrix_copy(m);
	    if (*pm == NULL) {
		err = E_ALLOC;
	    }
	    user_var_delete_by_name(name, NULL);
	} else if (type == GRETL_TYPE_SERIES) {
	    *pv = current_series_index(dset, name);
	} else {
	    err = E_TYPES; /* FIXME cleanup? */
	}
    }

    free(genstr);

    return err;
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
	*err = E_ARGS;
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

/* Extract the next conversion spec from *pfmt, find and evaluate the
   corresponding elements in *pargs, and print the result. Note
   that @t should be negative (and therefore ignored) unless
   we're doing the "genr markers" special thing.
*/

static int print_arg (char **pfmt, char **pargs, DATASET *dset,
		      int t, PRN *prn)
{
    const char *intconv = "dxul";
    char *fmt = NULL;
    char *arg = NULL;
    char *str = NULL;
    gretl_matrix *m = NULL;
    int free_m = 0;
    int series_v = -1;
    int free_v = 0;
    double x = NADBL;
    int flen = 0, alen = 0;
    int wstar = 0, pstar = 0;
    int wid = 0, prec = 0;
    int got_scalar = 0;
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
	    x = printf_get_int(arg, dset, &err);
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
	    x = printf_get_int(arg, dset, &err);
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
	int v = -1;

	if (fc == 's') {
	    /* must be printing a string */
	    str = printf_get_string(arg, dset, t, &err);
	} else if (numeric_string(arg)) {
	    /* printing a scalar value */
	    x = dot_atof(arg);
	    got_scalar = 1;
	} else if (gretl_is_scalar(arg)) {
	    /* printing a named scalar */
	    x = gretl_scalar_get_value(arg, NULL);
	    got_scalar = 1;
	} else if ((m = get_matrix_by_name(arg)) != NULL) {
	    ; /* printing a named matrix */
	} else if ((v = current_series_index(dset, arg)) >= 0) {
	    /* printing a named series, unless t >= 0 */
	    if (t < 0) {
		series_v = v; 
	    } else {
		x = dset->Z[v][t];
		got_scalar = 1;
	    }
	} else {
	    /* try treating as generated value */
	    err = gen_arg_val(arg, dset, &x, &m, &v, &got_scalar);
	    if (m != NULL) {
		free_m = 1;
	    } else if (v > 0) {
		series_v = v;
		free_v = 1;
	    }
	} 
	*pargs += alen;
    }

    free(arg);

    if (err) {
	goto bailout;
    } 

    if (got_scalar && na(x)) {
	fc = fmt[flen - 1] = 's';
	str = gretl_strdup("NA");
    }	

    /* do the actual printing */

    if (m != NULL) {
	/* printing a matrix */
	if (!wstar) wid = -1;
	if (!pstar) prec = -1;
	gretl_matrix_print_with_format(m, fmt, wid, prec, prn);
    } else if (series_v >= 0) {
	/* printing a series */
	printf_series(series_v, dset, fmt, wid, prec, 
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
    } else if (got_scalar) {
	/* printing a (non-missing) scalar value */
	if (!isnan(x) && !isinf(x)) {
	    x *= (1 + 0x1p-52);
	}
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
    
    if (free_m) {
	gretl_matrix_free(m);
    } else if (free_v) {
	dataset_drop_last_variables(dset, 1);
    }
    
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

#define SPRINTF_COMPAT 1

static user_var *get_stringvar_target (const char *targ, 
				       DATASET *dset,
				       PRN *prn,
				       int *err)
{
    user_var *uvar = NULL;

    if (*targ == '\0') {
	*err = E_PARSE;
    } else {
	uvar = get_user_var_by_name(targ);
	if (uvar != NULL && 
	    user_var_get_type(uvar) != GRETL_TYPE_STRING) {
	    *err = E_TYPES;
	    return NULL;
	}
#if SPRINTF_COMPAT
	if (uvar == NULL) {
	    char genline[64];

	    sprintf(genline, "string %s = \"\"", targ);
	    *err = generate(genline, dset, OPT_NONE, prn);
	    if (!*err) {
		uvar = get_user_var_of_type_by_name(targ, GRETL_TYPE_STRING);
		
	    }
	}
#endif	
	if (uvar == NULL) {
	    gretl_errmsg_sprintf(_("%s: not a string variable"), targ);
	    *err = E_DATA;
	}
    }

    return uvar;
}

/* supports both printf and sprintf */

static int real_do_printf (const char *line, DATASET *dset, 
			   PRN *inprn, int t)
{
    PRN *prn = NULL;
    char targ[VNAMELEN];
    char *format = NULL;
    char *args = NULL;
    user_var *uvar = NULL;
    int sp, err = 0;

    gretl_error_clear();

    *targ = '\0';

    err = split_printf_line(line, targ, &sp, &format, &args);
    if (err) {
	return err;
    }

    if (sp) {
	/* sprintf: we need a target string variable */
	uvar = get_stringvar_target(targ, dset, inprn, &err);
	if (err) {
	    return err;
	}
    }

    /* Even for printf we'll buffer the output locally in 
       case there's an error part way through the printing.
    */
    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
    if (err) {
	return err;
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
		err = print_arg(&p, &q, dset, t, prn);
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

    if (!err) {
	const char *buf = gretl_print_get_buffer(prn);

	if (sp) {
	    char *tmp = gretl_strdup(buf);

	    if (tmp == NULL) {
		err = E_ALLOC;
	    } else {
		user_var_replace_value(uvar, tmp);
	    }
	} else {
	    pputs(inprn, buf);
	}
    } else if (!sp) {
	pputc(inprn, '\n');
    }

    gretl_print_destroy(prn);
    free(format);
    free(args);

    return err;
}

/**
 * do_printf:
 * @line: command line.
 * @dset: dataset struct.
 * @prn: printing struct.
 *
 * Implements a somewhat limited version of C's printf()
 * for use in gretl scripts.
 *
 * Returns: 0 on success, non-zero on error.
 */

int do_printf (const char *line, DATASET *dset, PRN *prn)
{
    return real_do_printf(line, dset, prn, -1);
}

/* below: sscanf apparatus */

struct bracket_scan {
    int reverse;
    char *chrs;
};

static int bscan_escape (char *s)
{
    int i, n = strlen(s);

    for (i=0; i<n; i++) {
	if (s[i] == '\\' && (i == 0 || s[i-1] != '\\')) {
	    if (s[i+1] == 'n') {
		s[i] = '\n';
		shift_string_left(s + i + 1, 1);
		i++;
	    } else if (s[i+1] == 't') {
		s[i] = '\t';
		shift_string_left(s + i + 1, 1);
		i++;
	    }		
	}
    }

    return 0;
}

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

    bscan_escape(bscan->chrs);

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

#if PSDEBUG
    fprintf(stderr, "scan_string: *psrc='%s', n=%d\n", *psrc, n);
    int m = strcspn(*psrc, "\n");
    fprintf(stderr, "m=%d (chrs='%s')\n", m, bscan->chrs);
#endif

    if (width != WIDTH_UNSPEC && n > width) {
	n = width;
    }

    if (targ != NULL) {
	user_var *uvar;

	uvar = get_user_var_of_type_by_name(targ, GRETL_TYPE_STRING);

	if (uvar == NULL) {
	    gretl_errmsg_sprintf(_("%s: not a string variable"), targ);
	    err = E_DATA;
	} else if (n > 0) {
	    char *conv = gretl_strndup(*psrc, n);

	    if (conv == NULL) {
		err = E_ALLOC;
	    } else {
		user_var_replace_value(uvar, conv);
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

static int scan_scalar (char *targ, char **psrc,
			int fc, int width, DATASET *dset, 
			int *ns)
{
    char *endp = NULL;
    long k = 0;
    double x = 0;
    int err = 0;

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
	char *genline;

	if (fc == 'd') {
	    genline = gretl_strdup_printf("%s=%ld", targ, k);
	} else {
	    genline = gretl_strdup_printf("%s=%.16g", targ, x);
	}
	err = generate(genline, dset, OPT_P, NULL);
	free(genline);
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
		     DATASET *dset, int *ns)
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
#if PSDEBUG
	fprintf(stderr, "calling scan_string\n");
#endif
	err = scan_string(targ, psrc, wid, &bscan, ns);
    } else if (!err && fc == 'm') {
	err = scan_matrix(targ, psrc, wid, ns);
    } else if (!err) {
	err = scan_scalar(targ, psrc, fc, wid, dset, ns);
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

static char *get_literal_or_stringvar (const char **src, int *err)
{
    const char *s = *src;
    char *ret = NULL;
    char *tmp = NULL;
    int literal = 0;
    int svar = 0;
    int n = 0;

    if (*s == '"') {
	literal = 1;
	s++;
	n = get_literal_length(s);
    } else if (*s == '@') {
	svar = 1;
	s++;
	n = gretl_namechar_spn(s);
    } else {
	n = gretl_charpos(',', s);
	if (n == gretl_namechar_spn(s)) {
	    svar = 1;
	}
    }

    if (n == 0) {
	/* nothing there */
	*err = E_PARSE;
	return NULL;
    }

    tmp = gretl_strndup(s, n);
    if (tmp == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    if (literal) {
	ret = tmp;
	s++;
    } else if (svar) {
	const char *p = get_string_by_name(tmp);

	if (p == NULL) {
	    *err = E_UNKVAR;
	} else {
	    ret = gretl_strdup(p);
	    if (ret == NULL) {
		*err = E_ALLOC;
	    }
	}
	free(tmp);
    } else {
	ret = generate_string(tmp, NULL, err);
    }

    /* advance the caller's pointer */
    *src = s + n;

    return ret;
}

/* Eat white space; detect comma and eat if found; eat more white
   space.  Return 1 if comma was found else 0.
*/

int eat_comma (const char **src) 
{
    const char *s = *src;
    int ret = 0;

    s += strspn(s, " ");

    if (*s == ',') {
	ret = 1;
	s++;
	s += strspn(s, " ");
    }

    *src = s;

    return ret;
}

/* split line into source string, format and arguments: 
   acceptable sources are string literals or string variables
*/

static int 
split_scanf_line (const char *s, char **psrc, char **pfmt, 
		  char **args)
{
    int gotcom, err = 0;

    if (!strncmp(s, "sscanf ", 7)) {
	s += 7;
    } 

    s += strspn(s, " ");

    /* get the source string */
    *psrc = get_literal_or_stringvar(&s, &err);
    if (err) {
	return err;
    }
    
    /* skip to format string */
    gotcom = eat_comma(&s);
    if (!gotcom) {
	/* empty format */
	return E_PARSE;
    }

    /* get format string */
    *pfmt = get_literal_or_stringvar(&s, &err);
    if (err) {
	return err;
    }

    /* skip to first arg */
    gotcom = eat_comma(&s);
    if (!gotcom) {
	/* empty args */
	return E_PARSE;
    }

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
 * @dset: dataset struct.
 * @prn: printing struct.
 *
 * Implements a somewhat limited version of C's sscanf()
 * for use in gretl scripts.
 *
 * Returns: 0 on success, non-zero on error.
 */

int do_sscanf (const char *line, DATASET *dset, PRN *prn)
{
    char *r, *p, *q;
    char *src = NULL;
    char *format = NULL;
    char *args = NULL;
    int err = 0;

    nscan = 0;
    gretl_error_clear();

#if PSDEBUG
    fprintf(stderr, "do_sscanf: line = '%s'\n", line);
#endif

    err = split_scanf_line(line, &src, &format, &args);
    if (err) {
	return err;
    }

#if PSDEBUG
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
	    err = scan_arg(&r, &p, &q, dset, &nscan);
	} else {
	    break;
	}
    }

    free(src);
    free(format);
    free(args);

    if (!err && gretl_messages_on() && !gretl_looping_quietly()) {
	pprintf(prn, "Number of items successfully scanned = %d\n", 
		nscan);
    }

    return err;
}

/* The originating command is of form:

     genr markers = format, arg1,...

   We're assuming that we're just getting the part after
   the equals sign here, in the variable @s.
*/

int generate_obs_markers (const char *s, DATASET *dset)
{
    PRN *prn;
    int t, err = 0;

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);

    if (err) {
	return err;
    }

    if (dset->S == NULL) {
	err = dataset_allocate_obs_markers(dset);
    }

    if (!err) {
	const char *buf;

	for (t=0; t<dset->n && !err; t++) {
	    gretl_print_reset_buffer(prn);
	    err = real_do_printf(s, dset, prn, t);
	    if (!err) {
		buf = gretl_print_get_buffer(prn);
		dset->S[t][0] = '\0';
		strncat(dset->S[t], buf, OBSLEN - 1);
	    }
	}
    }

    gretl_print_destroy(prn);
	
    return err;
}
